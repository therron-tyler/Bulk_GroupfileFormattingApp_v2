# Tyler Therron
# Winter Lab Macrophage Genomics - Rheumatology Department - Feinberg School of Medicine
# Group File Maker App for Winter Lab Database of Bulk RNA-sequencing Database (v2)
# 2025/11/05

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(colourpicker)
  library(jsonlite)
})

options(shiny.maxRequestSize = 1000 * 1024^2)

`%||%` <- function(x, y) {
  if (is.null(x) || (is.character(x) && length(x) == 0)) y else x
}

# ----------------------- Helpers ----------------------- #
fallback_palette <- function(n){
  hues <- seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

infer_id_cols <- function(df){
  cand <- names(df)
  patt <- c(
    "^symbol$", "^gene_?symbol$", "^genesymbol$",
    "^ensembl(_?gene)?(_?id)?$", "^ensembl(\\.id)?$"
  )
  cand[sapply(cand, function(nm){
    any(grepl(paste(patt, collapse = "|"), nm, ignore.case = TRUE))
  })]
}

numeric_sample_cols <- function(df, id_cols){
  nums <- names(df)[vapply(df, is.numeric, logical(1))]
  setdiff(nums, id_cols)
}

# ----------------------- UI ----------------------- #
ui <- fluidPage(
  titlePanel("Winter Lab Bulk RNA-seq: Groupfile Builder (v2)"),
  tags$hr(),
  
  sidebarLayout(
    sidebarPanel(width = 4,
                 h4("1) Upload dataset"),
                 fileInput("expression_dataframe", "Choose CSV file",
                           accept = c(".csv","text/csv","text/plain")),
                 helpText("Any of these header styles are fine:",
                          tags$ul(
                            tags$li("Both: 'EnsemblID' + 'Gene_Symbol' (or 'Symbol', 'GeneSymbol')"),
                            tags$li("Only 'EnsemblID'"),
                            tags$li("Only 'Gene_Symbol' / 'GeneSymbol' / 'Symbol'"),
                            tags$li("Neither (we'll still detect numeric sample columns)")
                          )
                 ),
                 uiOutput("id_col_detect"),
                 tags$hr(),
                 
                 h4("2) Choose sample columns"),
                 uiOutput("sample_picker"),
                 actionButton("toggle_all", "Select/Deselect All"),
                 tags$hr(),
                 
                 h4("3) Define groups & colors"),
                 numericInput("n_groups", "Number of groups", value = 2, min = 1, step = 1),
                 uiOutput("group_name_color_ui"),
                 actionButton("auto_colors", "Auto-assign colors"),
                 tags$hr(),
                 
                 h4("4) Grouping method"),
                 radioButtons(
                   "group_mode", NULL,
                   choices = c("Token-based (auto)" = "token",
                               "Manual table"       = "manual"),
                   selected = "token"
                 ),
                 
                 conditionalPanel(
                   condition = "input.group_mode == 'token'",
                   uiOutput("token_help"),
                   textInput("token_delim", "Delimiter regex", value = "[-_.]+"),
                   uiOutput("token_index_ui"),
                   checkboxInput("token_drop_singletons", "Hide tokens with only 1 sample", value = FALSE),
                   actionButton("apply_token", "Apply token grouping")
                 ),
                 
                 tags$hr(),
                 downloadButton("download_groupfile", "Download groupfile (Sample,Group,Color)")
    ),
    
    mainPanel(width = 8,
              tabsetPanel(id = "tabs",
                          tabPanel("Preview data",
                                   br(),
                                   uiOutput("template_hint"),
                                   DTOutput("head_table")
                          ),
                          tabPanel("Assignment table",
                                   br(),
                                   helpText("Assign groups with token-based auto-fill, then fine-tune here. Multi-select rows and apply a group in bulk."),
                                   fluidRow(
                                     column(6, selectInput("bulk_group", "Set selected rows to group", choices = NULL)),
                                     column(6, actionButton("apply_bulk", "Apply to selected"))
                                   ),
                                   DTOutput("assign_dt")
                          )
              )
    )
  )
)

# ----------------------- Server ----------------------- #
server <- function(input, output, session){
  
  # ---- Data load ----
  raw_df <- reactive({
    req(input$expression_dataframe)
    readr::read_csv(input$expression_dataframe$datapath, show_col_types = FALSE)
  })
  
  output$head_table <- renderDT({
    req(raw_df())
    datatable(head(raw_df(), 10),
              options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # ---- ID columns ----
  output$id_col_detect <- renderUI({
    req(raw_df())
    df  <- raw_df()
    ids <- infer_id_cols(df)
    if (!length(ids)) {
      div(style = "margin-top:-10px;",
          helpText("No explicit gene ID columns detected. That's ok."))
    } else {
      checkboxGroupInput(
        "id_cols",
        "Detected gene ID columns (optional):",
        choices = ids,
        selected = ids
      )
    }
  })
  
  # ---- Sample columns ----
  output$sample_picker <- renderUI({
    req(raw_df())
    id_cols  <- input$id_cols %||% infer_id_cols(raw_df())
    num_cols <- numeric_sample_cols(raw_df(), id_cols)
    choices  <- if (length(num_cols)) num_cols else setdiff(names(raw_df()), id_cols)
    
    selectizeInput(
      "sample_cols", "Sample columns",
      choices = choices,
      selected = choices,
      multiple = TRUE,
      options = list(plugins = list("remove_button"))
    )
  })
  
  observeEvent(input$toggle_all, {
    req(raw_df())
    id_cols  <- input$id_cols %||% infer_id_cols(raw_df())
    num_cols <- numeric_sample_cols(raw_df(), id_cols)
    offered  <- if (length(num_cols)) num_cols else setdiff(names(raw_df()), id_cols)
    
    curr <- input$sample_cols %||% character(0)
    if (length(curr) < length(offered)) {
      updateSelectizeInput(session, "sample_cols", selected = offered, server = TRUE)
    } else {
      updateSelectizeInput(session, "sample_cols", selected = character(0), server = TRUE)
    }
  })
  
  # ---- Groups: names & colors ----
  groups_rv <- reactiveVal(
    tibble(
      Group = paste0("Group ", 1:2),
      Color = c("#E86FDE", "#3CC9B9")
    )
  )
  
  observeEvent(input$n_groups, {
    ng  <- input$n_groups %||% 1
    cur <- groups_rv()
    if (nrow(cur) < ng) {
      add <- tibble(
        Group = paste0("Group ", (nrow(cur) + 1):ng),
        Color = fallback_palette(ng - nrow(cur))
      )
      groups_rv(bind_rows(cur, add))
    } else if (nrow(cur) > ng) {
      groups_rv(cur[seq_len(ng), ])
    }
  }, ignoreInit = TRUE)
  
  output$group_name_color_ui <- renderUI({
    df <- groups_rv()
    req(nrow(df) > 0)
    tagList(
      lapply(seq_len(nrow(df)), function(i){
        fluidRow(
          column(
            8,
            textInput(
              paste0("grp_name_", i),
              label = paste("Group", i, "name"),
              value = df$Group[i]
            )
          ),
          column(
            4,
            colourInput(
              paste0("grp_col_", i),
              label = "Color",
              value = df$Color[i],
              showColour = "both"
            )
          )
        )
      })
    )
  })
  
  # sync group table from inputs (NULL-safe)
  observe({
    req(input$n_groups)
    ng  <- input$n_groups
    cur <- groups_rv()
    fallback_cols <- fallback_palette(max(ng, 1))
    
    Group <- character(ng)
    Color <- character(ng)
    
    for (i in seq_len(ng)) {
      g_in <- input[[paste0("grp_name_", i)]]
      c_in <- input[[paste0("grp_col_", i)]]
      
      Group[i] <- if (!is.null(g_in) && nzchar(g_in)) {
        g_in
      } else if (i <= nrow(cur)) {
        cur$Group[i]
      } else {
        paste0("Group ", i)
      }
      
      Color[i] <- if (!is.null(c_in) && nzchar(c_in)) {
        c_in
      } else if (i <= nrow(cur) && nzchar(cur$Color[i])) {
        cur$Color[i]
      } else {
        fallback_cols[i]
      }
    }
    
    groups_rv(tibble(Group = Group, Color = Color))
  })
  
  observeEvent(input$auto_colors, {
    df <- groups_rv()
    req(nrow(df) > 0)
    df$Color <- fallback_palette(nrow(df))
    groups_rv(df)
    for (i in seq_len(nrow(df))) {
      updateColourInput(session, paste0("grp_col_", i), value = df$Color[i])
    }
  })
  
  # ---- Assignment table state ----
  assign_df <- reactiveVal(
    tibble(Sample = character(), Group = character(), Color = character())
  )
  
  final_group_df <- reactive({
    df <- assign_df()
    req(nrow(df) > 0)
    
    # Join to current group definitions to ensure colors are in sync
    out <- df %>%
      left_join(groups_rv(), by = "Group") %>%
      transmute(
        Sample,
        Group,
        Color = ifelse(!is.na(Color.y) & nzchar(Color.y),
                       Color.y,
                       Color.x)  # fall back if needed
      )
    
    out
  })
  
  observe({
    req(raw_df())
    samples <- input$sample_cols
    req(length(samples) > 0)
    df <- tibble(
      Sample = samples,
      Group  = NA_character_,
      Color  = NA_character_
    )
    assign_df(df)
    updateSelectInput(session, "bulk_group", choices = groups_rv()$Group)
  })
  
  # ---- Token-based helper text & auto-delim ----
  output$token_help <- renderUI({
    df <- assign_df()
    req(nrow(df) > 0)
    s <- paste(df$Sample, collapse = " ")
    txt <- if (grepl("_", s)) {
      "We detected your sample names use <code>_</code> as a separator. Using that to auto-group. <span style='font-size:11px;color:#777'>(Advanced: edit pattern)</span>"
    } else if (grepl("-", s)) {
      "We detected your sample names use <code>-</code> as a separator. Using that to auto-group. <span style='font-size:11px;color:#777'>(Advanced: edit pattern)</span>"
    } else if (grepl("\\.", s)) {
      "We detected your sample names use <code>.</code> as a separator. Using that to auto-group. <span style='font-size:11px;color:#777'>(Advanced: edit pattern)</span>"
    } else {
      "No clear separator detected. Defaulting to split on <code>_</code>, <code>-</code>, or <code>.</code>. <span style='font-size:11px;color:#777'>(Advanced: edit pattern)</span>"
    }
    HTML(txt)
  })
  
  observeEvent(assign_df(), {
    df <- assign_df()
    if (nrow(df) == 0) return()
    s <- paste(df$Sample, collapse = " ")
    guess <- if (grepl("_", s)) {
      "[-_]+"
    } else if (grepl("-", s)) {
      "[-]+"
    } else if (grepl("\\.", s)) {
      "[.]+"
    } else {
      "[-_.]+"
    }
    cur <- input$token_delim
    if (is.null(cur) || cur == "" || cur == "[-_.]+") {
      updateTextInput(session, "token_delim", value = guess)
    }
  }, ignoreInit = TRUE)
  
  # ---- Token-based grouping ----
  tokens_tbl <- reactive({
    df <- assign_df()
    req(nrow(df) > 0)
    del <- input$token_delim %||% "[-_.]+"
    tibble(Sample = df$Sample) %>%
      mutate(
        .split   = str_split(Sample, pattern = del),
        n_tokens = lengths(.split)
      )
  })
  
  output$token_index_ui <- renderUI({
    tbl <- tokens_tbl()
    req(nrow(tbl) > 0)
    max_tok <- suppressWarnings(max(tbl$n_tokens, na.rm = TRUE))
    if (!is.finite(max_tok) || max_tok < 1) max_tok <- 1
    sliderInput(
      "token_idx", "Token index to group by",
      min = 1, max = max_tok, value = 1, step = 1
    )
  })
  
  observeEvent(input$apply_token, {
    df_assign <- assign_df()
    req(nrow(df_assign) > 0)
    
    idx <- input$token_idx %||% 1
    del <- input$token_delim %||% "[-_.]+"
    
    split_list <- str_split(df_assign$Sample, pattern = del)
    token <- vapply(split_list, function(v){
      if (length(v) >= idx) v[[idx]] else ""
    }, character(1))
    
    if (isTRUE(input$token_drop_singletons)) {
      tab  <- sort(table(token), decreasing = TRUE)
      keep <- names(tab)[tab > 1]
      token[!(token %in% keep)] <- NA_character_
    }
    
    df_assign$Group <- token
    
    grp_df <- groups_rv()
    uniq   <- sort(unique(na.omit(df_assign$Group)))
    if (length(uniq) > 0) {
      if (nrow(grp_df) < length(uniq)) {
        add <- tibble(
          Group = paste0("Group ", (nrow(grp_df) + 1):length(uniq)),
          Color = fallback_palette(length(uniq) - nrow(grp_df))
        )
        grp_df <- bind_rows(grp_df, add)
      }
      color_map <- setNames(fallback_palette(length(uniq)), uniq)
      overlap <- intersect(uniq, grp_df$Group)
      color_map[overlap] <- grp_df$Color[match(overlap, grp_df$Group)]
      df_assign$Color <- ifelse(
        is.na(df_assign$Group),
        NA_character_,
        color_map[df_assign$Group]
      )
    }
    
    assign_df(df_assign)
    updateSelectInput(
      session, "bulk_group",
      choices = sort(unique(na.omit(df_assign$Group)))
    )
    updateTabsetPanel(session, "tabs", selected = "Assignment table")
  })
  
  output$assign_dt <- renderDT({
    df <- final_group_df()
    req(nrow(df) > 0)
    grp_names <- groups_rv()$Group
    
    dt <- datatable(
      df,
      selection = "multiple",
      rownames  = FALSE,
      options   = list(
        scrollX   = TRUE,
        pageLength = 12,
        columnDefs = list(
          list(
            targets = 1,
            render  = DT::JS(sprintf(
              "function(data,type,row,meta){
               var v = data || '';
               var opts = %s;
               var sel = '';
               opts.forEach(function(o){
                 sel += '<option ' + (o==v?'selected':'') + '>' + o + '</option>';
               });
               return type === 'display'
                 ? '<select class=\"dt-group\">' + sel + '</select>'
                 : data;
             }",
             jsonlite::toJSON(grp_names)
            ))
          )
        )
      )
    )
    
    dt %>%
      formatStyle(
        "Color",
        target = "cell",
        backgroundColor = styleEqual(unique(df$Color), unique(df$Color))
      )
  }, server = FALSE)
  
  # bind js handler once
  observe({
    session$sendCustomMessage("bindGroupSelect", list(id = "assign_dt"))
  })
  
  # Handle per-row dropdown change
  observeEvent(input$assign_change, {
    ch <- input$assign_change
    req(ch$row, ch$value)
    df <- assign_df()
    req(ch$row <= nrow(df))
    df$Group[ch$row] <- ch$value
    
    df <- df %>%
      left_join(groups_rv(), by = "Group") %>%
      transmute(
        Sample,
        Group,
        Color = ifelse(!is.na(Color.y) & nzchar(Color.y),
                       Color.y,
                       Color.x %||% Color.y)
      )
    
    assign_df(df)
    replaceData(dataTableProxy("assign_dt"), df,
                resetPaging = FALSE, rownames = FALSE)
  })
  
  observeEvent(input$apply_bulk, {
    rows <- input$assign_dt_rows_selected
    grp  <- input$bulk_group
    req(length(rows) >= 1, !is.null(grp), nzchar(grp))
    
    df <- assign_df()
    req(nrow(df) > 0)
    
    # update selected rows
    df$Group[rows] <- grp
    
    # fill colors from groups_rv
    df <- df %>%
      left_join(groups_rv(), by = "Group") %>%
      transmute(
        Sample,
        Group,
        Color = ifelse(!is.na(Color.y) & nzchar(Color.y),
                       Color.y,
                       Color.x %||% Color.y)
      )
    
    assign_df(df)
    # no replaceData() needed: renderDT() above re-runs automatically from assign_df()
  })
  
  output$download_groupfile <- downloadHandler(
    filename = function() {
      paste0("groupfile_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv")
    },
    content = function(file) {
      df <- final_group_df()
      req(nrow(df) > 0)
      
      # Validation: all assigned?
      if (any(is.na(df$Group) | df$Group == "")) {
        showModal(modalDialog(
          title = "Unassigned samples",
          "Some samples do not have a Group yet. Please assign all samples before downloading.",
          easyClose = TRUE
        ))
        return()
      }
      
      if (any(is.na(df$Color) | df$Color == "")) {
        showModal(modalDialog(
          title = "Missing colors",
          "Some groups do not have a Color defined. Please set colors for all groups before downloading.",
          easyClose = TRUE
        ))
        return()
      }
      
      readr::write_csv(df, file)
    }
  )
  # ---- Template hint ----
  output$template_hint <- renderUI({
    HTML(paste0(
      "<b>Note:</b> We'll automatically detect sample columns as numeric columns in your file. ",
      "If your table includes gene ID columns like <code>EnsemblID</code> and/or <code>Gene_Symbol</code>, ",
      "we'll ignore them for grouping; one or both are acceptable."
    ))
  })
}

# JS handler for the in-table group dropdown
js <- "
Shiny.addCustomMessageHandler('bindGroupSelect', function(opts){
  var id  = opts.id;
  var tbl = $('#' + id + ' table');
  tbl.on('change', 'select.dt-group', function(){
    var rowIdx = $(this).closest('tr').index();
    var val    = $(this).val();
    Shiny.setInputValue('assign_change',
      {row: rowIdx + 1, value: val},
      {priority: 'event'});
  });
});
"

shinyApp(
  ui = tagList(ui, tags$script(HTML(js))),
  server = server
)