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
  library(plotly)
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
                 fileInput(
                   "expression_dataframe", "Choose CSV file",
                   accept = c(".csv","text/csv","text/plain")
                 ),
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
                 tags$hr(),
                 
                 h4("3) Define groups & colors"),
                 numericInput("n_groups", "Number of groups", value = 2, min = 1, step = 1),
                 uiOutput("group_name_color_ui"),
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
                   textInput("token_delim", "Delimiter regex", value = "[_\\.-]+"),
                   helpText("For each group below, choose one or more tokens. A sample is assigned to that group if its name contains all selected tokens."),
                   uiOutput("token_group_selectors"),
                   # checkboxInput("token_drop_singletons", "Hide groups with only 1 sample", value = FALSE),
                   actionButton("apply_token", "Apply token grouping")
                 ),
                 
                 tags$hr(),
                 downloadButton("download_groupfile", "Download groupfile (Sample,Group,Color)")
    ),
    
    mainPanel(width = 8,
              tabsetPanel(id = "tabs",
                          tabPanel("Assignment table",
                                   br(),
                                   tags$p(
                                     "Assign groups with token-based auto-fill, then fine-tune here. Multi-select rows and apply a group in bulk.",
                                     style = "font-size:16px; color:#000000; font-weight:600; margin-bottom:8px;"
                                   ),
                                   fluidRow(
                                     column(6, selectInput("bulk_group", "Set selected rows to group", choices = NULL)),
                                     column(6, actionButton("apply_bulk", "Apply to selected"))
                                   ),
                                   DTOutput("assign_dt")
                          ),
                          tabPanel("Visualization preview",
                                   br(),
                                   helpText("Quick QC: boxplots of top expressed genes using your current group assignments and colors."),
                                   uiOutput("viz_gene_ui"),
                                   plotlyOutput("viz_plot", height = "420px")
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
  
  # ---- Groups: names & colors ----
  
  groups_rv <- reactiveVal(
    tibble(
      Group = paste0("Group ", 1:2),
      Color = c("#E86FDE", "#3CC9B9")
    )
  )
  
  # Adjust number of groups
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
  
  # UI for group names/colors
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
  
  # Sync group labels/colors from UI back into groups_rv (no resizing here)
  observe({
    df <- groups_rv()
    ng <- nrow(df)
    if (!ng) return()
    
    new_G <- df$Group
    new_C <- df$Color
    
    for (i in seq_len(ng)) {
      g_in <- input[[paste0("grp_name_", i)]]
      c_in <- input[[paste0("grp_col_", i)]]
      
      if (!is.null(g_in) && nzchar(g_in)) new_G[i] <- g_in
      if (!is.null(c_in) && nzchar(c_in)) new_C[i] <- c_in
    }
    
    groups_rv(tibble(Group = new_G, Color = new_C))
  })
  
  # Auto colors
  observeEvent(input$auto_colors, {
    df <- groups_rv(); req(nrow(df) > 0)
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
  
  # token_choices <- reactive({
  #   df <- assign_df()
  #   req(nrow(df) > 0)
  #   
  #   del <- input$token_delim %||% "[_\\.-]+"
  #   
  #   toks <- str_split(df$Sample, pattern = del)
  #   vals <- unique(unlist(toks))
  #   vals <- vals[!is.na(vals) & vals != ""]
  #   sort(vals)
  # })
  
  token_choices <- reactive({
    df <- assign_df()
    req(nrow(df) > 0)
    
    del <- input$token_delim %||% "[_\\.-]+"
    
    toks_list <- str_split(df$Sample, pattern = del)
    
    # flatten, clean
    toks <- unlist(lapply(toks_list, function(v){
      v[!is.na(v) & v != ""]
    }))
    
    if (!length(toks)) return(character(0))
    
    # count occurrences
    tab <- table(toks)
    
    # keep only tokens that appear in >= 2 samples
    keep <- names(tab)[tab > 1]
    
    sort(keep)
  })
  
  # Initialize assign_df when samples change
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
  
  # Final mapping used for viz/download
  final_group_df <- reactive({
    df <- assign_df()
    req(nrow(df) > 0)
    
    # Join to current group definitions for authoritative colors
    out <- df %>%
      left_join(groups_rv(), by = "Group") %>%
      mutate(
        Color = ifelse(!is.na(Color.y) & nzchar(Color.y), Color.y, Color.x)
      ) %>%
      select(Sample, Group, Color)
    
    out
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
  
  # ---- Tokenization table ----
  
  tokens_tbl <- reactive({
    df <- assign_df()
    req(nrow(df) > 0)
    
    del <- input$token_delim %||% "[-_.]+"
    
    splits <- str_split(df$Sample, pattern = del)
    max_tokens <- max(lengths(splits))
    
    padded <- lapply(splits, function(v){
      length(v) <- max_tokens
      v
    })
    mat <- do.call(rbind, padded)
    
    out <- as_tibble(mat, .name_repair = ~ paste0("T", seq_along(.)))
    out$Sample <- df$Sample
    out
  })
  
  # Token position selector
  output$token_dims_ui <- renderUI({
    tt <- tokens_tbl()
    req(nrow(tt) > 0)
    
    tcols <- grep("^T[0-9]+$", names(tt), value = TRUE)
    if (!length(tcols)) return(NULL)
    
    choices_list <- lapply(tcols, function(tc){
      vals <- unique(na.omit(tt[[tc]]))
      vals <- vals[vals != ""]
      if (length(vals) <= 1) return(NULL)  # skip non-informative
      
      lbl <- paste0(
        "Position ", sub("T", "", tc), " (",
        paste(head(vals, 4), collapse = " / "),
        if (length(vals) > 4) " / ..." else "",
        ")"
      )
      setNames(tc, lbl)
    })
    
    choices_vec <- do.call(c, choices_list)
    
    if (is.null(choices_vec) || !length(choices_vec)) {
      return(helpText("No informative token positions detected. Adjust the delimiter or sample names."))
    }
    
    selectInput(
      "token_dims",
      "Combine which token positions into groups?",
      choices = choices_vec,
      multiple = TRUE,
      selected = choices_vec[1]
    )
  })
  
  # Apply token-based grouping
  observeEvent(input$apply_token, {
    df <- assign_df()
    req(nrow(df) > 0)
    
    grps <- groups_rv()
    req(nrow(grps) > 0)
    
    del <- input$token_delim %||% "[_\\.-]+"
    
    # Split each sample into its tokens
    toks_list <- str_split(df$Sample, pattern = del)
    sample_tokens <- lapply(toks_list, function(v){
      unique(v[!is.na(v) & v != ""])
    })
    
    # Collect selected tokens per group (one rule per group)
    grp_tokens <- lapply(seq_len(nrow(grps)), function(i){
      sel <- input[[paste0("grp_tokens_", i)]]
      sel <- sel[!is.na(sel) & sel != ""]
      sel
    })
    names(grp_tokens) <- grps$Group
    
    # If no group has tokens, bail nicely
    if (all(vapply(grp_tokens, length, integer(1)) == 0)) {
      showModal(modalDialog(
        title = "No token rules defined",
        "Please select one or more tokens for at least one group before applying.",
        easyClose = TRUE
      ))
      return()
    }
    
    # Assign each sample to at most one group:
    # a sample is in group G if it contains ALL tokens for G.
    # If multiple groups match, first (top-most) group wins.
    assigned_group <- rep(NA_character_, length(sample_tokens))
    
    for (g in names(grp_tokens)) {
      toks <- grp_tokens[[g]]
      if (!length(toks)) next
      
      hit <- vapply(sample_tokens, function(st){
        all(toks %in% st)
      }, logical(1))
      
      # Only fill unassigned samples
      assignable <- is.na(assigned_group) & hit
      assigned_group[assignable] <- g
    }
    
    # Optionally drop groups with only a single sample
    # if (isTRUE(input$token_drop_singletons)) {
    #   tab <- table(assigned_group)
    #   keep <- names(tab)[tab > 1]
    #   assigned_group[!(assigned_group %in% keep)] <- NA_character_
    # }
    
    tab <- table(assigned_group)
    keep <- names(tab)[tab > 1]
    assigned_group[!(assigned_group %in% keep)] <- NA_character_
    
    # Apply to df
    df$Group <- assigned_group
    
    # Map colors: exactly one color per group from groups_rv
    color_map <- setNames(grps$Color, grps$Group)
    df$Color <- ifelse(
      is.na(df$Group),
      NA_character_,
      unname(color_map[df$Group])
    )
    
    assign_df(df)
    
    # Update bulk group dropdown to only show used groups
    used_groups <- sort(unique(na.omit(df$Group)))
    updateSelectInput(session, "bulk_group", choices = used_groups)
    
    # Jump to table so user can see the effect
    updateTabsetPanel(session, "tabs", selected = "Assignment table")
  })
  
  # ---- Assignment table (manual edit + bulk) ----
  
  output$assign_dt <- renderDT({
    df <- assign_df()
    req(nrow(df) > 0)
    
    # Display-only copy
    display_df <- df
    display_df$Group[is.na(display_df$Group) | display_df$Group == ""] <- "Unassigned"
    
    datatable(
      display_df,
      selection = "multiple",
      rownames  = FALSE,
      options   = list(
        scrollX    = TRUE,
        pageLength = 200
      )
    ) %>%
      formatStyle(
        "Color",
        target = "cell",
        backgroundColor = styleEqual(
          unique(df$Color),
          unique(df$Color)
        )
      ) %>%
      formatStyle(
        "Group",
        color = styleEqual("Unassigned", "grey40"),
        fontStyle = styleEqual("Unassigned", "italic")
      )
  }, server = FALSE)
  
  # Bulk apply
  observeEvent(input$apply_bulk, {
    rows <- input$assign_dt_rows_selected
    grp  <- input$bulk_group
    req(length(rows) >= 1, !is.null(grp), nzchar(grp))
    
    df <- assign_df()
    req(nrow(df) > 0)
    
    df$Group[rows] <- grp
    
    df <- df %>%
      left_join(groups_rv(), by = "Group") %>%
      transmute(
        Sample,
        Group,
        Color = ifelse(!is.na(Color.y) & nzchar(Color.y),
                       Color.y,
                       Color.x)
      )
    
    assign_df(df)
  })
  
  # ---- Visualization preview ----
  
  expr_long <- reactive({
    df <- raw_df()
    req(df)
    
    id_cols  <- input$id_cols %||% infer_id_cols(df)
    samples  <- input$sample_cols
    req(length(samples) > 0)
    
    keep_cols <- unique(c(id_cols, samples))
    df <- df[, keep_cols, drop = FALSE]
    
    gene_col <- intersect(
      c("Gene_Symbol", "GeneSymbol", "Symbol", "gene", "Gene"),
      names(df)
    )[1]
    if (is.na(gene_col)) {
      gene_col <- intersect(
        c("EnsemblID", "ensembl_gene_id", "ensembl_id"),
        names(df)
      )[1]
    }
    
    df_long <- df %>%
      pivot_longer(
        cols = all_of(samples),
        names_to  = "Sample",
        values_to = "Count"
      )
    
    df_long$Gene <- if (!is.na(gene_col)) df_long[[gene_col]] else seq_len(nrow(df_long))
    
    df_long
  })
  
  high_expr_data <- reactive({
    expr <- expr_long()
    map  <- final_group_df()
    
    df <- expr %>%
      inner_join(map, by = "Sample") %>%
      filter(!is.na(Group), !is.na(Color))
    
    req(nrow(df) > 0)
    
    df$Count <- suppressWarnings(as.numeric(df$Count))
    
    top_genes <- df %>%
      group_by(Gene) %>%
      summarize(mean_expr = mean(Count, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_expr)) %>%
      slice_head(n = 12) %>%
      pull(Gene)
    
    list(df = df, top_genes = top_genes)
  })
  
  output$viz_gene_ui <- renderUI({
    he <- high_expr_data()
    choices <- he$top_genes
    req(length(choices) > 0)
    
    selectInput(
      "viz_gene",
      "Choose gene to preview (top expressed):",
      choices = choices,
      selected = choices[1]
    )
  })
  
  output$viz_plot <- renderPlotly({
    he <- high_expr_data()
    df <- he$df
    req(nrow(df) > 0)
    
    gene <- input$viz_gene %||% he$top_genes[1]
    req(gene)
    
    map <- final_group_df()
    color_map <- setNames(map$Color, map$Group)
    
    gdat <- df %>% filter(Gene == gene)
    req(nrow(gdat) > 0)
    
    groups <- unique(gdat$Group)
    
    p <- plot_ly()
    for (g in groups) {
      sub <- gdat %>% filter(Group == g)
      col <- unname(color_map[g])
      
      hover_txt <- paste0(
        "<b>Sample</b>: ", sub$Sample,
        "<br><b>Group</b>: ", sub$Group,
        "<br><b>Expr</b>: ", signif(sub$Count, 4)
      )
      
      p <- p %>%
        add_trace(
          data = sub,
          x = ~Group,
          y = ~Count,
          type = "box",
          name = g,
          marker = list(color = col),
          line   = list(color = col),
          boxpoints = "all",
          jitter    = 0.3,
          pointpos  = 0,
          text      = hover_txt,
          hoverinfo = "text",
          showlegend = FALSE
        )
    }
    
    p %>%
      layout(
        title = list(
          text = paste0("Expression preview for <b>", gene, "</b>"),
          x = 0
        ),
        xaxis = list(title = ""),
        yaxis = list(
          title = "Expression (units from uploaded file)",
          zeroline = TRUE
        ),
        boxmode = "group",
        margin = list(t = 60, b = 40, l = 60, r = 20)
      )
  })
  
  output$token_group_selectors <- renderUI({
    req(input$group_mode == "token")
    choices <- token_choices()
    req(length(choices) > 0)
    
    grps <- groups_rv()
    req(nrow(grps) > 0)
    
    tagList(
      lapply(seq_len(nrow(grps)), function(i){
        div(
          tags$b(grps$Group[i]),
          selectizeInput(
            inputId = paste0("grp_tokens_", i),
            label   = NULL,
            choices = choices,
            multiple = TRUE,
            options = list(
              plugins = list("remove_button"),
              placeholder = "Select tokens that define this group"
            )
          ),
          tags$hr(style = "margin:4px 0;")
        )
      })
    )
  })
  
  # ---- Download groupfile ----
  
  output$download_groupfile <- downloadHandler(
    filename = function() {
      paste0("groupfile_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv")
    },
    content = function(file) {
      df <- final_group_df()
      req(nrow(df) > 0)
      
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
}

shinyApp(
  ui = ui,
  server = server
)