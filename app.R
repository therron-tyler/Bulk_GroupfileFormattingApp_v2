# Tyler Therron
# Winter Lab Macrophage Genomics - Rheumatology Department - Feinberg School of Medicine
# Group File Maker App for Winter Lab Database of Bulk RNA-sequencing Database (v2)
# 2025/11/05

# Key upgrades vs. your v1
# - Accepts inputs with EnsemblID + Gene Symbol, EnsemblID-only, Gene Symbol-only, or neither.
# - Robust sample detection (numeric columns by default; you can override).
# - Faster grouping without drag-and-drop: three options
#   (A) Token-based auto-grouping from sample names (split by delimiters like '_' '-' '.')
#   (B) Regex rule builder (pattern → group name → color)
#   (C) Manual table edit with dropdowns (bulk assign via multi-select)
# - Import an external mapping CSV (Sample,Group,Color) and/or export your current mapping.
# - Live preview of group balance; color palette auto-assignment with easy override.

suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(colourpicker)
})

options(shiny.maxRequestSize = 1000 * 1024^2)

# ----------------------- Helpers ----------------------- #
shorten <- function(x, n=40){
  ifelse(nchar(x) > n, paste0(substr(x,1,n-3), "..."), x)
}

fallback_palette <- function(n){
  # Distinct-ish palette without extra deps
  hues <- seq(15, 375, length = n + 1)
  grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}

infer_id_cols <- function(df){
  cand <- names(df)
  # typical names to treat as gene-id columns (non-sample)
  patt <- c("^symbol$", "^gene_?symbol$", "^ensembl_gene_id$")
  id_cols <- cand[sapply(cand, function(nm){ any(grepl(paste(patt, collapse='|'), nm, ignore.case = TRUE)) })]
  id_cols
}

numeric_sample_cols <- function(df, id_cols){
  nums <- names(df)[vapply(df, is.numeric, logical(1))]
  setdiff(nums, id_cols)
}

apply_regex_rules <- function(assign_df, rules){
  # rules: data.frame(pattern, group, color) in order
  if(is.null(rules) || nrow(rules) == 0) return(assign_df)
  for(i in seq_len(nrow(rules))){
    pat <- rules$pattern[i]
    grp <- rules$group[i]
    col <- rules$color[i]
    if(!isTRUE(nzchar(pat)) || !isTRUE(nzchar(grp))) next
    idx <- grepl(pat, assign_df$Sample, ignore.case = TRUE)
    if(!any(idx)) next
    assign_df$Group[idx] <- grp
    if(isTRUE(nzchar(col))){
      assign_df$Color[idx] <- col
    }
  }
  assign_df
}

`%||%` <- function(x, y) if (is.null(x) || (is.character(x) && length(x) == 0)) y else x

# ----------------------- UI ----------------------- #
ui <- fluidPage(
  titlePanel("Winter Lab Bulk RNA‑seq: Groupfile Builder (v2)"),
  tags$hr(),
  
  sidebarLayout(
    sidebarPanel(width = 4,
                 h4("1) Upload dataset"),
                 fileInput("expression_dataframe", "Choose CSV file", accept = c(".csv","text/csv","text/plain")),
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
                 radioButtons("group_mode", NULL, inline = FALSE,
                              choices = c("Token-based (auto)" = "token",
                                          "Regex rules" = "regex",
                                          "Manual table" = "manual",
                                          "Import mapping CSV" = "import"),
                              selected = "token"),
                 
                 conditionalPanel(
                   condition = "input.group_mode == 'token'",
                   textInput("token_delim", "Delimiter regex", value = "[-_.]+"),
                   uiOutput("token_index_ui"),
                   checkboxInput("token_drop_singletons", "Hide tokens with only 1 sample", value = FALSE),
                   actionButton("apply_token", "Apply token grouping")
                 ),
                 
                 conditionalPanel(
                   condition = "input.group_mode == 'regex'",
                   uiOutput("regex_rule_rows"),
                   actionButton("add_rule", "+ Add rule"),
                   actionButton("clear_rules", "Clear all rules"),
                   actionButton("apply_regex", "Apply regex rules")
                 ),
                 
                 conditionalPanel(
                   condition = "input.group_mode == 'import'",
                   fileInput("mapping_csv", "Upload mapping CSV (Sample,Group,Color)", accept = ".csv"),
                   actionButton("apply_import", "Apply imported mapping")
                 ),
                 
                 tags$hr(),
                 downloadButton("download_mapping", "Download current mapping CSV"),
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
                                   helpText("Assign groups: choose a method at left, then fine-tune here. You can also multi-select rows and set group in bulk."),
                                   fluidRow(
                                     column(6, selectInput("bulk_group", "Set selected rows to group", choices = NULL)),
                                     column(6, actionButton("apply_bulk", "Apply to selected"))
                                   ),
                                   DTOutput("assign_dt")
                          ),
                          tabPanel("Group summary",
                                   br(),
                                   plotOutput("group_bar", height = 320),
                                   DTOutput("group_counts")
                          )
              )
    )
  )
)

# ----------------------- Server ----------------------- #
server <- function(input, output, session){
  
  # Reactive: read main data
  raw_df <- reactive({
    req(input$expression_dataframe)
    readr::read_csv(input$expression_dataframe$datapath, show_col_types = FALSE)
  })
  
  # Show a compact head preview
  output$head_table <- renderDT({
    req(raw_df())
    DT::datatable(head(raw_df(), 10), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # ID column detection message & choice (purely informational; mapping doesn't require IDs)
  output$id_col_detect <- renderUI({
    req(raw_df())
    df <- raw_df()
    ids <- infer_id_cols(df)
    if(length(ids) == 0){
      div(style = "margin-top:-10px;", helpText("No explicit gene ID columns detected. That's ok."))
    } else {
      checkboxGroupInput("id_cols", "Detected gene ID columns (optional):",
                         choices = ids, selected = ids)
    }
  })
  
  # Sample column picker (defaults to numeric non-ID columns)
  output$sample_picker <- renderUI({
    req(raw_df())
    id_cols <- input$id_cols %||% infer_id_cols(raw_df())
    num_cols <- numeric_sample_cols(raw_df(), id_cols)
    # Fallback: if no numeric columns, offer all non-ID columns
    choices <- if(length(num_cols)) num_cols else setdiff(names(raw_df()), id_cols)
    selectizeInput("sample_cols", "Sample columns", choices = choices, selected = choices,
                   multiple = TRUE, options = list(plugins = list("remove_button")))
  })
  
  observeEvent(input$toggle_all, {
    req(input$sample_cols)
    curr <- input$sample_cols
    # simpler: just toggle between none and all currently offered
    offered <- isolate({ req(raw_df()); id_cols <- input$id_cols %||% infer_id_cols(raw_df());
    num_cols <- numeric_sample_cols(raw_df(), id_cols);
    if(length(num_cols)) num_cols else setdiff(names(raw_df()), id_cols) })
    if(length(curr) < length(offered)){
      updateSelectizeInput(session, "sample_cols", selected = offered, server = TRUE)
    } else {
      updateSelectizeInput(session, "sample_cols", selected = character(0), server = TRUE)
    }
  })
  
  # Group names & colors
  groups_rv <- reactiveVal(tibble(Group = paste0("Group ", seq_len(2)), Color = c("#E86FDE", "#3CC9B9")))
  
  observeEvent(input$n_groups, {
    ng <- input$n_groups %||% 1
    cur <- groups_rv()
    if(nrow(cur) < ng){
      add <- tibble(Group = paste0("Group ", (nrow(cur)+1):ng), Color = fallback_palette(ng - nrow(cur)))
      groups_rv(bind_rows(cur, add))
    } else if(nrow(cur) > ng){
      groups_rv(cur %>% slice(1:ng))
    }
  }, ignoreInit = TRUE)
  
  output$group_name_color_ui <- renderUI({
    df <- groups_rv()
    req(nrow(df) >= 1)
    tagList(
      lapply(seq_len(nrow(df)), function(i){
        fluidRow(
          column(8, textInput(paste0("grp_name_", i), label = paste("Group", i, "name"), value = df$Group[i])),
          column(4, colourInput(paste0("grp_col_", i), label = "Color", value = df$Color[i], showColour = "both"))
        )
      })
    )
  })
  

  observe({
    req(input$n_groups)
    ng <- input$n_groups
    
    # get current groups for fallback
    cur <- groups_rv()
    fallback_cols <- fallback_palette(max(ng, 1))
    
    Group <- character(ng)
    Color <- character(ng)
    
    for (i in seq_len(ng)) {
      # fallback to existing value or "Group i"
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
    df <- groups_rv(); req(nrow(df) > 0)
    df$Color <- fallback_palette(nrow(df))
    groups_rv(df)
    # push back to inputs
    for(i in seq_len(nrow(df))){
      updateColourInput(session, paste0("grp_col_", i), value = df$Color[i])
    }
  })
  
  # Assignment data.frame (Sample, Group, Color)
  assign_df <- reactiveVal(tibble(Sample = character(), Group = character(), Color = character()))
  
  # Build initial assignment table when sample cols change
  observe({
    req(raw_df())
    samples <- input$sample_cols
    req(length(samples) >= 1)
    df <- tibble(Sample = samples, Group = NA_character_, Color = NA_character_)
    assign_df(df)
    updateSelectInput(session, "bulk_group", choices = groups_rv()$Group)
  })
  
  # ----- Token-based grouping ----- #
  tokens_tbl <- reactive({
    req(assign_df())
    del <- input$token_delim %||% "[-_.]+"
    tibble(Sample = assign_df()$Sample) %>%
      mutate(.split = str_split(Sample, pattern = del, simplify = FALSE),
             n_tokens = lengths(.split))
  })
  
  output$token_index_ui <- renderUI({
    req(tokens_tbl())
    max_tok <- max(tokens_tbl()$n_tokens)
    sliderInput("token_idx", "Token index to group by", min = 1, max = max(1, max_tok), value = 1, step = 1)
  })
  
  observeEvent(input$apply_token, {
    req(tokens_tbl())
    idx <- input$token_idx %||% 1
    del <- input$token_delim %||% "[-_.]+"
    tmp <- tibble(Sample = assign_df()$Sample)
    split_list <- str_split(tmp$Sample, pattern = del)
    token <- vapply(split_list, function(v){ if(length(v) >= idx) v[[idx]] else "" }, character(1))
    df <- assign_df()
    
    if(isTRUE(input$token_drop_singletons)){
      tab <- sort(table(token), decreasing = TRUE)
      keep <- names(tab)[tab > 1]
      token[!(token %in% keep)] <- NA_character_
    }
    
    # Map tokens to current groups if names match; else create/rename groups
    df$Group <- token
    
    # Assign colors: match group names to groups_rv where possible; otherwise palette
    grp_df <- groups_rv()
    uniq <- sort(unique(na.omit(df$Group)))
    if(length(uniq) > 0){
      # If not enough groups in config, expand silently
      if(nrow(grp_df) < length(uniq)){
        add <- tibble(Group = paste0("Group ", (nrow(grp_df)+1):(length(uniq))),
                      Color = fallback_palette(length(uniq) - nrow(grp_df)))
        grp_df <- bind_rows(grp_df, add)
      }
      # Align by name where possible, else fill sequentially
      color_map <- setNames(rep(fallback_palette(length(uniq)), length(uniq)), uniq)
      # overwrite with configured ones where names match
      overlap <- intersect(uniq, grp_df$Group)
      color_map[overlap] <- grp_df$Color[match(overlap, grp_df$Group)]
      df$Color <- ifelse(is.na(df$Group), NA_character_, color_map[df$Group])
    }
    
    assign_df(df)
    updateSelectInput(session, "bulk_group", choices = sort(unique(na.omit(df$Group))))
    updateTabsetPanel(session, "tabs", selected = "Assignment table")
  })
  
  # ----- Regex rules ----- #
  rules <- reactiveVal(tibble(pattern = character(), group = character(), color = character()))
  
  output$regex_rule_rows <- renderUI({
    r <- rules()
    if(nrow(r) == 0){
      fluidRow(column(12, helpText("Add one or more rules. Order matters; later rules can overwrite earlier ones.")))
    }
    tagList(
      lapply(seq_len(nrow(r)), function(i){
        fluidRow(
          column(5, textInput(paste0("rule_pat_", i), "regex pattern", value = r$pattern[i], placeholder = "e.g., CD11clo|CD11chi")),
          column(4, textInput(paste0("rule_grp_", i), "group name", value = r$group[i])),
          column(3, colourInput(paste0("rule_col_", i), "color", value = r$color[i]))
        )
      })
    )
  })
  
  observeEvent(input$add_rule, {
    r <- rules()
    rules(bind_rows(r, tibble(pattern = "", group = "", color = "")))
  })
  observeEvent(input$clear_rules, { rules(tibble(pattern = character(), group = character(), color = character())) })
  
  # Keep rules in sync with inputs
  observe({
    r <- rules(); if(nrow(r) == 0) return()
    for(i in seq_len(nrow(r))){
      r$pattern[i] <- input[[paste0("rule_pat_", i)]] %||% r$pattern[i]
      r$group[i]   <- input[[paste0("rule_grp_", i)]] %||% r$group[i]
      r$color[i]   <- input[[paste0("rule_col_", i)]] %||% r$color[i]
    }
    rules(r)
  })
  
  observeEvent(input$apply_regex, {
    df <- assign_df(); req(nrow(df) > 0)
    r <- rules(); req(nrow(r) > 0)
    df <- apply_regex_rules(df, r)
    # fill colors from groups_rv if missing
    grp_df <- groups_rv()
    df <- df %>% left_join(grp_df, by = c("Group" = "Group"), keep = FALSE, suffix = c("", ".cfg")) %>%
      mutate(Color = ifelse(!is.na(Color) & nzchar(Color), Color, Color.cfg)) %>%
      select(Sample, Group, Color)
    assign_df(df)
    updateSelectInput(session, "bulk_group", choices = groups_rv()$Group)
    updateTabsetPanel(session, "tabs", selected = "Assignment table")
  })
  
  # ----- Import mapping CSV ----- #
  observeEvent(input$apply_import, {
    req(input$mapping_csv)
    map <- readr::read_csv(input$mapping_csv$datapath, show_col_types = FALSE)
    stopifnot(all(c("Sample","Group","Color") %in% names(map)))
    df <- assign_df(); req(nrow(df) > 0)
    df <- df %>% select(Sample) %>% left_join(map, by = "Sample")
    assign_df(df)
    updateSelectInput(session, "bulk_group", choices = sort(unique(na.omit(df$Group))))
    updateTabsetPanel(session, "tabs", selected = "Assignment table")
  })
  
  # ----- Manual table + bulk set ----- #
  output$assign_dt <- renderDT({
    df <- assign_df(); req(nrow(df) > 0)
    grp_names <- groups_rv()$Group
    # build a select editor for group column
    dt <- datatable(df, selection = "multiple", rownames = FALSE, options = list(scrollX = TRUE, pageLength = 12,
                                                                                 columnDefs = list(list(targets = 1, render = DT::JS(
                                                                                   sprintf("function(data,type,row,meta){var v=data||'';var opts=%s;var sel='';opts.forEach(function(o){sel+=('<option '+(o==v?'selected':'')+'>'+o+'</option>')});return type==='display'?('<select class=\"dt-group\">'+sel+'</select>'):data;}", jsonlite::toJSON(grp_names))
                                                                                 )))))
    dt %>% formatStyle("Color", target = "cell",
                       backgroundColor = styleEqual(groups_rv()$Color, groups_rv()$Color))
  }, server = FALSE)
  
  # Capture group selection changes from the select dropdowns in DT
  observe({
    sess <- session
    proxy <- dataTableProxy("assign_dt")
    # js will send changes via Shiny.onInputChange('assign_change', {row:..., value:...})
  })
  
  # JavaScript to handle select changes inside the DT
  observe({
    session$sendCustomMessage("bindGroupSelect", list(id = "assign_dt"))
  })
  
  # Receive changes
  observeEvent(input$assign_change, {
    ch <- input$assign_change; req(ch$row, ch$value)
    df <- assign_df(); req(ch$row <= nrow(df))
    df$Group[ch$row] <- ch$value
    # color by configured table
    df <- df %>% left_join(groups_rv(), by = c("Group" = "Group")) %>% mutate(Color = Color.y %||% Color.x) %>%
      transmute(Sample, Group, Color)
    assign_df(df)
    replaceData(dataTableProxy("assign_dt"), df, resetPaging = FALSE, rownames = FALSE)
  })
  
  # Bulk apply
  observeEvent(input$apply_bulk, {
    req(input$assign_dt_rows_selected)
    grp <- input$bulk_group %||% return(NULL)
    df <- assign_df()
    idx <- input$assign_dt_rows_selected
    df$Group[idx] <- grp
    df <- df %>% left_join(groups_rv(), by = c("Group" = "Group")) %>% mutate(Color = Color.y %||% Color.x) %>%
      transmute(Sample, Group, Color)
    assign_df(df)
    replaceData(dataTableProxy("assign_dt"), df, resetPaging = FALSE, rownames = FALSE)
  })
  
  # Group summary
  output$group_counts <- renderDT({
    df <- assign_df(); req(nrow(df) > 0)
    cnt <- df %>% count(Group, sort = TRUE) %>% mutate(Color = groups_rv()$Color[match(Group, groups_rv()$Group)])
    datatable(cnt, options = list(dom = 't', pageLength = 100), rownames = FALSE)
  })
  
  output$group_bar <- renderPlot({
    df <- assign_df(); req(nrow(df) > 0)
    cnt <- df %>% count(Group) %>% arrange(desc(n))
    par(mar = c(5,10,2,1))
    barplot(cnt$n, names.arg = cnt$Group, horiz = TRUE, las = 1)
  })
  
  # Downloads
  output$download_mapping <- downloadHandler(
    filename = function(){ paste0("mapping_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv") },
    content = function(file){
      df <- assign_df(); req(nrow(df) > 0)
      write_csv(df, file)
    }
  )
  
  output$download_groupfile <- downloadHandler(
    filename = function(){ paste0("groupfile_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv") },
    content = function(file){
      df <- assign_df(); req(nrow(df) > 0)
      # Ensure required columns and order
      out <- df %>% transmute(Sample, Group, Color)
      # Basic validation
      if(any(!nzchar(out$Sample))) stop("Empty sample names present.")
      if(any(!nzchar(out$Group)))  stop("Some samples are unassigned to any group.")
      if(any(!nzchar(out$Color)))  stop("Some samples are missing colors.")
      write_csv(out, file)
    }
  )
  
  # Template hint (explains ID handling)
  output$template_hint <- renderUI({
    HTML(
      paste0(
        "<b>Note:</b> We'll automatically detect sample columns as numeric columns in your file. ",
        "If your table includes gene ID columns like <code>EnsemblID</code> and/or <code>Gene_Symbol</code>, ",
        "we'll ignore them for grouping; both one or both are acceptable."
      )
    )
  })
  
}

# Widget bindings for DT select-in-cell
js <- "
Shiny.addCustomMessageHandler('bindGroupSelect', function(opts){
  var id = opts.id;
  var tbl = $('#' + id + ' table');
  tbl.on('change', 'select.dt-group', function(){
    var rowIdx = $(this).closest('tr').index();
    var val = $(this).val();
    Shiny.setInputValue('assign_change', {row: rowIdx + 1, value: val}, {priority: 'event'});
  });
});
"

# Launch
shinyApp(
  ui = tagList(ui, tags$script(HTML(js))),
  server = server
)
