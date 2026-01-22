#setwd("/Users/shivam/Desktop/POM_Lab/experiments/Server/app_5_withINFg/")
# Install required packages if not already installed
if(!requireNamespace("shiny", quietly = TRUE)) install.packages("shiny")
if(!requireNamespace("shinythemes", quietly = TRUE)) install.packages("shinythemes")
if(!requireNamespace("dtw", quietly = TRUE)) install.packages("dtw")
if(!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
if(!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if(!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if(!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if(!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if(!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
if(!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if(!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")

# Load Packages
lapply(c("shiny", "shinythemes","dplyr", "tidyverse", 
         "ggplot2", "ggpubr", "stringr", "gridExtra", "dtw", 
         "matrixStats", "data.table"), library, character.only = TRUE)

# Load GSE22307 dataset
data_long <- read.csv("GSE22307_data_long.csv")


# Load GSE131032 dataset
data_long2 <- read.csv("GSE131032_data_long.csv")

# Load GSE282271_data_long

data_long3<-read.csv("GSE282271_data_long.csv")

# Load Inflammation Genes
GIN_26 <- read.csv("dge_update_25jun24.txt", header = FALSE)
GIN26_vals <- str_to_title(GIN_26$V1)

# Function to get inflammation data for any dataset
get_inflam_data <- function(data) {
  data[data$gene %in% GIN26_vals, ]
}




# Plotting function (returns ggplot object)

plot_spline_comparison <- function(gene_list, dl, time_col, custom_label = NULL) {
  
  # -------------------------
  # Basic checks
  # -------------------------
  if (is.null(gene_list) || length(gene_list) == 0) {
    stop("Provide at least one gene.")
  }
  
  required_cols <- c("gene", "samples", time_col, "Zscore_expression")
  if (!all(required_cols %in% colnames(dl))) {
    stop(
      paste(
        "Dataset must contain columns:",
        paste(required_cols, collapse = ", ")
      )
    )
  }
  
  # -------------------------
  # Prepare test expression
  # -------------------------
  plot_data <- dl %>%
    filter(
      gene %in% gene_list,
      !is.na(.data[[time_col]])
    ) %>%
    group_by(samples, .data[[time_col]]) %>%
    summarise(
      Zscore_expression = mean(Zscore_expression, na.rm = TRUE),
      .groups = "drop"
    )
  
  # -------------------------
  # Prepare inflammation genes
  # -------------------------
  infl_data <- get_inflam_data(dl) %>%
    filter(!is.na(.data[[time_col]]))
  
  # -------------------------
  # Label handling
  # -------------------------
  label_text <- if (!is.null(custom_label) && custom_label != "") {
    custom_label
  } else if (length(gene_list) == 1) {
    gene_list
  } else {
    paste0("Avg(", length(gene_list), " genes)")
  }
  
  # -------------------------
  # Safe annotation positions
  # -------------------------
  time_vals <- dl[[time_col]]
  time_vals <- time_vals[!is.na(time_vals)]
  has_time <- length(time_vals) > 0
  x_max <- if (has_time) max(time_vals) else NA
  
  y_vals <- dl$Zscore_expression
  y_vals <- y_vals[!is.na(y_vals)]
  has_y <- length(y_vals) > 0
  y_max <- if (has_y) max(y_vals) else NA
  
  # -------------------------
  # Base plot
  # -------------------------
  p <- ggplot() +
    geom_smooth(
      data = plot_data,
      aes(x = .data[[time_col]], y = Zscore_expression),
      formula = y ~ s(x, bs = "cs", k = 4),
      method = "gam",
      color = "darkgreen",
      fill = "lightgreen",
      se = TRUE
    ) +
    geom_smooth(
      data = infl_data,
      aes(x = .data[[time_col]], y = Zscore_expression),
      formula = y ~ s(x, bs = "cs", k = 4),
      method = "gam",
      color = "red",
      fill = "pink",
      se = TRUE
    ) +
    theme_classic() +
    theme(
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    labs(
      x = ifelse(
        time_col == "hours_of_IFNg_treatment",
        "Hours of IFNγ treatment",
        "Day of DSS treatment"
      ),
      y = "Z-score expression",
      title = paste("Spline:", label_text)
    )
  
  # -------------------------
  # Safe annotations
  # -------------------------
  if (has_time && has_y) {
    p <- p +
      annotate(
        "text",
        x = x_max,
        y = y_max * 0.95,
        label = "Inflammation gene expression",
        color = "red",
        hjust = 1,
        fontface = "bold",
        size = 10
      ) +
      annotate(
        "text",
        x = x_max,
        y = y_max * 0.75,
        label = label_text,
        color = "darkgreen",
        hjust = 1,
        fontface = "bold",
        size = 10
      )
  }
  return(p)
}



# Function to calculate spline cureve stats

calculate_expression_stats <- function(gene_list, dl, time_col) {
  
  inf_genes <- get_inflam_data(dl)
  
  # Filter test genes
  plot_data <- dl %>%
    filter(gene %in% gene_list, !is.na(.data[[time_col]])) %>%
    group_by(samples, .data[[time_col]]) %>%
    summarise(Zscore_expression = mean(Zscore_expression, na.rm = TRUE), .groups = "drop")
  
  # Mean spline per time
  test_spline <- plot_data %>%
    group_by(.data[[time_col]]) %>%
    summarise(expr_test = mean(Zscore_expression, na.rm = TRUE), .groups = "drop")
  
  infl_spline <- inf_genes %>%
    filter(!is.na(.data[[time_col]])) %>%
    group_by(.data[[time_col]]) %>%
    summarise(expr_inflam = mean(Zscore_expression, na.rm = TRUE), .groups = "drop")
  
  # Merge on time
  merged <- inner_join(test_spline, infl_spline, by = time_col)
  
  # If no overlap, return informative table
  if (nrow(merged) < 2) {
    return(
      data.frame(
        Parameter = "Error",
        Value = paste(
          "No overlapping time points between test genes and inflammation genes for",
          time_col
        )
      )
    )
  }
  
  # ---- SAFE STATISTICS ----
  mae <- mean(abs(merged$expr_test - merged$expr_inflam))
  global_var_test <- var(merged$expr_test)
  global_var_inflam <- var(merged$expr_inflam)
  
  dtw_dist <- dtw::dtw(merged$expr_test, merged$expr_inflam)$distance
  
  var_each_timepoint <- apply(
    merged[, c("expr_test", "expr_inflam")],
    1,
    var
  )
  
  spearman_rho <- suppressWarnings(
    cor(merged$expr_test, merged$expr_inflam, method = "spearman")
  )
  
  data.frame(
    Parameter = c(
      "Mean Absolute Error",
      "Global Variance (Test)",
      "Global Variance (Inflammation)",
      "DTW Distance",
      paste0("Variance (", merged[[time_col]], ")"),
      "Spearman Correlation (ρ)"
    ),
    Value = c(
      round(mae, 4),
      round(global_var_test, 4),
      round(global_var_inflam, 4),
      round(dtw_dist, 4),
      round(var_each_timepoint, 4),
      round(spearman_rho, 4)
    )
  )
}




# Shiny UI
ui <- fluidPage(
  theme = shinythemes::shinytheme("cosmo"),
  titlePanel("Gene Expression Spline Viewer (Multiple datasets)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset_choice", "Choose dataset:",
                  choices = c("GSE22307 (DSS treatment)" = "data_long", 
                              "GSE131032 (DSS treatment and recovery)" = "data_long2",
                              "GSE282271 (INFg treatment)"="data_long3"),
                  selected = "data_long"),
      radioButtons("mode", "Input mode:", choices = c("Single Gene", "Upload Gene List"), selected = "Single Gene"),
      conditionalPanel(
        condition = "input.mode == 'Single Gene'",
        uiOutput("gene_selector_ui")
      ),
      conditionalPanel(
        condition = "input.mode == 'Upload Gene List'",
        fileInput("gene_file", "Upload CSV (first column = gene names)", accept = c(".csv")),
        textInput("custom_label", "Custom label for test expression (optional)", value = "Average Expression")
      ),
      actionButton("refresh_choices", "Refresh gene list"),
      br(), br(),
      downloadButton("downloadPlot", "Download Plot (PNG)"),
      downloadButton("downloadStats", "Download Stats (CSV)")
    ),
    mainPanel(
      plotOutput("splinePlot", width = "800px", height = "600px"),
      h4("Statistical Comparison"),
      tableOutput("statsTable")
    )
  )
)

# Shiny Server
server <- function(input, output, session) {
  # reactive dataset
  ds <- reactive({
    if (input$dataset_choice == "data_long") {data_long} else if (input$dataset_choice == "data_long2") {data_long2} else {data_long3}
    
  })
  
  # time point selector
  time_var <- reactive({
    if (input$dataset_choice == "data_long3") {
      "hours_of_IFNg_treatment"
    } else {
      "day_of_DSS_treatment"
    }
  })
  
  # gene selector UI (populated from chosen dataset)
  output$gene_selector_ui <- renderUI({
    df <- ds()
    genes <- unique(df$gene)
    selectInput("selected_gene", "Select gene:", choices = genes, selected = genes[1], multiple = FALSE, selectize = TRUE)
  })
  
  # allow refresh (useful if datasets changed externally)
  observeEvent(input$refresh_choices, {
    updateSelectInput(session, "selected_gene", choices = unique(ds()$gene))
  })
  
  # reactive: list of genes to analyze
  selected_genes <- reactive({
    if (input$mode == "Single Gene") {
      req(input$selected_gene)
      return(c(input$selected_gene))
    } else {
      req(input$gene_file)
      genes <- read.csv(input$gene_file$datapath, header = TRUE, stringsAsFactors = FALSE)[[1]]
      # If user supplied a custom label we won't change gene_list, the label is used for plot title
      return(as.character(genes))
    }
  })
  
  # reactive: computed plot
  plot_obj <- reactive({
    req(selected_genes())
    dl <- ds()
    # guard for missing expected columns
    if (!"Zscore_expression" %in% colnames(dl)) {
      stop("The selected dataset must contain 'Zscore_expression' column.")
    }
    p <- plot_spline_comparison(
      selected_genes(),
      dl,
      time_col = time_var(),
      custom_label = if (input$mode == "Upload Gene List") input$custom_label else NULL
    )
    # if custom label provided for uploaded lists, replace title
    if (input$mode == "Upload Gene List" && !is.null(input$custom_label) && input$custom_label != "") {
      p <- p + ggtitle(paste("Spline:", input$custom_label))
    }
    return(p)
  })
  
  # render plot
  output$splinePlot <- renderPlot({
    p <- plot_obj()
    print(p)
  })
  
  # reactive: stats table
  stats_tbl <- reactive({
    req(selected_genes())
    dl <- ds()
    calculate_expression_stats(selected_genes(), dl, time_var())
  })
  
  output$statsTable <- renderTable({
    stats_tbl()
  })
  
  # downloads
  output$downloadStats <- downloadHandler(
    filename = function() {
      suffix <- if (input$mode == "Single Gene") input$selected_gene else "uploaded_genes"
      paste0("stats_", input$dataset_choice, "_", suffix, ".csv")
    },
    content = function(file) {
      write.csv(stats_tbl(), file, row.names = FALSE)
    }
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      suffix <- if (input$mode == "Single Gene") input$selected_gene else "uploaded_genes"
      paste0("spline_", input$dataset_choice, "_", suffix, ".png")
    },
    content = function(file) {
      ggsave(filename = file, plot = plot_obj(), device = "png", width = 8, height = 6, dpi = 300)
    }
  )
  # SESSION TIMEOUT CODE
  # Reactive timer and reset trigger
  timeout_timer <- reactiveVal(Sys.time())
  
  # Watch input for activity (reset timeout)
  observe({
    input$selected_gene
    input$gene_file
    input$custom_label
    input$mode
    
    # Reset inactivity timer
    timeout_timer(Sys.time())
  })
  
  # Session killer if 1 minutes passed without activity
  observe({
    invalidateLater(10000, session)  # Check every 10 seconds
    if (difftime(Sys.time(), timeout_timer(), units = "secs") > 600) {
      session$close()
    }
  })
}



# Run app
shinyApp(ui = ui, server = server)


