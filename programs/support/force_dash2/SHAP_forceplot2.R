# =========================================================
#  Libraries
# =========================================================
library(shiny)
library(ggplot2)
library(dplyr)
library(plotly)
library(forcats)
library(readxl)
library(scales)

# =========================================================
#  Load Data and Mapping
# =========================================================
test       <- readRDS("test.rds")
test_baked       <- readRDS("test_baked.rds")
shap_values_shap_avg      <- readRDS("shap_values_shap_avg.rds")
dummy_list <- read_excel("dummy_namelist.xlsx",col_names = FALSE)

colnames(dummy_list) <- c("dummy", "group", "display","rename")

group_display_map <- dummy_list %>%
  select(group, display) %>%
  distinct()

test_baked_clean <- test_baked %>%
  select(-any_of("ar_primary"))

stopifnot(ncol(shap_values_shap_avg) == ncol(test_baked_clean))
stopifnot(all(colnames(shap_values_shap_avg) == colnames(test_baked_clean)))

# =========================================================
#  Identify variable types
# =========================================================
cont_vars <- c("pir", "LBXEOPCT", "nlr", "LBXIGE", "crp_mg_l", "cot_ng_ml")

cat_vars <- setdiff(
  names(test)[sapply(test, function(x) is.factor(x) || is.character(x))],
  c("ar_primary", cont_vars)
)

# =========================================================
#  Helper: Safe label mapping
# =========================================================
make_safe_labels <- function(vars) {
  dummy_map <- dummy_list %>%
    mutate(group_std = tolower(group),
           display_std = display)
  
  labels <- sapply(vars, function(v) {
    v_std <- tolower(v)
    match_row <- dummy_map %>%
      filter(group_std == v_std)
    
    if (nrow(match_row) > 0 && !is.na(match_row$display_std[1])) {
      return(match_row$display_std[1])
    } else {
      return(v)
    }
  })
  
  setNames(vars, labels)
}

# =========================================================
#  Helper: Custom ribbon plotting without overlap
# =========================================================
geom_force_ribbon <- function(df, fill_pos = "#ef8a62", fill_neg = "#67a9cf") {
  df$sign_grp <- cumsum(c(0, diff(sign(df$cum_shap))) != 0)
  segs <- split(df, df$sign_grp)
  p <- ggplot()
  for (seg in segs) {
    if (mean(seg$cum_shap, na.rm = TRUE) > 0) {
      p <- p + geom_ribbon(
        data = seg,
        aes(x = x_value, ymin = 0, ymax = cum_shap),
        fill = fill_pos, alpha = 0.7
      )
    } else {
      p <- p + geom_ribbon(
        data = seg,
        aes(x = x_value, ymin = cum_shap, ymax = 0),
        fill = fill_neg, alpha = 0.7
      )
    }
  }
  return(p)
}

# =========================================================
#  UI
# =========================================================
ui <- fluidPage(
  titlePanel("Grouped SHAP Visualization Dashboard"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "feature_display", "Select Feature Group:",
        choices = setNames(group_display_map$group, group_display_map$display),
        selected = group_display_map$group[1]
      ),
      selectInput(
        "x_var_cont", "Select Continuous X-axis Variable:",
        choices = make_safe_labels(cont_vars),
        selected = "cot_ng_ml"
      ),
      selectInput(
        "x_var_cat", "Select Categorical X-axis Variable:",
        choices = make_safe_labels(cat_vars),
        selected = cat_vars[1]
      ),
      helpText("① Force-like SHAP = cumulative trajectory along samples\n
                ② Aggregated Trend = mean SHAP per X\n
                ③ Categorical Variables = mean SHAP by factor levels\n
                ④ Summary Importance = top mean |SHAP| features")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("① Force-like SHAP", plotlyOutput("forcePlot", height = "750px")),
        tabPanel("② Aggregated Trend", plotlyOutput("aggPlot", height = "750px")),
        tabPanel("③ Categorical Variables", plotlyOutput("catPlot", height = "750px")),
        tabPanel("④ Summary Importance", plotlyOutput("summaryPlot", height = "750px"))
      )
    )
  )
)

# =========================================================
#  SERVER
# =========================================================
server <- function(input, output, session) {
  
  # Group-level SHAP mean across dummies
  df_filtered <- reactive({
    group_name <- input$feature_display
    related_dummies <- dummy_list$dummy[dummy_list$group == group_name]
    shap_df <- shap_values_shap_avg[, related_dummies[related_dummies %in% colnames(shap_values_shap_avg)], drop = FALSE]
    y_agg <- rowMeans(shap_df, na.rm = TRUE)
    tibble(observation = seq_along(y_agg), y_value = y_agg)
  })
  
  # =====================================================
  # (1) Force-like SHAP Plot
  # =====================================================
  output$forcePlot <- renderPlotly({
    df_plot <- df_filtered()
    req(input$x_var_cont)
    x_var <- input$x_var_cont
    df_plot$x_value <- test_baked_clean[[x_var]]
    
    display_name <- group_display_map$display[group_display_map$group == input$feature_display]
    if (length(display_name) == 0) display_name <- input$feature_display
    
    x_display_label <- names(make_safe_labels(x_var))
    
    df_cont <- df_plot %>%
      arrange(x_value) %>%
      mutate(cum_shap = cumsum(y_value))
    
    max_point <- df_cont %>% slice(which.max(cum_shap))
    min_point <- df_cont %>% slice(which.min(cum_shap))
    
    y_min <- min(df_cont$cum_shap, na.rm = TRUE)
    y_max <- max(df_cont$cum_shap, na.rm = TRUE)
    y_pad <- 0.25 * (y_max - y_min)
    
    p <- geom_force_ribbon(df_cont) +
      geom_line(data = df_cont, aes(x = x_value, y = cum_shap),
                color = "gray40", linewidth = 0.4) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      geom_point(data = max_point, aes(x = x_value, y = cum_shap),
                 color = "#e41a1c", size = 3) +
      geom_point(data = min_point, aes(x = x_value, y = cum_shap),
                 color = "#377eb8", size = 3) +
      labs(
        title = paste0("Force-like cumulative SHAP for ", display_name),
        x = x_display_label, y = "Cumulative SHAP"
      ) +
      coord_cartesian(ylim = c(y_min - y_pad, y_max + y_pad)) +
      theme_minimal(base_size = 16)
    
    ggplotly(p)
  })

  # =====================================================
  # (2) Aggregated Trend Plot (with red/blue points)
  # =====================================================
  output$aggPlot <- renderPlotly({
    df_plot <- df_filtered()
    req(input$x_var_cont)
    x_var <- input$x_var_cont
    df_plot$x_value <- test_baked_clean[[x_var]]
    
    display_name <- group_display_map$display[group_display_map$group == input$feature_display]
    if (length(display_name) == 0) display_name <- input$feature_display
    
    x_display_label <- names(make_safe_labels(x_var))
    
    df_agg <- df_plot %>%
      group_by(x_value) %>%
      summarise(mean_y = mean(y_value, na.rm = TRUE), .groups = "drop") %>%
      arrange(x_value) %>%
      mutate(cum_shap = cumsum(mean_y))
    
    max_point <- df_agg %>% slice(which.max(cum_shap))
    min_point <- df_agg %>% slice(which.min(cum_shap))
    
    y_min <- min(df_agg$cum_shap, na.rm = TRUE)
    y_max <- max(df_agg$cum_shap, na.rm = TRUE)
    y_pad <- 0.25 * (y_max - y_min)
    
    p <- geom_force_ribbon(df_agg) +
      geom_line(data = df_agg, aes(x = x_value, y = cum_shap),
                color = "gray40", linewidth = 0.4) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      geom_point(data = max_point, aes(x = x_value, y = cum_shap),
                 color = "#e41a1c", size = 3) +
      geom_point(data = min_point, aes(x = x_value, y = cum_shap),
                 color = "#377eb8", size = 3) +
      labs(
        title = paste0("Aggregated cumulative SHAP trend for ", display_name),
        x = x_display_label, y = "Cumulative Mean SHAP"
      ) +
      coord_cartesian(ylim = c(y_min - y_pad, y_max + y_pad)) +
      theme_minimal(base_size = 16)
    
    ggplotly(p)
  })
  
  # =====================================================
  # (3) Categorical Variables Plot
  # =====================================================
  output$catPlot <- renderPlotly({
    req(input$x_var_cat)
    cat_var <- input$x_var_cat
    if (!cat_var %in% names(test)) return(NULL)
    cat_values <- as.factor(test[[cat_var]])
    if (length(unique(cat_values)) > 20) return(NULL)
    
    cat_df <- df_filtered() %>%
      mutate(x_value = cat_values) %>%
      group_by(x_value) %>%
      summarise(mean_y = mean(y_value, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        Sign = ifelse(mean_y > 0, "Positive", "Negative"),
        mean_y_fmt = sprintf("%.6f", mean_y)
      ) %>%
      arrange(mean_y)
    
    display_name <- group_display_map$display[group_display_map$group == input$feature_display]
    if (length(display_name) == 0) display_name <- input$feature_display
    
    y_display_label <- names(make_safe_labels(cat_var))
    
    p <- ggplot(cat_df, aes(
      y = reorder(x_value, mean_y),
      x = mean_y,
      fill = Sign,
      text = paste0(
        display_name, ": ", x_value, "<br>",
        "Mean SHAP: ", mean_y_fmt, "<br>",
        "Sign: ", Sign
      )
    )) +
      geom_col(width = 0.65, alpha = 0.85) +
      scale_fill_manual(values = c("Positive" = "#ef8a62", "Negative" = "#67a9cf")) +
      geom_text(aes(label = mean_y_fmt),
                hjust = ifelse(cat_df$mean_y > 0, -0.15, 1.1),
                color = "black", size = 4, fontface = "bold") +
      labs(
        title = paste0("Mean SHAP by ", display_name),
        x = "Mean SHAP",
        y = y_display_label
      ) +
      theme_minimal(base_size = 16) +
      coord_cartesian(xlim = c(min(cat_df$mean_y)*1.3, max(cat_df$mean_y)*1.3)) +
      theme(
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 12)
      )
    
    ggplotly(p, tooltip = "text") %>%
      layout(hoverlabel = list(font = list(size = 13)))
  })
  
  # =====================================================
  # (4) Summary Importance Plot (with Excel display names)
  # =====================================================
  output$summaryPlot <- renderPlotly({
    mean_abs_shap <- round(colMeans(abs(shap_values_shap_avg)), 4)
    
    df_sum <- tibble(
      feature = names(mean_abs_shap),
      importance = mean_abs_shap
    ) %>%
      arrange(desc(importance)) %>%
      slice(1:25)
    
    rename_map<-c(
      "allergy_any_Yes_x_cot_ng_ml" = "Any Allergy: Yes × Cotinine (ng/mL)",
      "race_Non.Hispanic.White_x_LBXIGE" = "Race: Non-Hispanic White × Total IgE",
      "allergy_any_Yes_x_race_Non.Hispanic.White" = "Any Allergy × Race: Non-Hispanic White",
      "allergy_any_Yes_x_race_Non.Hispanic.Black" = "Any Allergy × Race: Non-Hispanic Black"
    )
    dummy_labels <- read_excel("dummy_namelist.xlsx")
    colnames(dummy_labels)[1:4] <- c("dummy", "group", "display_old", "rename")
    
    df_sum <- df_sum %>%
      left_join(dummy_list %>% select(dummy, rename), by = c("feature" = "dummy")) %>%
      mutate(feature_display = ifelse(!is.na(rename), rename, feature)) %>%
      mutate(feature_display = ifelse(
        feature %in% names(rename_map),
        rename_map[feature],
        feature_display
      ))
    
    
    p <- ggplot(df_sum, aes(x = reorder(feature_display, importance), y = importance)) +
      geom_col(fill = "#3182bd") +
      coord_flip() +
      labs(
        title = "Top 25 Most Important Features (mean |SHAP|)",
        x = NULL,
        y = "Mean |SHAP| Value"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.text.y = element_text(size = 12)
      )
    
    ggplotly(p)
  })
  
}

# =========================================================
# Run App
# =========================================================
shinyApp(ui, server)
