# ------------------------------------------------------------------------------
# Purpose: Provide reusable plotting helpers for manuscript result figures.
# Inputs: Prepared paper-result summary tables and plotting selections.
# Outputs: ggplot objects used in public-health and cost-effectiveness scripts.
# Dependencies: ggplot2/scales plus pipeline-standard column naming conventions.
# Run Stage: Function library called by s5-s8 result scripts.
# ------------------------------------------------------------------------------

#-helpful functions---------------------------------------------------------------------------------------------------------------
fill_colors <- c("5-doses" = "#94623d", "7-doses" = "#5489d0")

# Function: normalize_analysis_labels
# Purpose: Apply shared label recoding used across manuscript result scripts.
# Inputs: data frame with `coverage_assumption` and/or `seasonal_ages` columns.
# Outputs: Data frame with standardized label values for plotting/tables.
# Assumptions: Uses project-standard label mapping for coverage and seasonal ages.

normalize_analysis_labels <- function(data) {
  if ("coverage_assumption" %in% names(data)) {
    data <- data |>
      mutate(
        coverage_assumption = case_when(
          coverage_assumption == "S2-enhanced-communication" ~ "S2-moderate-dropout",
          coverage_assumption == "S3-sustained-dropout-10perc" ~ "S3-sustained-dropout",
          TRUE ~ coverage_assumption
        )
      )
  }

  if ("seasonal_ages" %in% names(data)) {
    data <- data |>
      mutate(
        seasonal_ages = case_when(
          seasonal_ages == "5-months-3-years" ~ "dose 1: 5m - 36m",
          TRUE ~ seasonal_ages
        )
      )
  }

  data
}

# Function: plot_events_averted
# Purpose: Plot events averted by PfPR and scenario facets.
# Inputs: Prepared summary data and plotting filters.
# Outputs: ggplot object for events-averted summaries.
# Assumptions: Assumes input labels and metric column names follow pipeline conventions.

plot_events_averted <- function(
    data,
    coverage = "S1-optimal-campaigns",
    SMC = 0,
    RTSS = c("Seasonal", "Hybrid"),
    dosing_assumption = c("5-doses", "7-doses"),
    vaccine_model = "model-1",
    seasonal_ages = c("5-months-1-year", "5-months-3-years", "age-based-primary"),
    y_var = "cases_averted_per_100000_FVCdose_dependent_discounted",
    fill_var = "dosing_assumption",
    facet_formula = str_to_title(seasonality) ~ coverage_assumption + label_var,
    fill_colors = c("5-doses" = "#94623d", "7-doses" = "#5489d0")) {
  ggplot(
    data %>%
      filter(
        SMC == !!SMC,
        RTSS %in% !!RTSS,
        dosing_assumption %in% !!dosing_assumption,
        vaccine_model %in% !!vaccine_model,
        # seasonal_ages %in% !!seasonal_ages,
        coverage_assumption %in% !!coverage
      ),
    aes(
      x = as.factor(pfpr * 100),
      y = .data[[y_var]],
      fill = .data[[fill_var]]
    )
  ) +
    geom_col(position = "dodge") +
    facet_grid(facet_formula) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = fill_colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8),
      legend.position = "top"
    ) +
    labs(
      y = glue::glue("{y_var}"),
      x = "PfPR 2-10 (%)",
      fill = " "
    )
}

# plot_events_averted_ab <- function(
#     data,
#     coverage = "S1-optimal-campaigns",
#     SMC = 0,
#     RTSS = c("Seasonal", "Hybrid"),
#     dosing_assumption = c("5-doses", "7-doses"),
#     vaccine_model = "model-1",
#     seasonal_ages = c("5-months-1-year", "5-months-3-years", "age-based-primary"),
#     y_var = "cases_averted_per_100000_FVCdose_dependent_discounted",
#     fill_var = "dosing_assumption",
#     facet_formula = str_to_title(seasonality) ~ coverage_assumption + label_var,
#     fill_colors = c("5-doses" = "#94623d", "7-doses" = "#5489d0")) {
#   ggplot(
#     data %>%
#       filter(
#         SMC == !!SMC,
#         RTSS %in% !!RTSS,
#         dosing_assumption %in% !!dosing_assumption,
#         vaccine_model %in% !!vaccine_model,
#         # seasonal_ages %in% !!seasonal_ages,
#         coverage_assumption %in% !!coverage
#       ),
#     aes(
#       x = as.factor(age_lower),
#       y = .data[[y_var]],
#       fill = .data[[fill_var]]
#     )
#   ) +
#     geom_col(position = "dodge") +
#     facet_grid(facet_formula, scales = "free") +
#     scale_y_continuous(labels = scales::comma) +
#     scale_fill_manual(values = fill_colors) +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(size = 8),
#       legend.position = "top"
#     ) +
#     labs(
#       y = glue::glue("{y_var}"),
#       x = "Age group (years)",
#       fill = " "
#     )
# }

# Function: plot_events_averted_ab
# Purpose: Plot age-banded events averted with per-panel totals.
# Inputs: Age-stratified data and plotting filters.
# Outputs: ggplot object for age-banded events-averted summaries.
# Assumptions: Assumes age_lower and selected y-variable columns are present.

plot_events_averted_ab <- function(
    data,
    coverage = "S1-optimal-campaigns",
    SMC = 0,
    RTSS = c("Seasonal", "Hybrid"),
    dosing_assumption = c("5-doses", "7-doses"),
    vaccine_model = "model-1",
    seasonal_ages = c("5-months-1-year", "5-months-3-years", "age-based-primary"),
    y_var = "cases_averted_per_100000_FVCdose_dependent_discounted",
    fill_var = "dosing_assumption",
    facet_formula = str_to_title(seasonality) ~ coverage_assumption + label_var,
    fill_colors = c("5-doses" = "#94623d", "7-doses" = "#5489d0")) {
  # Filter data
  plot_data <- data %>%
    filter(
      SMC == !!SMC,
      RTSS %in% !!RTSS,
      dosing_assumption %in% !!dosing_assumption,
      vaccine_model %in% !!vaccine_model,
      coverage_assumption %in% !!coverage
    )

  # Compute totals by facet and dosing assumption
  text_data <- plot_data %>%
    group_by(seasonality, coverage_assumption, pfpr, label_var, dosing_assumption) %>%
    summarise(
      total_events = sum(.data[[y_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      ea_label = glue("EA: {format(round(total_events,0), big.mark=',')}")
    )

  # Split data for separate geoms to color match
  text_data_5 <- text_data %>% filter(dosing_assumption == "5-doses")
  text_data_7 <- text_data %>% filter(dosing_assumption == "7-doses")

  ggplot(
    plot_data,
    aes(
      x = as.factor(age_lower),
      y = .data[[y_var]],
      fill = .data[[fill_var]]
    )
  ) +
    geom_col(position = "dodge") +
    facet_grid(facet_formula, scales = "free") +
    # Add text for 5-doses
    geom_text(
      data = text_data_5,
      aes(
        x = Inf, y = Inf, label = ea_label
      ),
      hjust = 1.1, vjust = 1.4,
      color = fill_colors["5-doses"],
      size = 3.2,
      inherit.aes = FALSE
    ) +
    # Add text for 7-doses
    geom_text(
      data = text_data_7,
      aes(
        x = Inf, y = Inf, label = ea_label
      ),
      hjust = 1.1, vjust = 2.6,
      color = fill_colors["7-doses"],
      size = 3.2,
      inherit.aes = FALSE
    ) +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = fill_colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8),
      legend.position = "top"
    ) +
    labs(
      y = glue("{y_var}"),
      x = "Age group (years)",
      fill = " "
    )
}



# Function: plot_incremental_impact
# Purpose: Plot incremental impact comparisons between dose schedules.
# Inputs: Incremental summary data and outcome variable selection.
# Outputs: ggplot object for incremental impact outputs.
# Assumptions: Assumes precomputed incremental metrics in input data.

plot_incremental_impact <- function(
    data,
    y_var,
    smc_var = 0,
    plot_title,
    filename) {
  library(ggplot2)
  library(scales)

  p <- ggplot(
    data %>%
      filter(
        SMC == smc_var,
        dosing_assumption_to == "7-doses",
        dosing_assumption_from == "5-doses",
        pfpr > 0.01
      ),
    aes(
      x = as.factor(pfpr * 100),
      y = .data[[y_var]],
      color = RTSS,
      shape = vaccine_model
    )
  ) +
    geom_point(position = position_jitter(width = 0.01)) +
    geom_smooth(
      aes(
        group = interaction(vaccine_model, label_var, seasonality, coverage_assumption),
        linetype = vaccine_model
      ),
      se = FALSE,
      size = 1
    ) +
    scale_color_manual(values = c(
      "Hybrid" = "#550527",
      "Seasonal" = "#688E26"
    )) +
    facet_grid(seasonality ~ coverage_assumption) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme_minimal() +
    labs(
      x = "PfPR 2-10 (%)",
      y = "Incremental impact of 7-dose vs 5-dose schedule (%)",
      title = plot_title,
      col = "Delivery",
      shape = "Vaccine model",
      linetype = "Vaccine model"
    ) +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      legend.title.align = 0.5
    ) +
    guides(
      color = guide_legend(nrow = 3),
      shape = guide_legend(nrow = 2),
      linetype = guide_legend(nrow = 2)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1)

  return(p)
}

# Function: plot_incremental_impact_by_outcome
# Purpose: Plot incremental impact faceted by selected outcomes.
# Inputs: Incremental data and outcome/label filters.
# Outputs: ggplot object faceted by outcome and setting.
# Assumptions: Assumes standardized outcome columns in input frame.

plot_incremental_impact_by_outcome <- function(
    data,
    smc_var = 0,
    plot_title,
    filename = NULL) {
  library(ggplot2)
  library(scales)

  p <- ggplot(
    data %>%
      filter(
        SMC == smc_var,
        dosing_assumption_to == "7-doses",
        dosing_assumption_from == "5-doses"
      ),
    aes(
      x = as.factor(pfpr * 100),
      y = incremental_impact,
      color = RTSS
    )
  ) +
    geom_point(position = position_jitter(width = 0.03)) +
    scale_color_manual(values = c(
      "Hybrid" = "#550527",
      "Seasonal" = "#688E26"
    )) +
    facet_grid(str_to_title(seasonality) ~ outcome) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme_minimal() +
    labs(
      x = "PfPR 2-10 (%)",
      y = "% addidional outcomes averted",
      title = plot_title,
      col = "Delivery",
      shape = "Coverage"
    ) +
    theme(
      legend.position = "right" # ,
      # legend.box = "horizontal",
      # legend.title.align = 0.5
    ) +
    # guides(
    #   color = guide_legend(nrow = 3),
    #   shape = guide_legend(nrow = 2)
    # ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1)
  # scale_y_continuous(limits=c(-30, 105))

  if (!is.null(filename)) {
    ggsave(
      plot = p,
      filename = filename,
      width = 8,
      height = 5,
      dpi = 600
    )
  }

  return(p)
}

# plot_cost_effectiveness <- function(
#     data,
#     title_text,
#     y_var,
#     y_lab,
#     smc_value = 0,
#     cost_per_dose_values = c(4),
#     dosing_assumptions = c("5-doses", "7-doses")) {
#   y_var_enquo <- enquo(y_var)
#
#   ggplot(
#     data %>%
#       filter(
#         SMC == smc_value,
#         rtss_cost_per_dose %in% cost_per_dose_values,
#         dosing_assumption %in% dosing_assumptions
#       ),
#     aes(
#       x = as.factor(pfpr * 100),
#       y = !!y_var_enquo,
#       fill = dosing_assumption
#     )
#   ) +
#     geom_boxplot(
#       outlier.size = 0.5, position = position_dodge(width = 0.8),
#       alpha = 0.8
#     ) +
#     facet_grid(str_to_title(seasonality) ~ vaccine_model + paste0("$", rtss_cost_per_dose, " per dose")) +
#     scale_y_log10(labels = comma) +
#     scale_fill_manual(values = fill_colors) +
#     labs(
#       title = title_text,
#       x = "PfPR 2-10 (%)",
#       y = y_lab,
#       fill = "Schedule"
#     ) +
#     theme_minimal() +
#     theme(legend.position = "top")
# }

# Function: plot_cost_effectiveness
# Purpose: Plot cost-effectiveness relationship for selected scenarios.
# Inputs: Cost-effectiveness summary data and plotting options.
# Outputs: ggplot object for CE scatter/line displays.
# Assumptions: Assumes cost and effectiveness columns already computed.

plot_cost_effectiveness <- function(
    data,
    title_text,
    y_var,
    y_lab,
    colour_var,
    fill_var,
    group_var,
    facet_row,
    facet_col,
    smc_value = 0,
    cost_per_dose_values = c(4),
    dosing_assumptions = c("5-doses", "7-doses"),
    fill_colours) {
  ggplot(
    data %>%
      filter(
        SMC == smc_value,
        rtss_cost_per_dose %in% cost_per_dose_values,
        dosing_assumption %in% dosing_assumptions
      ),
    aes(
      x = as.factor(pfpr * 100),
      y = {{ y_var }},
      color = {{ colour_var }},
      fill = {{ fill_var }},
      group = {{ group_var }},
      linetype = vaccine_model,
      shape = coverage_assumption
    )
  ) +
    geom_line() +
    geom_point() +
    facet_grid(
      rows = vars(!!sym(facet_row)),
      cols = vars(!!sym(facet_col))
    ) +
    scale_y_log10(labels = scales::comma) +
    labs(
      title = title_text,
      x = "PfPR 2-10 (%)",
      y = y_lab,
      color = rlang::as_label(enquo(colour_var)),
      fill = rlang::as_label(enquo(fill_var))
    ) +
    theme_minimal() +
    theme(legend.position = "top") +
    scale_fill_manual(values = fill_colours) +
    scale_color_manual(values = fill_colours)
}

# Function: plot_diff_ICER_points
# Purpose: Plot pairwise ICER differences across settings.
# Inputs: ICER summary frame and scenario selectors.
# Outputs: ggplot object showing ICER deltas.
# Assumptions: Assumes precomputed ICER fields and consistent scenario labels.

plot_diff_ICER_points <- function(
    data,
    y_var,
    pfpr_threshold = 0.1,
    smc_value = 0,
    shape_var = coverage_assumption,
    color_var = label_var,
    alpha_var = paste0("$", rtss_cost_per_dose),
    facet_row = str_to_title(seasonality),
    facet_col = vaccine_model,
    facet_scale = "free_y",
    plot_title = "Difference in ICER (7-dose minus 5-dose) across PfPR levels",
    y_lab = "Incremental Cost-Effectiveness") {
  ggplot(
    data |> filter(pfpr >= pfpr_threshold, SMC == smc_value),
    aes(
      x = as.factor(pfpr * 100),
      y = {{ y_var }},
      shape = {{ shape_var }},
      color = {{ color_var }},
      alpha = {{ alpha_var }}
    )
  ) +
    geom_point(
      position = position_jitter(width = 0.2, height = 0),
      size = 2
    ) +
    facet_grid(
      rows = vars({{ facet_row }}),
      cols = vars({{ facet_col }}),
      scales = facet_scale
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = plot_title,
      x = "PfPR 2-10 (%)",
      y = y_lab,
      col = "",
      shape = "",
      alpha = "",
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(size = 8)
    )
}
