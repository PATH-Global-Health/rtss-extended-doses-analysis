# ------------------------------------------------------------------------------
# Purpose: Run CE presentation and table outputs using 5-dose schedule as baseline comparator.
# Inputs: Precomputed 5-dose-baseline cost output dataset.
# Outputs: 5-dose-baseline CE figure and summary table outputs.
# Dependencies: data_and_libraries.R and patchwork/table helpers.
# Run Stage: Pipeline step 8: 5-dose-baseline result outputs.
# Notes: Header wording corrected; no logic or output naming changes.
# ------------------------------------------------------------------------------

# source
source("./02.scripts/data_and_libraries.R")
source("./02.scripts/functions/f9-plotting-results.R")
library(patchwork)

#-READ IN AND PREP--------------------------------------------------------------------------------------------------------------------
cost_output_vaccine_5dose_baseline <- readRDS("./03.outputs/paper-results/cost_outputs_5dose_baseline.rds")

## remove the seasonal age groups that are not used in the analysis-------------
cost_output_vaccine_5dose_baseline <- cost_output_vaccine_5dose_baseline |> filter(seasonal_ages != "5-months-1-year")

## rename variables------------------------------------------------------------
cost_output_vaccine_5dose_baseline <- normalize_analysis_labels(cost_output_vaccine_5dose_baseline)

colors_labels <-
  c(
    "5-doses" = "#94623d",
    "7-doses" = "#5489d0"
  )

cost_output_vaccine_5dose_baseline <-
  cost_output_vaccine_5dose_baseline |>
  mutate(
    label = paste0(dosing_assumption, " | ", RTSS),
    seasonality = str_to_title(seasonality),
    RTSS = paste0(RTSS, " delivery"),
    pfpr_label = paste0("PfPR 2-10: ", pfpr * 100, "%"),
    rtss_cost_per_dose = round(as.numeric(rtss_cost_per_dose), 2)
  )

cost_output_vaccine_5dose_baseline$pfpr_label <- factor(
  cost_output_vaccine_5dose_baseline$pfpr_label,
  levels = c(
    "PfPR 2-10: 1%", "PfPR 2-10: 3%", "PfPR 2-10: 5%", "PfPR 2-10: 10%", "PfPR 2-10: 15%",
    "PfPR 2-10: 20%", "PfPR 2-10: 25%", "PfPR 2-10: 35%", "PfPR 2-10: 45%", "PfPR 2-10: 55%",
    "PfPR 2-10: 65%"
  )
)

## -Standardise the object name used below--------------------------------------
cost_outputs_5base <- cost_output_vaccine_5dose_baseline

#-SUMMARISE ICERs OVER UNCERTAINTY-----------------------------------------------------------------------------------------------
ce_summary_5base <-
  cost_outputs_5base |>
  group_by(
    RTSS,
    label,
    pfpr,
    pfpr_label,
    seasonality,
    coverage_assumption,
    SMC,
    rtss_cost_per_dose
  ) |>
  summarise(
    across(
      all_of(c(
        "cost_effectiveness_daly_discounted",
        "cost_effectiveness_case_discounted",
        "incremental_cost_discounted",
        "dalys_averted_discounted",
        "cases_averted_discounted"
      )),
      list(
        med = ~ median(.x, na.rm = TRUE),
        l95 = ~ quantile(.x, 0.025, na.rm = TRUE),
        u95 = ~ quantile(.x, 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  mutate(
    seasonality = str_to_title(seasonality),
    cost_per_dose = case_when(
      dplyr::near(rtss_cost_per_dose, 5.19) ~ "$4",
      dplyr::near(rtss_cost_per_dose, 2.78) ~ "$2",
      dplyr::near(rtss_cost_per_dose, 6.39) ~ "$5",
      TRUE ~ NA_character_
    )
  )

write.csv(ce_summary_5base, "5-dose-as-baseline.csv")


#-ICER by PfPR plot----------------------------------------------------------------------------
# Function: make_ce_plot_5base
# Purpose: Create CE comparison plots using 5-dose schedule as baseline.
# Inputs: summary_df and plotting controls.
# Outputs: ggplot object for 5-dose-baseline CE analysis.
# Assumptions: Assumes baseline-adjusted ICER metrics are precomputed.

make_ce_plot_5base <- function(summary_df,
                               outcome = c("cases", "dalys"),
                               cost_label = "$4",
                               smc_filter = 0,
                               coverage_filter = "S2-moderate-dropout") {
  outcome <- match.arg(outcome)

  if (outcome == "cases") {
    y_med_col <- "cost_effectiveness_case_discounted_med"
    y_l95_col <- "cost_effectiveness_case_discounted_l95"
    y_u95_col <- "cost_effectiveness_case_discounted_u95"
    y_label <- "Incremental ICER per case\n averted vs 5-dose (USD$)"
    subtitle <- "Uncomplicated cases"
  } else {
    y_med_col <- "cost_effectiveness_daly_discounted_med"
    y_l95_col <- "cost_effectiveness_daly_discounted_l95"
    y_u95_col <- "cost_effectiveness_daly_discounted_u95"
    y_label <- "Incremental ICER per DALY\n averted vs 5-dose (USD$)"
    subtitle <- "DALYs"
  }

  plot_data <-
    summary_df |>
    filter(
      SMC == smc_filter,
      pfpr > 0.01,
      coverage_assumption == coverage_filter,
      cost_per_dose == cost_label
    )

  ggplot(plot_data) +
    geom_ribbon(
      aes(
        x = pfpr * 100,
        ymin = .data[[y_l95_col]],
        ymax = .data[[y_u95_col]],
        group = RTSS,
        col = RTSS,
        fill = RTSS
      ),
      alpha = 0.25,
      colour = NA
    ) +
    geom_line(
      aes(
        x = pfpr * 100,
        y = .data[[y_med_col]],
        group = RTSS,
        col = RTSS,
        fill = RTSS
      ),
      linewidth = 0.5
    ) +
    geom_point(
      aes(
        x = pfpr * 100,
        y = .data[[y_med_col]],
        col = RTSS,
        fill = RTSS
      ),
      size = 0.7
    ) +
    facet_grid(. ~ seasonality, scales = "free_y") +
    scale_x_continuous(breaks = c(3, 5, 10, 15, 20, 25, 35, 45, 55, 65)) +
    theme_minimal(base_size = 10) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(
      x = "PfPR 2–10 (%)",
      y = y_label,
      subtitle = subtitle
    ) +
    scale_color_manual(values = c(
      "#C75687",
      "#BCBD8B"
    )) +
    scale_fill_manual(values = c(
      "#C75687",
      "#BCBD8B"
    ))
}

ce5_cases <- make_ce_plot_5base(ce_summary_5base, outcome = "cases", cost_label = "$4", smc_filter = 0)
ce5_dalys <- make_ce_plot_5base(ce_summary_5base, outcome = "dalys", cost_label = "$4", smc_filter = 0)

(ce5_cases / ce5_dalys) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "./05.plots/publication-plots/appendix/icer-7v5-plane.png",
  width = 8, height = 8, dpi = 600
)

#-SUMMARY TABLE---------------------------------------------------------------------------------------------
# Function: format_icer
# Purpose: Format ICER numeric values for human-readable table output.
# Inputs: x numeric vector.
# Outputs: Character vector of formatted ICER values.
# Assumptions: Assumes values are in expected cost-effectiveness units.

format_icer <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(
      abs(x) < 1,
      sprintf("%.2f", x), # values between -1 and 1 → 2 decimals
      sprintf("%.0f", x) # everything else → integer
    )
  )
}

# Function: make_ce_table
# Purpose: Generate a publication-ready ICER summary table.
# Inputs: summary_df and table formatting options.
# Outputs: Formatted table object/data frame for export.
# Assumptions: Assumes standardized scenario labels and columns.

make_ce_table <- function(summary_df,
                          cost_label = c("$2", "$4", "$5"),
                          smc_filter = 0,
                          coverage_filter = "S2-moderate-dropout") {
  summary_df |>
    filter(
      SMC %in% smc_filter,
      pfpr > 0.01,
      coverage_assumption %in% coverage_filter,
      cost_per_dose %in% cost_label
    ) |>
    mutate(
      seasonality = stringr::str_to_title(seasonality),
      pfpr_label  = paste0(pfpr * 100, "%"),
      label       = paste0("7-dose | ", RTSS)
    ) |>
    # format med / l95 / u95 as strings with the custom rules
    mutate(
      case_med_str = format_icer(cost_effectiveness_case_discounted_med),
      case_l95_str = format_icer(cost_effectiveness_case_discounted_l95),
      case_u95_str = format_icer(cost_effectiveness_case_discounted_u95),
      daly_med_str = format_icer(cost_effectiveness_daly_discounted_med),
      daly_l95_str = format_icer(cost_effectiveness_daly_discounted_l95),
      daly_u95_str = format_icer(cost_effectiveness_daly_discounted_u95),
      `ICER per case averted (USD$)` = sprintf(
        "%s (%s–%s)",
        case_med_str, case_l95_str, case_u95_str
      ),
      `ICER per DALY averted (USD$)` = sprintf(
        "%s (%s–%s)",
        daly_med_str, daly_l95_str, daly_u95_str
      )
    ) |>
    # keep key descriptors + pfpr + formatted strings
    select(
      label,
      seasonality,
      cost_per_dose,
      pfpr_label,
      SMC,
      coverage_assumption,
      `ICER per case averted (USD$)`,
      `ICER per DALY averted (USD$)`
    ) |>
    # one row per outcome per (label, seasonality, cost_per_dose)
    tidyr::pivot_longer(
      cols = c(
        `ICER per case averted (USD$)`,
        `ICER per DALY averted (USD$)`
      ),
      names_to = "outcome",
      values_to = "value"
    ) |>
    # wide on PfPR so each column is a PfPR level
    tidyr::pivot_wider(
      names_from  = pfpr_label,
      values_from = value
    ) |>
    arrange(seasonality, label, outcome)
}

# Table with all vaccine delivery costs SMC and coverage
ce_table_all <- make_ce_table(ce_summary_5base,
  cost_label = c("$2", "$4", "$5"),
  smc_filter = c(0, 0.75),
  coverage_filter = c("S1-optimal-campaigns", "S2-moderate-dropout", "S3-sustained-dropout")
)

writexl::write_xlsx(list(`ce-table-5-dose-bl` = ce_table_all),
  path = "./05.plots/publication-plots/tables/cost-effect-table-revisions-supp.xlsx"
)
