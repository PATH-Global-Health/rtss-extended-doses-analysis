# ------------------------------------------------------------------------------
# Purpose: Generate cost-effectiveness figures, ICER summaries, and frontier analyses.
# Inputs: Paper-results cost output dataset from step 4.
# Outputs: Main-text and appendix CE figures plus ICER summary tables.
# Dependencies: data_and_libraries.R and CE/statistical helper packages.
# Run Stage: Pipeline step 6: manuscript cost-effectiveness results.
# ------------------------------------------------------------------------------

source("./02.scripts/data_and_libraries.R")
source("./02.scripts/functions/f9-plotting-results.R")
library(patchwork)

#-READ IN AND PREP--------------------------------------------------------------------------------------------------------------------
cost_outputs <- readRDS("./03.outputs/paper-results/cost_outputs.rds")

## remove the seasonal age groups that are not used in the analysis-------------
cost_outputs <- cost_outputs |> filter(seasonal_ages != "5-months-1-year")

## rename variables-------------------------------------------------------------
cost_outputs <- normalize_analysis_labels(cost_outputs)

colors_labels <-
  c(
    "5-doses" = "#94623d",
    "7-doses" = "#5489d0"
  )

cost_outputs <-
  cost_outputs |>
  mutate(
    label = paste0(dosing_assumption, " | ", RTSS),
    seasonality = str_to_title(seasonality),
    RTSS = paste0(RTSS, " delivery"),
    pfpr_label = paste0("PfPR 2-10: ", pfpr * 100, "%"),
    rtss_cost_per_dose = round(as.numeric(rtss_cost_per_dose), 2)
  )

cost_outputs$pfpr_label <- factor(
  cost_outputs$pfpr_label,
  levels = c(
    "PfPR 2-10: 1%", "PfPR 2-10: 3%", "PfPR 2-10: 5%", "PfPR 2-10: 10%", "PfPR 2-10: 15%",
    "PfPR 2-10: 20%", "PfPR 2-10: 25%", "PfPR 2-10: 35%", "PfPR 2-10: 45%", "PfPR 2-10: 55%",
    "PfPR 2-10: 65%"
  )
)

#-INCREMENTAL COST VS INCREMENTAL IMPACT-------------------------------------------------------------------------------------

## cost vs impact for cases and DALYS scatter plot
## accross delivery, pfpr
# Function: make_cost_vs_impact_plot
# Purpose: Create cost-versus-impact comparison plots for selected outcomes.
# Inputs: Prepared summary data and plotting settings.
# Outputs: ggplot object for manuscript or appendix use.
# Assumptions: Assumes expected summary columns and factor labels are present.

make_cost_vs_impact_plot <- function(data,
                                     rtss_type,
                                     outcome = c("cases", "dalys"),
                                     colors_labels,
                                     cost_per_dose = 5.19) {
  outcome <- match.arg(outcome)

  # Choose x variable and label based on outcome type
  x_col <- if (outcome == "cases") "cases_averted_discounted" else "dalys_averted_discounted"
  x_label <- if (outcome == "cases") {
    "Uncomplicated cases averted per 100,000 population (thousands)"
  } else {
    "DALYs averted per 100,000 population (thousands)"
  }

  df_plot <-
    data |>
    filter(
      rtss_cost_per_dose == cost_per_dose,
      RTSS == rtss_type
    ) |>
    mutate(
      x_val = ((.data[[x_col]] / (n_u5 + n_o5)) * 100000) / 1000,
      y_val = ((incremental_cost_discounted / (n_u5 + n_o5)) * 100000) / 1000
    )

  ggplot(df_plot) +
    geom_point(
      aes(
        x   = x_val,
        y   = y_val,
        col = dosing_assumption
      ),
      alpha = 0.2
    ) +
    facet_wrap(. ~ pfpr_label) +
    theme_minimal() +
    scale_y_continuous(limits = c(0, NA)) +
    scale_color_manual(values = colors_labels) +
    labs(
      x        = x_label,
      y        = "Incremental cost per 100,000 population (thousands)",
      col      = "Dosing",
      subtitle = rtss_type
    ) +
    theme(axis.title = element_text(size = 8))
}

# individual plots
cost_vs_impact_seas_cases <-
  make_cost_vs_impact_plot(
    data          = cost_outputs,
    rtss_type     = "Seasonal delivery",
    outcome       = "cases",
    colors_labels = colors_labels
  )

cost_vs_impact_seas_dalys <-
  make_cost_vs_impact_plot(
    data          = cost_outputs,
    rtss_type     = "Seasonal delivery",
    outcome       = "dalys",
    colors_labels = colors_labels
  )

cost_vs_impact_hyb_cases <-
  make_cost_vs_impact_plot(
    data          = cost_outputs,
    rtss_type     = "Hybrid delivery",
    outcome       = "cases",
    colors_labels = colors_labels
  )

cost_vs_impact_hyb_dalys <-
  make_cost_vs_impact_plot(
    data          = cost_outputs,
    rtss_type     = "Hybrid delivery",
    outcome       = "dalys",
    colors_labels = colors_labels
  )

# combined plot
(cost_vs_impact_seas_cases + cost_vs_impact_seas_dalys) /
  (cost_vs_impact_hyb_cases + cost_vs_impact_hyb_dalys) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("05.plots/publication-plots/appendix/cost-vs-impact.png", width = 10, height = 8, dpi = 600)

#-ICER---------------------------------------------------------------------------------------------------------------

## -summarise data median and 95% range------------
ce_summary <-
  cost_outputs |>
  group_by(
    RTSS,
    dosing_assumption,
    pfpr,
    seasonality,
    coverage_assumption,
    SMC,
    rtss_cost_per_dose
  ) |>
  summarise(
    across(
      all_of(c(
        "cost_effectiveness_daly_discounted",
        "cost_effectiveness_case_discounted"
      )),
      list(
        med   = ~ median(.x, na.rm = TRUE),
        l95   = ~ quantile(.x, 0.025, na.rm = TRUE),
        u95   = ~ quantile(.x, 0.975, na.rm = TRUE)
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

# plotting function for ICER plot
# Function: make_ce_plot_main
# Purpose: Build main cost-effectiveness panel plots from summarized ICER inputs.
# Inputs: summary_df and selected outcome/format controls.
# Outputs: ggplot object for main CE figure output.
# Assumptions: Assumes summary data already aggregated by scenario and setting.

make_ce_plot_main <- function(summary_df,
                              outcome = c("cases", "dalys"),
                              cost_label = "$4",
                              smc_filter = 0,
                              coverage_filter = "S2-moderate-dropout") {
  outcome <- match.arg(outcome)

  # Map to correct columns / labels
  if (outcome == "cases") {
    y_med_col <- "cost_effectiveness_case_discounted_med"
    y_l95_col <- "cost_effectiveness_case_discounted_l95"
    y_u95_col <- "cost_effectiveness_case_discounted_u95"
    y_label <- "ICER per case averted (USD$)"
    subtitle <- "Uncomplicated cases"
  } else {
    y_med_col <- "cost_effectiveness_daly_discounted_med"
    y_l95_col <- "cost_effectiveness_daly_discounted_l95"
    y_u95_col <- "cost_effectiveness_daly_discounted_u95"
    y_label <- "ICER per DALY averted (USD$)"
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

  ggplot(
    plot_data,
    aes(
      x = (pfpr * 100),
      y = .data[[y_med_col]],
      color = dosing_assumption,
      fill = dosing_assumption,
      linetype = RTSS,
      group = interaction(dosing_assumption, RTSS)
    )
  ) +
    # 95% CI ribbon FIRST so it sits under the lines
    geom_ribbon(
      aes(
        ymin  = .data[[y_l95_col]],
        ymax  = .data[[y_u95_col]]
      ),
      alpha = 0.3,
      linetype = 0,
      colour = NA
    ) +
    # Median line
    geom_line(linewidth = 0.4) +
    # Median points
    geom_point(size = 0.6) +
    facet_grid(. ~ seasonality) +
    scale_x_continuous(breaks = c(3, 5, 10, 15, 20, 25, 35, 45, 55, 65)) +
    scale_y_continuous(limits=c(0,NA))+
    scale_color_manual(values = colors_labels) +
    scale_fill_manual(values = colors_labels) +
    scale_linetype_manual(values = c(
      "Hybrid delivery" = "solid",
      "Seasonal delivery" = "dashed"
    )) +
    theme_minimal(8) +
    labs(
      x = "PfPR 2–10 (%)",
      y = y_label,
      subtitle = subtitle,
      color = "Dosing Schedule",
      fill = "Dosing Schedule",
      linetype = "RTSS"
    )
}

# Create the two panels
ce_1 <- make_ce_plot_main(ce_summary, outcome = "cases", cost_label = "$4")
ce_2 <- make_ce_plot_main(ce_summary, outcome = "dalys", cost_label = "$4")

ce_results_main_text_plot <- ce_1 / ce_2 +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("./05.plots/publication-plots/main-text/ce-results-main-text.png", plot = ce_results_main_text_plot, width = 10, height = 7, dpi = 600, bg = "white")
ggsave("./05.plots/publication-plots/main-text/ce-results-main-text.pdf", plot = ce_results_main_text_plot, width = 10, height = 7, dpi = 600, bg = "white")

# Make this in table summary
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
      label       = paste0(dosing_assumption, " | ", RTSS)
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
ce_table_all <- make_ce_table(ce_summary,
  cost_label = c("$2", "$4", "$5"),
  smc_filter = c(0, 0.75),
  coverage_filter = c("S1-optimal-campaigns", "S2-moderate-dropout", "S3-sustained-dropout")
)

writexl::write_xlsx(list(`ce-table1` = ce_table_all),
  path = "./05.plots/publication-plots/tables/cost-effect-table-revisions.xlsx"
)

# how to plot for SMC and coverage too? Large facet plot for supplement
supp_costs_icer <-
  ggplot(
    ce_summary |> filter(rtss_cost_per_dose == 5.19)
  ) +
  geom_point(
    aes(
      x = pfpr * 100,
      y = cost_effectiveness_case_discounted_med,
      col = dosing_assumption,
      group = interaction(dosing_assumption, as.factor(SMC), coverage_assumption, RTSS)
    )
  ) +
  geom_line(
    aes(
      x = pfpr * 100,
      y = cost_effectiveness_case_discounted_med,
      col = dosing_assumption,
      linetype = as.factor(SMC),
      group = interaction(dosing_assumption, as.factor(SMC), coverage_assumption, RTSS)
    )
  ) +
  facet_grid(seasonality ~ RTSS + coverage_assumption) +
  scale_color_manual(values = colors_labels) +
  scale_fill_manual(values = colors_labels) +
  scale_x_continuous(breaks = c(3, 5, 10, 15, 20, 25, 35, 45, 55, 65)) +
  theme_minimal(8) +
  labs(
    x = "PfPR 2–10 (%)",
    y = "ICER per case averted (USD$)",
    color = "Dosing Schedule",
    linetype = "SMC coverage"
  )

supp_dalys_icer <-
  ggplot(
    ce_summary |> filter(rtss_cost_per_dose == 5.19)
  ) +
  geom_point(
    aes(
      x = pfpr * 100,
      y = cost_effectiveness_daly_discounted_med,
      col = dosing_assumption,
      group = interaction(dosing_assumption, as.factor(SMC), coverage_assumption, RTSS)
    )
  ) +
  geom_line(
    aes(
      x = pfpr * 100,
      y = cost_effectiveness_daly_discounted_med,
      col = dosing_assumption,
      linetype = as.factor(SMC),
      group = interaction(dosing_assumption, as.factor(SMC), coverage_assumption, RTSS)
    )
  ) +
  facet_grid(seasonality ~ RTSS + coverage_assumption) +
  scale_color_manual(values = colors_labels) +
  scale_fill_manual(values = colors_labels) +
  scale_x_continuous(breaks = c(3, 5, 10, 15, 20, 25, 35, 45, 55, 65)) +
  theme_minimal(8) +
  labs(
    x = "PfPR 2–10 (%)",
    y = "ICER per DALY averted (USD$)",
    color = "Dosing Schedule",
    linetype = "SMC coverage"
  )

supp_costs_icer / supp_dalys_icer +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("./05.plots/publication-plots/appendix/ce-results-supp-smc-coverage.png", width = 10, height = 8, dpi = 600)


#-EXTENDED DOMINANCE ANALYSIS------------------------------------------------------------

# median summaries
# Function: summarise_frontier_inputs
# Purpose: Prepare frontier-analysis inputs grouped by setting and strategy.
# Inputs: cost_outputs plus grouping/selection settings.
# Outputs: Summarized frontier input table.
# Assumptions: Assumes dominance-analysis columns exist in input data.

summarise_frontier_inputs <- function(cost_outputs,
                                      cost_per_dose = 5.19,
                                      smc_filter = 0,
                                      coverage_filter = "S2-moderate-dropout",
                                      effect = c("dalys", "cases")) {
  effect <- match.arg(effect)

  # choose the underlying effect column
  effect_col <- if (effect == "dalys") {
    "dalys_averted_discounted"
  } else {
    "cases_averted_discounted"
  }

  cost_outputs %>%
    filter(
      SMC == smc_filter,
      coverage_assumption == coverage_filter,
      rtss_cost_per_dose == cost_per_dose,
      pfpr > 0.01
    ) %>%
    group_by(
      pfpr,
      seasonality,
      RTSS,
      dosing_assumption
    ) %>%
    summarise(
      effect_med = median(
        (.data[[effect_col]] / (n_u5 + n_o5)) * 1e5 / 1e3,
        na.rm = TRUE
      ),
      cost_med = median(
        (incremental_cost_discounted / (n_u5 + n_o5)) * 1e5 / 1e3,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      seasonality = str_to_title(seasonality),
      pfpr_label  = paste0(pfpr * 100, "%"),
      strategy    = paste(dosing_assumption, RTSS, sep = " | "),
      effect_type = effect
    )
}


# Helper 2: Run dampack::calculate_icers() within each setting
#           Returns dominance status + ICERs for each strategy
# Function: run_icers_by_setting
# Purpose: Run ICER calculations for each epidemiologic setting subset.
# Inputs: frontier_summary table.
# Outputs: Combined ICER results table by setting.
# Assumptions: Assumes dampack-style ICER functions are available.

run_icers_by_setting <- function(frontier_summary) {
  frontier_summary %>%
    group_by(pfpr, pfpr_label, seasonality) %>%
    group_modify(~ {
      df <- .x

      icer_obj <- dampack::calculate_icers(
        cost       = df$cost_med,
        effect     = df$effect_med,
        strategies = df$strategy
      )

      as.data.frame(icer_obj) %>%
        rename(
          strategy     = Strategy,
          cost_med     = Cost,
          effect_med   = Effect,
          delta_cost   = Inc_Cost,
          delta_effect = Inc_Effect,
          icer         = ICER,
          status       = Status
        ) %>%
        left_join(
          df %>%
            select(strategy, RTSS, dosing_assumption, effect_type) %>%
            distinct(),
          by = "strategy"
        )
    }) %>%
    ungroup() %>%
    mutate(
      frontier_flag = if_else(status == "ND",
        "Efficient Frontier", "Weakly Dominated"
      ),
      frontier_flag = factor(
        frontier_flag,
        levels = c("Efficient Frontier", "Weakly Dominated")
      )
    )
}


# Helper 3: Plot ICERs in a facet grid
#           label = "frontier", "all", "none", "weak"
# Function: plot_icers_faceted
# Purpose: Plot ICER outcomes faceted across seasonality/PfPR settings.
# Inputs: icer_df and plotting options.
# Outputs: Faceted ggplot of ICER and dominance results.
# Assumptions: Assumes ICER output columns and labels are standardized.

plot_icers_faceted <- function(icer_df,
                               label = c("frontier", "all", "none", "weak"),
                               effect_label = "DALYs averted per 100,000 population (thousands)") {
  label <- match.arg(label)

  # Which points to label?
  label_data <- switch(label,
    frontier = dplyr::filter(icer_df, frontier_flag == "Efficient Frontier"),
    all      = icer_df,
    weak     = dplyr::filter(icer_df, frontier_flag == "Weakly Dominated"),
    none     = NULL
  )

  # Data for connecting the frontier line within each facet
  frontier_lines <-
    icer_df %>%
    dplyr::filter(frontier_flag == "Efficient Frontier") %>%
    dplyr::arrange(pfpr, seasonality, effect_med)

  ggplot(icer_df, aes(x = effect_med, y = cost_med)) +
    # Efficient frontier connecting line
    geom_path(
      data = frontier_lines,
      aes(group = interaction(pfpr, seasonality)),
      lineend = "round"
    ) +
    # All points
    geom_point(aes(col = frontier_flag)) +

    # Optional labels
    {
      if (!is.null(label_data)) {
        ggrepel::geom_text_repel(
          data = label_data,
          aes(label = strategy),
          size = 2,
          box.padding = unit(0.3, "lines"),
          point.padding = unit(0.5, "lines"),
          force_pull = 0.2,
          force = 2,
          min.segment.length = 3,
          label.size = 0,
          fill = NA
        )
      } else {
        NULL
      }
    } +
    facet_wrap(seasonality ~ pfpr, scales = "free", ncol = 10) +
    theme_minimal() +
    labs(
      x   = effect_label,
      y   = "Incremental cost per 100,000 population (thousands)",
      col = " "
    ) +
    scale_color_manual(values = c("black", "#db652a"))
}


# Run the full pipeline
#   - assumes cost_outputs already exists in the workspace


# DALY-based frontier
frontier_dalys <-
  summarise_frontier_inputs(
    cost_outputs,
    cost_per_dose   = 5.19, # change if needed
    smc_filter      = 0,
    coverage_filter = "S2-moderate-dropout",
    effect          = "dalys"
  )

icers_dalys <- run_icers_by_setting(frontier_dalys)

p_frontier_dalys <- plot_icers_faceted(
  icers_dalys,
  label        = "all", # "frontier", "all", "none", or "weak"
  effect_label = "DALYs averted per 100,000 population (thousands)"
)

# Case-based frontier
frontier_cases <-
  summarise_frontier_inputs(
    cost_outputs,
    cost_per_dose   = 5.19,
    smc_filter      = 0,
    coverage_filter = "S2-moderate-dropout",
    effect          = "cases"
  )

icers_cases <- run_icers_by_setting(frontier_cases)

p_frontier_cases <- plot_icers_faceted(
  icers_cases,
  label        = "all",
  effect_label = "Cases averted per 100,000 population (thousands)"
)

# Print plots (if running interactively)
print(p_frontier_dalys)
print(p_frontier_cases)

p_frontier_dalys / p_frontier_cases +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme_minimal(8) &
  theme(legend.position = "bottom")

ggsave("./05.plots/publication-plots/appendix/extended-dominance-core-analysis.png", width = 16, height = 8, dpi = 600)

