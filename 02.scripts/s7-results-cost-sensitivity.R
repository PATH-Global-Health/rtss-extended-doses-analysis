# ------------------------------------------------------------------------------
# Purpose: Run one-way sensitivity analyses for cost-effectiveness assumptions.
# Inputs: Processed output dataset and sensitivity assumption grid.
# Outputs: Sensitivity CE planes and tornado plots for appendix materials.
# Dependencies: data_and_libraries.R with plotting and tidyverse utilities.
# Run Stage: Pipeline step 7: sensitivity analyses.
# ------------------------------------------------------------------------------

source("./02.scripts/data_and_libraries.R")
source("./02.scripts/functions/f9-plotting-results.R")

#-READ IN PROCESSED OUTPUTS-----------------------------------------------------
processed_out <- readRDS("./03.outputs/paper-results/processed_outputs.rds")


#-Core function: compute cost & CE for a given assumption set-------------------
# Function: compute_cost_effectiveness
# Purpose: Compute costs and cost-effectiveness metrics for one assumption set.
# Inputs: processed_out and unit-cost assumption parameters.
# Outputs: Data frame with CE outcomes under specified assumptions.
# Assumptions: Assumes processed output schema from s4-processing-model-runs.R.

compute_cost_effectiveness <- function(
    processed_out,
    # vaccine price per dose (excluding consumables)
    rtss_price_per_dose = 4,
    # consumables per dose (keep fixed for sensitivity, update if needed)
    rtss_consumables_per_dose = 1.19,
    # delivery cost per dose
    delivery_cost_seasonal = 1.93,
    delivery_cost_hybrid = 1.22,
    # case management multiplier
    cm_mult = 1,
    # scenario label (for tracking sensitivity scenario)
    scenario_id = "baseline") {
  # ---- 1a. Build cost components ------------------------------------

  # Total vaccine product cost per dose (price + consumables)
  rtss_cost_per_dose <- rtss_price_per_dose

  # SMC and treatment cost data (as in central analysis)
  SMCcost <- 1.07 # per dose including delivery

  RDT <- 0.57 + (0.57 * 0.15)
  AL_adult <- 0.49 * 24
  AL_child <- 0.49 * 12
  outpatient <- 2.19
  inpatient <- 26.06

  # Apply case-management multiplier to the full CM package
  TREATcost_adult <- cm_mult * (RDT + AL_adult + outpatient)
  TREATcost_child <- cm_mult * (RDT + AL_child + outpatient)

  SEVcost_adult <- cm_mult * (RDT + AL_adult + inpatient)
  SEVcost_child <- cm_mult * (RDT + AL_child + inpatient)

  # ---- 1b. Split processed_out into <10 and 0–5 ---------------------

  o5_output <-
    processed_out |>
    filter(age_category == "5+") |>
    select(
      ID, runname, drawID, EIR, warmup, sim_length, population, pfpr,
      seasonality, speciesprop, treatment, SMC, RTSS, dosing_assumption,
      coverage_assumption, seasonal_ages, vaccine_model,
      n_o5 = n,
      o5_cases = cases,
      o5_severe = severe_cases,
      o5_deaths = deaths,
      o5_dalys = daly,
      o5_cases_discounted = cases_discounted,
      o5_severe_discounted = severe_cases_discounted,
      o5_deaths_discounted = deaths_discounted,
      o5_dalys_discounted = daly_discounted
    )

  u5_output <-
    processed_out |>
    filter(age_category == "0-5") |>
    select(
      ID, runname, drawID, EIR, warmup, sim_length, population, pfpr,
      seasonality, speciesprop, treatment, SMC, RTSS, dosing_assumption,
      coverage_assumption, seasonal_ages, vaccine_model,
      n_u5 = n,
      u5_cases = cases,
      u5_severe = severe_cases,
      u5_deaths = deaths,
      u5_dalys = daly,
      u5_cases_discounted = cases_discounted,
      u5_severe_discounted = severe_cases_discounted,
      u5_deaths_discounted = deaths_discounted,
      u5_dalys_discounted = daly_discounted,
      dose1:dose7,
      n_smc_treated,
      dose1_discounted:dose7_discounted,
      n_smc_treated_discounted
    )

  cost_output <-
    left_join(
      u5_output, o5_output,
      by = c(
        "ID", "runname", "drawID", "EIR", "warmup", "sim_length", "population", "pfpr",
        "seasonality", "speciesprop", "treatment", "SMC", "RTSS", "dosing_assumption",
        "coverage_assumption", "seasonal_ages", "vaccine_model"
      )
    )

  # ---- 1c. Vaccine cost per RTSS schedule ---------------------------

  vaccine_costs <-
    tibble(
      RTSS = c("Seasonal", "Hybrid", "none")
    ) |>
    mutate(
      delivery_cost = case_when(
        RTSS == "Seasonal" ~ delivery_cost_seasonal,
        RTSS == "Hybrid" ~ delivery_cost_hybrid,
        RTSS == "none" ~ 0
      ),
      rtss_cost_per_dose = if_else(RTSS == "none", 0, rtss_cost_per_dose)
    )

  cost_output <-
    cost_output |>
    left_join(vaccine_costs, by = "RTSS")

  # ---- 1d. Cost variables -------------------------------------------

  cost_output <-
    cost_output |>
    mutate(
      # uncomplicated treatment costs
      cost_clinical =
        (((o5_cases - o5_severe) * treatment * TREATcost_adult) +
          ((u5_cases - u5_severe) * treatment * TREATcost_child)) * 0.77,

      # severe treatment costs
      cost_severe =
        ((o5_severe * treatment * SEVcost_adult) +
          (u5_severe * treatment * SEVcost_child)) * 0.77,

      # SMC
      cost_smc = n_smc_treated * SMCcost,

      # Vaccine (undiscounted)
      cost_vaccine =
        (dose1 + dose2 + dose3 + dose4 + dose5 + dose6 + dose7) *
          (rtss_cost_per_dose + delivery_cost),
      cost_total = cost_clinical + cost_severe + cost_smc + cost_vaccine,
      cost_total_u5 =
        ((u5_cases - u5_severe) * treatment * TREATcost_child) * 0.77 +
          (u5_severe * treatment * SEVcost_child) * 0.77 +
          cost_smc + cost_vaccine,

      # DISCOUNTED costs
      cost_clinical_discounted =
        (((o5_cases_discounted - o5_severe_discounted) * treatment * TREATcost_adult) +
          ((u5_cases_discounted - u5_severe_discounted) * treatment * TREATcost_child)) * 0.77,
      cost_severe_discounted =
        ((o5_severe_discounted * treatment * SEVcost_adult) +
          (u5_severe_discounted * treatment * SEVcost_child)) * 0.77,
      cost_smc_discounted = n_smc_treated_discounted * SMCcost,
      cost_vaccine_discounted =
        (dose1_discounted + dose2_discounted + dose3_discounted + dose4_discounted +
          dose5_discounted + dose6_discounted + dose7_discounted) *
          (rtss_cost_per_dose + delivery_cost),
      cost_total_discounted =
        cost_clinical_discounted + cost_severe_discounted +
          cost_smc_discounted + cost_vaccine_discounted,
      cost_total_u5_discounted =
        ((u5_cases_discounted - u5_severe_discounted) * treatment * TREATcost_child) * 0.77 +
          (u5_severe_discounted * treatment * SEVcost_child) * 0.77 +
          cost_smc_discounted + cost_vaccine_discounted
    )

  # ---- 1e. Baseline (no vaccine) ------------------------------------

  cost_output_baseline <-
    cost_output |>
    filter(vaccine_model == "no-vaccine") |>
    select(
      drawID, EIR, warmup, sim_length, population, pfpr,
      seasonality, speciesprop, treatment, SMC,
      o5_dalys_baseline = o5_dalys,
      o5_cases_baseline = o5_cases,
      o5_severe_baseline = o5_severe,
      o5_deaths_baseline = o5_deaths,
      o5_dalys_discounted_baseline = o5_dalys_discounted,
      o5_cases_discounted_baseline = o5_cases_discounted,
      o5_severe_discounted_baseline = o5_severe_discounted,
      o5_deaths_discounted_baseline = o5_deaths_discounted,
      u5_dalys_baseline = u5_dalys,
      u5_cases_baseline = u5_cases,
      u5_severe_baseline = u5_severe,
      u5_deaths_baseline = u5_deaths,
      u5_dalys_discounted_baseline = u5_dalys_discounted,
      u5_cases_discounted_baseline = u5_cases_discounted,
      u5_severe_discounted_baseline = u5_severe_discounted,
      u5_deaths_discounted_baseline = u5_deaths_discounted,
      cost_clinical_baseline = cost_clinical,
      cost_severe_baseline = cost_severe,
      cost_smc_baseline = cost_smc,
      cost_vaccine_baseline = cost_vaccine,
      cost_total_baseline = cost_total,
      cost_total_u5_baseline = cost_total_u5,
      cost_clinical_discounted_baseline = cost_clinical_discounted,
      cost_severe_discounted_baseline = cost_severe_discounted,
      cost_smc_discounted_baseline = cost_smc_discounted,
      cost_vaccine_discounted_baseline = cost_vaccine_discounted,
      cost_total_discounted_baseline = cost_total_discounted,
      cost_total_u5_discounted_baseline = cost_total_u5_discounted
    )

  # ---- 1f. Vaccine scenarios with CE --------------------------------

  cost_output_vaccine <-
    cost_output |>
    filter(vaccine_model != "no-vaccine") |>
    left_join(
      cost_output_baseline,
      by = c(
        "drawID", "EIR", "warmup", "sim_length", "population", "pfpr",
        "seasonality", "speciesprop", "treatment", "SMC"
      )
    ) |>
    mutate(
      # undiscounted
      cases_averted = (o5_cases_baseline + u5_cases_baseline) -
        (o5_cases + u5_cases),
      dalys_averted = (o5_dalys_baseline + u5_dalys_baseline) -
        (o5_dalys + u5_dalys),
      incremental_cost = (cost_total - cost_total_baseline),
      cost_saved = (cost_clinical_baseline + cost_severe_baseline) -
        (cost_clinical + cost_severe),
      cost_vax_diff = cost_vaccine - cost_vaccine_baseline,
      cost_effectiveness_daly =
        incremental_cost / dalys_averted,
      cost_effectiveness_case =
        incremental_cost / cases_averted,
      cost_effectiveness_daly_u5 =
        incremental_cost / (u5_dalys_baseline - u5_dalys),
      cost_effectiveness_case_u5 =
        incremental_cost / (u5_cases_baseline - u5_cases),

      # discounted
      cases_averted_discounted =
        (o5_cases_discounted_baseline + u5_cases_discounted_baseline) -
          (o5_cases_discounted + u5_cases_discounted),
      dalys_averted_discounted =
        (o5_dalys_discounted_baseline + u5_dalys_discounted_baseline) -
          (o5_dalys_discounted + u5_dalys_discounted),
      incremental_cost_discounted =
        (cost_total_discounted - cost_total_discounted_baseline),
      cost_saved_discounted =
        (cost_clinical_discounted_baseline + cost_severe_discounted_baseline) -
          (cost_clinical_discounted + cost_severe_discounted),
      cost_vax_diff_discounted =
        cost_vaccine_discounted - cost_vaccine_discounted_baseline,
      cost_effectiveness_daly_discounted =
        incremental_cost_discounted / dalys_averted_discounted,
      cost_effectiveness_case_discounted =
        incremental_cost_discounted / cases_averted_discounted,
      cost_effectiveness_daly_u5_discounted =
        incremental_cost_discounted / (u5_dalys_discounted_baseline - u5_dalys_discounted),
      cost_effectiveness_case_u5_discounted =
        incremental_cost_discounted / (u5_cases_discounted_baseline - u5_cases_discounted)
    )

  cost_output_vaccine
}


#-Define sensitivity scenarios--------------------------------------------------
sens_scenarios <- tibble::tribble(
  ~scenario_id, ~parameter, ~level, ~rtss_price_per_dose, ~delivery_cost_seasonal, ~delivery_cost_hybrid, ~cm_mult,
  "baseline", "Baseline", "central", 4 + 1.19, 1.93, 1.22, 1.0,
  "vax_low", "Vaccine price", "low", 2 + 0.78, 1.93, 1.22, 1.0,
  "vax_high", "Vaccine price", "high", 5 + 1.39, 1.93, 1.22, 1.0,
  "deliv_low", "Delivery cost", "low", 4 + 1.19, 1.37, 0.82, 1.0,
  "deliv_high", "Delivery cost", "high", 4 + 1.19, 2.49, 1.61, 1.0,
  "cm_low", "Case mgmt cost", "low", 4 + 1.19, 1.93, 1.22, 0.5,
  "cm_high", "Case mgmt cost", "high", 4 + 1.19, 1.93, 1.22, 2.0
)

#-Run all scenarios-------------------------------------------------------------
ce_sens <-
  sens_scenarios |>
  mutate(
    data = pmap(
      list(
        scenario_id,
        rtss_price_per_dose,
        delivery_cost_seasonal,
        delivery_cost_hybrid,
        cm_mult
      ),
      ~ compute_cost_effectiveness(
        processed_out              = processed_out,
        scenario_id                = ..1,
        rtss_price_per_dose        = ..2,
        delivery_cost_seasonal     = ..3,
        delivery_cost_hybrid       = ..4,
        cm_mult                    = ..5
      )
    )
  ) |>
  select(scenario_id, parameter, level, data) |>
  unnest(cols = data) |>
  filter(seasonal_ages != "5-months-1-year") |>
  normalize_analysis_labels() |>
  mutate(
    seasonality = str_to_title(seasonality),
    RTSS = paste0(RTSS, " delivery"),
    label = paste0(dosing_assumption, " | ", RTSS),
    sens_id = paste0(parameter, ": ", level),
    pfpr_label = paste0("PfPR ", pfpr * 100, "%")
  )

#-Cost-effectiveness planes (median CE) for all sensitivity scenarios-----------
ce_sens$pfpr_label <- factor(
  ce_sens$pfpr_label,
  levels = c("PfPR 3%", "PfPR 20%", "PfPR 65%")
)

pfpr_vals <- c(0.03, 0.30, 0.65)

sen_pl_1 <-
  ggplot(
    ce_sens |> filter(pfpr %in% pfpr_vals)
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(
    aes(
      x = ((dalys_averted_discounted / (n_u5 + n_o5)) * 100000) / 1000,
      y = ((incremental_cost_discounted / (n_u5 + n_o5)) * 100000) / 1000,
      colour = sens_id
    ),
    alpha = 0.8
  ) +
  facet_grid(pfpr_label ~ label) +
  labs(
    x = "DALYs averted per 100,000 population (thousands)",
    y = "Incremental cost per 100,000\n population (thousands)",
    colour = " "
  ) +
  theme_minimal() +
  scale_colour_viridis_d()

sen_plane_2 <-
  ggplot(
    ce_sens |> filter(pfpr %in% pfpr_vals)
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(
    aes(
      x = ((cases_averted_discounted / (n_u5 + n_o5)) * 100000) / 1000,
      y = ((incremental_cost_discounted / (n_u5 + n_o5)) * 100000) / 1000,
      colour = sens_id
    ),
    alpha = 0.8
  ) +
  facet_grid(pfpr_label ~ label) +
  labs(
    x = "Uncomplicated cases averted per 100,000 population (thousands)",
    y = "Incremental cost per 100,000\n population (thousands)",
    colour = " "
  ) +
  theme_minimal() +
  scale_colour_viridis_d()

sen_plane_2 / sen_pl_1 +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

ggsave("./05.plots/publication-plots/appendix/cost-sens-planes.png", width = 10, height = 8, dpi = 600)

#-Summarise ICERs for chosen PfPR & settings------------------------------------
pfpr_vals <- c(0.03, 0.20, 0.65)

icers_sens <-
  ce_sens |>
  filter(
    pfpr %in% pfpr_vals,
    SMC == 0,
    coverage_assumption == "S2-moderate-dropout"
  ) |>
  group_by(
    scenario_id, parameter, level,
    pfpr, seasonality, RTSS, dosing_assumption,
    SMC, coverage_assumption, label, sens_id, pfpr_label
  ) |>
  summarise(
    cost_per_daly_med = median(cost_effectiveness_daly_discounted, na.rm = TRUE),
    cost_per_case_med = median(cost_effectiveness_case_discounted, na.rm = TRUE),
    incremental_cost_med = median(incremental_cost_discounted, na.rm = TRUE),
    dalys_averted_med = median(dalys_averted_discounted, na.rm = TRUE),
    .groups = "drop"
  )

#-Tornado plots-----------------------------------------------------------------
# Function: make_tornado_single_pfpr
# Purpose: Construct tornado-chart summary for one PfPR sensitivity setting.
# Inputs: sensitivity summary data and selected PfPR/outcome.
# Outputs: ggplot tornado chart for scenario comparison.
# Assumptions: Assumes one-way sensitivity outputs already computed.

make_tornado_single_pfpr <- function(
    pfpr_val,
    seasonality_val,
    outcome = c("case", "daly"),
    main_title = NULL) {
  outcome <- match.arg(outcome)

  if (outcome == "case") {
    baseline_col <- sym("baseline_icer_case")
    value_col <- sym("cost_per_case_med")
    x_lab <- "Cost per uncomplicated case averted (2025 USD)"
  } else {
    baseline_col <- sym("baseline_icer_daly")
    value_col <- sym("cost_per_daly_med")
    x_lab <- "Cost per DALY averted (2025 USD)"
  }

  df_panel <-
    tornado_df |>
    filter(
      SMC == 0,
      pfpr == pfpr_val,
      coverage_assumption == "S2-moderate-dropout",
      seasonality == seasonality_val
    )

  base_panel <-
    baseline_df |>
    filter(
      SMC == 0,
      pfpr == pfpr_val,
      coverage_assumption == "S2-moderate-dropout",
      seasonality == seasonality_val
    )

  pfpr_lab <- unique(df_panel$pfpr_label)

  ggplot(df_panel) +
    # baseline line
    geom_vline(
      data = base_panel,
      aes(xintercept = !!baseline_col),
      inherit.aes = FALSE,
      colour = "black"
    ) +
    # tornado bars
    geom_segment(
      aes(
        x = !!baseline_col,
        xend = !!value_col,
        y = parameter,
        yend = parameter,
        colour = bound
      ),
      linewidth = 6,
      lineend = "butt"
    ) +
    facet_grid(~label, scales = "free_x") +
    scale_colour_manual(
      values = c(
        "Lower bound" = "#D6ED17",
        "Upper bound" = "#606060"
      ),
      name = NULL
    ) +
    labs(
      x = x_lab,
      y = NULL,
      title = main_title,
      subtitle = pfpr_lab[1]
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 8)
    )
}

# Seasonal – one row per PfPR
seas_case_003 <- make_tornado_single_pfpr(0.03, "Seasonal",
  outcome = "case",
  main_title = "Seasonal transmission"
)
seas_case_020 <- make_tornado_single_pfpr(0.20, "Seasonal",
  outcome = "case"
)
seas_case_065 <- make_tornado_single_pfpr(0.65, "Seasonal",
  outcome = "case"
)

# Highly seasonal – one row per PfPR
high_case_003 <- make_tornado_single_pfpr(0.03, "Highly Seasonal",
  outcome = "case",
  main_title = "Highly seasonal transmission"
)
high_case_020 <- make_tornado_single_pfpr(0.20, "Highly Seasonal",
  outcome = "case"
)
high_case_065 <- make_tornado_single_pfpr(0.65, "Highly Seasonal",
  outcome = "case"
)

# Stack rows, then put Seasonal and Highly Seasonal side-by-side
case_tornado_final <-
  (seas_case_003 / seas_case_020 / seas_case_065) |
    (high_case_003 / high_case_020 / high_case_065)

case_tornado_final <-
  case_tornado_final +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

print(case_tornado_final)

ggsave("./05.plots/publication-plots/appendix/cost-sens-tornado-cases.png", width = 15, height = 8, dpi = 600)

# Seasonal – DALYs
seas_daly_003 <- make_tornado_single_pfpr(0.03, "Seasonal",
  outcome = "daly",
  main_title = "Seasonal transmission"
)
seas_daly_020 <- make_tornado_single_pfpr(0.20, "Seasonal",
  outcome = "daly"
)
seas_daly_065 <- make_tornado_single_pfpr(0.65, "Seasonal",
  outcome = "daly"
)

# Highly seasonal – DALYs
high_daly_003 <- make_tornado_single_pfpr(0.03, "Highly Seasonal",
  outcome = "daly",
  main_title = "Highly seasonal transmission"
)
high_daly_020 <- make_tornado_single_pfpr(0.20, "Highly Seasonal",
  outcome = "daly"
)
high_daly_065 <- make_tornado_single_pfpr(0.65, "Highly Seasonal",
  outcome = "daly"
)

daly_tornado_final <-
  (seas_daly_003 / seas_daly_020 / seas_daly_065) |
    (high_daly_003 / high_daly_020 / high_daly_065)

daly_tornado_final <-
  daly_tornado_final +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

print(daly_tornado_final)
ggsave("./05.plots/publication-plots/appendix/cost-sens-tornado-dalys.png", width = 15, height = 8, dpi = 600)
