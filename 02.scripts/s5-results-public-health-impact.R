# ------------------------------------------------------------------------------
# Purpose: Generate public-health impact figures and summary tables for manuscript outputs.
# Inputs: Paper-results RDS outputs produced in step 4.
# Outputs: Main-text and appendix public-health plots/tables.
# Dependencies: data_and_libraries.R and plotting helpers in f9-plotting-results.R.
# Run Stage: Pipeline step 5: manuscript public-health results.
# ------------------------------------------------------------------------------

# Load required libraries
source("./02.scripts/data_and_libraries.R")
source("./02.scripts/functions/f9-plotting-results.R")

library(patchwork)

#-READ IN AND PREP--------------------------------------------------------------------------------------------------------------------
averted_outputs <- readRDS("./03.outputs/paper-results/averted_outputs.rds")
averted_outputs_ab <- readRDS("./03.outputs/paper-results/averted_outputs_age-based.rds")
incremental_outputs <- readRDS("./03.outputs/paper-results/incremental_outputs.rds")
cost_outputs <- readRDS("./03.outputs/paper-results/cost_outputs.rds")

## remove the seasonal age groups that are not used in the analysis-------------
averted_outputs     <- averted_outputs |>  filter(seasonal_ages != "5-months-1-year")
averted_outputs_ab  <- averted_outputs_ab |>  filter(seasonal_ages != "5-months-1-year")
incremental_outputs <- incremental_outputs |>  filter(seasonal_ages != "5-months-1-year")
cost_outputs        <- cost_outputs |>  filter(seasonal_ages != "5-months-1-year")

## rename variables-------------------------------------------------------------
averted_outputs <- normalize_analysis_labels(averted_outputs)
averted_outputs_ab <- normalize_analysis_labels(averted_outputs_ab)
incremental_outputs <- normalize_analysis_labels(incremental_outputs)
cost_outputs <- normalize_analysis_labels(cost_outputs)

#-PART 1 PUBLIC HEALTH IMPACT---------------------------------------------------------------------------------------------------------

# grouping variables to summarise impact by
group_vars <- c(
 "RTSS", "dosing_assumption", "coverage_assumption",
 "seasonality", "pfpr", "SMC", "age_category"
)

# variables we want to summarise over
averted_outcome_vars <- c(
  "cases_averted_discounted",
  "severe_cases_averted_discounted",
  "deaths_averted_discounted"
)

## All results with S2 moderate coverage

##-AVERTED OUTCOMES NO SMC UNDER 5S---------------------------------------------------------------------------

####-Summarised over our uncertainty in intervention models-----------
ph_impact_u5 <-
  averted_outputs |>
  filter(
    age_category == "0-5", # 0-5 year olds main text
    coverage_assumption == "S2-moderate-dropout"
  )  |>
  group_by(across(all_of(group_vars))) |>
  summarise(
    across(
      all_of(averted_outcome_vars),
      list(
        med   = ~ median(.x, na.rm = TRUE),
        l95   = ~ quantile(.x, 0.025, na.rm = TRUE),
        u95   = ~ quantile(.x, 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  mutate(label = paste0(dosing_assumption, " | ", RTSS ),
         seasonality = str_to_title(seasonality),
         RTSS = paste0(RTSS, " delivery")
         )

colors_labels <-
  c("5-doses" = "#94623d",
    "7-doses" = "#5489d0"
  )

###-Plot cases averted with error bars for no SMC @ baseline----------
cases <-
  ggplot(ph_impact_u5 |> filter(SMC == 0)) +
  geom_col(
    aes(
      x = as.factor(pfpr*100),
      y = cases_averted_discounted_med,
      fill = dosing_assumption
    ),
    position = "dodge"
  ) +
  geom_errorbar(
    aes(
      x = as.factor(pfpr*100),
      ymin = cases_averted_discounted_l95,
      ymax = cases_averted_discounted_u95,
      group = label
    ),
    position = "dodge"
  ) +
  labs(
    x = "PfPR2-10 (%)",
    y = "Cases averted in children under 5",
    fill = "Dosing assumption"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(RTSS ~ seasonality, nrow=1) +
  scale_fill_manual(
   values =  colors_labels
  ) +
  theme(legend.position = "top",
        axis.text = element_text(size=8),
        axis.title = element_text(size = 8))

###-Plot deaths averted with error bars for no SMC @ baseline--------
deaths <-
  ggplot(ph_impact_u5 |> filter(SMC == 0)) +
  geom_col(
    aes(
      x = as.factor(pfpr*100),
      y = deaths_averted_discounted_med,
      fill = dosing_assumption
    ),
    position = "dodge"
  ) +
  geom_errorbar(
    aes(
      x = as.factor(pfpr*100),
      ymin = deaths_averted_discounted_l95,
      ymax = deaths_averted_discounted_u95,
      group = label
    ),
    position = "dodge"
  ) +
  labs(
    x = "PfPR2-10 (%)",
    y = "Deaths averted in children under 5",
    fill = "Dosing assumption"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(RTSS ~ seasonality, nrow=1) +
  scale_fill_manual(
    values =  colors_labels
  ) +
  theme(legend.position = "none",
        axis.text = element_text(size=8),
        axis.title = element_text(size = 8))

### Save
public_health_impact_plot <- cases / deaths + plot_annotation(tag_levels = "A")

ggsave("./05.plots/publication-plots/main-text/public-health-impact-plot1.png", plot = public_health_impact_plot, width = 8, height = 6, dpi = 600, bg = "white")
ggsave("./05.plots/publication-plots/main-text/public-health-impact-plot1.pdf", plot = public_health_impact_plot, width = 8, height = 6, dpi = 600, bg = "white")

###-Summary table for appendix and pull to main text----------------
### proportion of cases and deaths averted
### absolute cases and deaths averted
### data point is 30% pfpr no SMC scenario
### Under 5 years, moderate dropout
ph_impact_u5_table <-
  averted_outputs |>
  filter(
    age_category == "0-5", # 0-5 year olds main text
    coverage_assumption == "S2-moderate-dropout",
    SMC == 0
  )  |>
  mutate(perc_cases_averted_discounted = (cases_averted_discounted / cases_discounted_baseline) * 100,
         perc_deaths_averted_discounted = (deaths_averted_discounted / deaths_discounted_baseline) * 100
  ) |>
  group_by(across(all_of(group_vars))) |>
  summarise(
    across(
      c(
        "cases_averted_discounted",
        "deaths_averted_discounted",
        "cases_averted_per_100000_FVC3_discounted",
        "deaths_averted_per_100000_FVC3_discounted",
        "perc_cases_averted_discounted",
        "perc_deaths_averted_discounted"
      ),
      list(
        med   = ~ median(.x, na.rm = TRUE),
        l95   = ~ quantile(.x, 0.025, na.rm = TRUE),
        u95   = ~ quantile(.x, 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  mutate(label = paste0(dosing_assumption, " | ", RTSS ),
         seasonality = str_to_title(seasonality),
         RTSS = paste0(RTSS, " delivery")
  ) |>
  mutate(
    # ROUND ALL NON-PERCENTAGE SUMMARY COLS TO 0 DECIMALS
    across(
      matches("^(cases_|deaths_).*_(med|l95|u95)$"),
      \(x) round(x, 0)
    ),

    # ROUND ALL PERCENTAGE SUMMARY COLS TO 1 DECIMAL
    across(
      matches("^perc_.*_(med|l95|u95)$"),
      \(x) round(x, 1)
    )
  ) |>
  # now make the summary within a single cell median values (lo - high)
  # cases and deaths averted and % reductions different cols
  # build "median (l95–u95)" strings for each metric
  mutate(
    `Unncomplicated cases averted in children younger than 5 years` = sprintf(
      "%d (%d–%d)",
      cases_averted_discounted_med,
      cases_averted_discounted_l95,
      cases_averted_discounted_u95
    ),
    `Deaths averted in children younger than 5 years` = sprintf(
      "%d (%d–%d)",
      deaths_averted_discounted_med,
      deaths_averted_discounted_l95,
      deaths_averted_discounted_u95
    ),
    `Unncomplicated cases averted per 100 000 fully vaccinated children` = sprintf(
      "%d (%d–%d)",
      cases_averted_per_100000_FVC3_discounted_med,
      cases_averted_per_100000_FVC3_discounted_l95,
      cases_averted_per_100000_FVC3_discounted_u95
    ),
    `Deaths averted per 100 000 fully vaccinated children` = sprintf(
      "%d (%d–%d)",
      deaths_averted_per_100000_FVC3_discounted_med,
      deaths_averted_per_100000_FVC3_discounted_l95,
      deaths_averted_per_100000_FVC3_discounted_u95
    ),

    `Proportion of unncomplicated cases averted in children younger than 5 years` = sprintf(
      "%.1f (%.1f–%.1f)",
      perc_cases_averted_discounted_med,
      perc_cases_averted_discounted_l95,
      perc_cases_averted_discounted_u95
    ),
    `Proportion of deaths averted in children younger than 5 years` = sprintf(
      "%.1f (%.1f–%.1f)",
      perc_deaths_averted_discounted_med,
      perc_deaths_averted_discounted_l95,
      perc_deaths_averted_discounted_u95
    )

  ) |>
  # select key variables
  select(
    all_of(group_vars),
    label, seasonality, RTSS,
    `Proportion of unncomplicated cases averted in children younger than 5 years`,
    `Proportion of deaths averted in children younger than 5 years`,
    `Unncomplicated cases averted in children younger than 5 years`,
    `Deaths averted in children younger than 5 years`,
    `Unncomplicated cases averted per 100 000 fully vaccinated children`,
    `Deaths averted per 100 000 fully vaccinated children`

  ) |>
  # pivot wider PfPR to columns
  mutate(
    pfpr = paste0(pfpr*100, "%")
  ) |>
  # pivot LONGER over outcome measures
  tidyr::pivot_longer(
    cols = c(
      `Proportion of unncomplicated cases averted in children younger than 5 years`,
      `Proportion of deaths averted in children younger than 5 years`,
      `Unncomplicated cases averted in children younger than 5 years`,
      `Deaths averted in children younger than 5 years`,
      `Unncomplicated cases averted per 100 000 fully vaccinated children`,
      `Deaths averted per 100 000 fully vaccinated children`
    ),
    names_to  = "outcome",
    values_to = "value"
  ) |>
  # now pivot WIDER so PfPR becomes columns
  tidyr::pivot_wider(
    names_from  = pfpr,
    values_from = value
  )

# save this to excel .xlsx format
writexl::write_xlsx(list(`u5-no-smc-moderate-drop-out` = ph_impact_u5_table),
                     path = "./05.plots/publication-plots/tables/supplementary-tables.xlsx")


##-INCREMENTAL IMPACT NO SMC UNDER 5s------------------------------------------------------------------------

####-Summarised over our uncertainty in intervention models, long format-----------
incremental_outputs <-
  averted_outputs |>
  mutate(
    dosing_num = case_when(
      dosing_assumption == "3-doses" ~ 3,
      dosing_assumption == "4-doses" ~ 4,
      dosing_assumption == "5-doses" ~ 5,
      dosing_assumption == "6-doses" ~ 6,
      dosing_assumption == "7-doses" ~ 7
    )
  ) |>
  inner_join(
    averted_outputs |>
      mutate(
        dosing_num = case_when(
          dosing_assumption == "3-doses" ~ 3,
          dosing_assumption == "4-doses" ~ 4,
          dosing_assumption == "5-doses" ~ 5,
          dosing_assumption == "6-doses" ~ 6,
          dosing_assumption == "7-doses" ~ 7
        )
      ),
    by = c(
      "RTSS", "seasonality", "pfpr", "SMC", "vaccine_model", "drawID", "seasonal_ages", "coverage_assumption",
      "age_category"
    ),
    suffix = c("_from", "_to")
  ) |>
  filter(dosing_num_to > dosing_num_from) |> # only compare to higher doses
  mutate(
    # CASES
    incremental_cases = cases_averted_discounted_to - cases_averted_discounted_from,
    incremental_percent_cases = (incremental_cases / cases_averted_discounted_from) * 100,

    # DEATHS
    incremental_deaths = deaths_averted_discounted_to - deaths_averted_discounted_from,
    incremental_percent_deaths = (incremental_deaths / deaths_averted_discounted_from) * 100,
  ) |>
  distinct() |>
  select(
    RTSS,
    seasonality,
    pfpr,
    vaccine_model,
    seasonal_ages,
    coverage_assumption,
    age_category,
    SMC,
    dosing_assumption_from,
    dosing_assumption_to,
    incremental_percent_cases,
    incremental_percent_deaths
  ) |>
  filter(
    dosing_assumption_to == "7-doses",
    dosing_assumption_from == "5-doses"
  ) |>
  group_by(
    RTSS,
    seasonality,
    pfpr,
    seasonal_ages,
    coverage_assumption,
    age_category,
    SMC,
    dosing_assumption_from,
    dosing_assumption_to,
  ) |>
  summarise(
    across(
      all_of( c("incremental_percent_cases",
              "incremental_percent_deaths")),
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
    RTSS = paste0(RTSS, " delivery")
  )

### long format for plotting
incremental_long <-
  incremental_outputs |>
  filter(
    pfpr > 0.01
  ) |>
  mutate(
    pfpr_percent = pfpr * 100  # nicer x-scale
  ) |>
  # pivot med / l95 / u95 for cases & deaths into long form
  pivot_longer(
    cols = c(
      incremental_percent_cases_med,
      incremental_percent_cases_l95,
      incremental_percent_cases_u95,
      incremental_percent_deaths_med,
      incremental_percent_deaths_l95,
      incremental_percent_deaths_u95
    ),
    names_to   = c("outcome", "stat"),
    names_pattern = "incremental_percent_(cases|deaths)_(med|l95|u95)",
    values_to  = "value"
  ) |>
  # wide again on stat so we have med, l95, u95 columns
  pivot_wider(
    names_from  = stat,
    values_from = value
  ) |>
  mutate(
    outcome = dplyr::recode(
      outcome,
      cases  = "Cases averted",
      deaths = "Deaths averted"
    )
  )

###-Plot incremental % cases and deaths averted for no SMC @ baseline U5s----------
incremental_impact_plot <- ggplot(
  incremental_long |>
    filter(
      age_category == "0-5",
      SMC == 0,
      coverage_assumption == "S2-moderate-dropout",
      ),
  aes(
    x = as.factor(pfpr_percent),
    y = med,
    color = RTSS,
    fill  = RTSS,
    linetype = RTSS
  )
) +
  geom_point() +
  # geom_smooth(
  #   se = FALSE,
  #   show.legend = FALSE
  # ) +
  facet_grid(
    . ~ seasonality  + outcome
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 60)
  ) +
  labs(
    x = "PfPR 2–10 (%)",
    y = "Incremental impact (%)",
    color   = " ",
    fill    = " "
  ) +
  theme_minimal() +
  theme(
    legend.position   = "top",
    legend.box        = "horizontal",
    legend.title.align = 0.5
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.7
  ) +
  scale_color_manual(values = c(
    "Hybrid delivery" = "#C75687",
    "Seasonal delivery" = "#BCBD8B"
  )) +
  scale_fill_manual(values = c(
    "Hybrid delivery" = "#C75687",
    "Seasonal delivery" = "#BCBD8B"
  ))

ggsave("./05.plots/publication-plots/main-text/incrremental-impact-plot1.png", plot = incremental_impact_plot, width = 8, height = 3.5, dpi = 600)
ggsave("./05.plots/publication-plots/main-text/incrremental-impact-plot1.pdf", plot = incremental_impact_plot, width = 8, height = 3.5, dpi = 600)

###-Summary table for appendix and pull to main text----------------
### incremental % cases and deaths averted
### proportion of cases and deaths averted
### Under 5 years, moderate dropout
ph_incremental_u5_table <-
  incremental_outputs |>
  filter(
    age_category == "0-5", # 0-5 year olds main text
    coverage_assumption == "S2-moderate-dropout",
    SMC == 0
  ) |>
  mutate(
    # ROUND ALL PERCENTAGE SUMMARY COLS TO 1 DECIMAL
    across(
      matches("^(incremental_percent_).*_(med|l95|u95)$"),
      \(x) round(x, 1)
    )
  ) |>
  # now make the summary within a single cell median values (lo - high)
  # cases and deaths averted and % reductions different cols
  # build "median (l95–u95)" strings for each metric
  mutate(
    `Incremental proportion of unncomplicated cases averted in children younger than 5 years (5 to 7 doses)` = sprintf(
      "%.1f (%.1f–%.1f)",
      incremental_percent_cases_med,
      incremental_percent_cases_l95,
      incremental_percent_cases_u95
    ),
    `Incremental proportion of deaths averted in children younger than 5 years (5 to 7 doses)` = sprintf(
      "%.1f (%.1f–%.1f)",
      incremental_percent_deaths_med,
      incremental_percent_deaths_l95,
      incremental_percent_deaths_u95
    )
  ) |>
  # select key variables
  select(
    RTSS,
    coverage_assumption,
    seasonality,
    pfpr,
    SMC,
    age_category,
    `Incremental proportion of unncomplicated cases averted in children younger than 5 years (5 to 7 doses)`,
    `Incremental proportion of deaths averted in children younger than 5 years (5 to 7 doses)`
  ) |>
  # pivot wider PfPR to columns
  mutate(
    pfpr = paste0(pfpr*100, "%")
  ) |>
  # pivot LONGER over outcome measures
  tidyr::pivot_longer(
    cols = c(
      `Incremental proportion of unncomplicated cases averted in children younger than 5 years (5 to 7 doses)`,
      `Incremental proportion of deaths averted in children younger than 5 years (5 to 7 doses)`
    ),
    names_to  = "outcome",
    values_to = "value"
  ) |>
  # now pivot WIDER so PfPR becomes columns
  tidyr::pivot_wider(
    names_from  = pfpr,
    values_from = value
  )

# save this to excel .xlsx format
writexl::write_xlsx(list(`Table S2-1` = ph_impact_u5_table,
                         `Table S2-2` = ph_incremental_u5_table),
                    path = "./05.plots/publication-plots/tables/supplementary-tables.xlsx")

##-IMPACT OF SMC AND COVERAGE--------------------------------------------------------------------------

###-Impact on incremental relationship-------------
incremental_long <-
  incremental_long |>
  mutate(SMC = case_when(
    SMC == 0 ~ "No SMC",
    SMC == 0.75 ~ "SMC (75% coverage)",
    TRUE ~ as.character(SMC)
  ))

ggplot(
  incremental_long |>
    filter(
      age_category == "0-5"
    ),
  aes(
    x = as.factor(pfpr_percent),
    y = med,
    color = coverage_assumption,
    fill  = coverage_assumption
  )
) +
  geom_point() +
  facet_grid(
    RTSS ~ seasonality  + outcome + SMC
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 40)
  ) +
  labs(
    x = "PfPR 2–10 (%)",
    y = "Incremental impact (%)",
    color   = " ",
    fill    = " "
  ) +
  theme_minimal() +
  theme(
    legend.position   = "top",
    legend.box        = "horizontal",
    legend.title.align = 0.5
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.7
  )

ggsave("./05.plots/publication-plots/appendix/incrremental-impact-plot-smc-coverage.png", width = 15, height = 6, dpi = 600)

###-Incremental summary table for appendix-------------
ph_incremental_u5_table_smc_and_cov <-
  incremental_outputs |>
  filter(
    age_category == "0-5", # 0-5 year olds main text
  ) |>
  mutate(
    # ROUND ALL PERCENTAGE SUMMARY COLS TO 1 DECIMAL
    across(
      matches("^(incremental_percent_).*_(med|l95|u95)$"),
      \(x) round(x, 1)
    )
  ) |>
  # now make the summary within a single cell median values (lo - high)
  # cases and deaths averted and % reductions different cols
  # build "median (l95–u95)" strings for each metric
  mutate(
    `Incremental proportion of unncomplicated cases averted in children younger than 5 years (5 to 7 doses)` = sprintf(
      "%.1f (%.1f–%.1f)",
      incremental_percent_cases_med,
      incremental_percent_cases_l95,
      incremental_percent_cases_u95
    ),
    `Incremental proportion of deaths averted in children younger than 5 years (5 to 7 doses)` = sprintf(
      "%.1f (%.1f–%.1f)",
      incremental_percent_deaths_med,
      incremental_percent_deaths_l95,
      incremental_percent_deaths_u95
    )
  ) |>
  # select key variables
  select(
    RTSS,
    coverage_assumption,
    seasonality,
    pfpr,
    SMC,
    age_category,
    `Incremental proportion of unncomplicated cases averted in children younger than 5 years (5 to 7 doses)`,
    `Incremental proportion of deaths averted in children younger than 5 years (5 to 7 doses)`
  ) |>
  # pivot wider PfPR to columns
  mutate(
    pfpr = paste0(pfpr*100, "%")
  ) |>
  # pivot LONGER over outcome measures
  tidyr::pivot_longer(
    cols = c(
      `Incremental proportion of unncomplicated cases averted in children younger than 5 years (5 to 7 doses)`,
      `Incremental proportion of deaths averted in children younger than 5 years (5 to 7 doses)`
    ),
    names_to  = "outcome",
    values_to = "value"
  ) |>
  # now pivot WIDER so PfPR becomes columns
  tidyr::pivot_wider(
    names_from  = pfpr,
    values_from = value
  )

# save this to excel .xlsx format
writexl::write_xlsx(list(`Table S2-1` = ph_impact_u5_table,
                         `Table S2-2` = ph_incremental_u5_table,
                         `Table S2-3` = ph_incremental_u5_table_smc_and_cov),
                    path = "./05.plots/publication-plots/tables/supplementary-tables.xlsx")


##-INFLUENCE OF EXTENDING OUT AGE RANGE TO UNDER 10S-------------------------------------------------------------------------

###-Public health impact table------------------------------------------
ph_impact_u10 <-
  averted_outputs |>
  filter(
    age_category == "<10", # 0-5 year olds main text
    coverage_assumption == "S2-moderate-dropout"
  )  |>
  group_by(across(all_of(group_vars))) |>
  summarise(
    across(
      all_of(averted_outcome_vars),
      list(
        med   = ~ median(.x, na.rm = TRUE),
        l95   = ~ quantile(.x, 0.025, na.rm = TRUE),
        u95   = ~ quantile(.x, 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  mutate(label = paste0(dosing_assumption, " | ", RTSS ),
         seasonality = str_to_title(seasonality),
         RTSS = paste0(RTSS, " delivery")
  )

ph_impact_u10_table <-
  averted_outputs |>
  filter(
    age_category == "<10", # 0-5 year olds main text
    coverage_assumption == "S2-moderate-dropout",
    SMC == 0
  )  |>
  mutate(perc_cases_averted_discounted = (cases_averted_discounted / cases_discounted_baseline) * 100,
         perc_deaths_averted_discounted = (deaths_averted_discounted / deaths_discounted_baseline) * 100
  ) |>
  group_by(across(all_of(group_vars))) |>
  summarise(
    across(
      c(
        "cases_averted_discounted",
        "deaths_averted_discounted",
        "cases_averted_per_100000_FVC3_discounted",
        "deaths_averted_per_100000_FVC3_discounted",
        "perc_cases_averted_discounted",
        "perc_deaths_averted_discounted"
      ),
      list(
        med   = ~ median(.x, na.rm = TRUE),
        l95   = ~ quantile(.x, 0.025, na.rm = TRUE),
        u95   = ~ quantile(.x, 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  mutate(label = paste0(dosing_assumption, " | ", RTSS ),
         seasonality = str_to_title(seasonality),
         RTSS = paste0(RTSS, " delivery")
  ) |>
  mutate(
    # ROUND ALL NON-PERCENTAGE SUMMARY COLS TO 0 DECIMALS
    across(
      matches("^(cases_|deaths_).*_(med|l95|u95)$"),
      \(x) round(x, 0)
    ),

    # ROUND ALL PERCENTAGE SUMMARY COLS TO 1 DECIMAL
    across(
      matches("^perc_.*_(med|l95|u95)$"),
      \(x) round(x, 1)
    )
  ) |>
  # now make the summary within a single cell median values (lo - high)
  # cases and deaths averted and % reductions different cols
  # build "median (l95–u95)" strings for each metric
  mutate(
    `Unncomplicated cases averted in children younger than 10 years` = sprintf(
      "%d (%d–%d)",
      cases_averted_discounted_med,
      cases_averted_discounted_l95,
      cases_averted_discounted_u95
    ),
    `Deaths averted in children younger than 10 years` = sprintf(
      "%d (%d–%d)",
      deaths_averted_discounted_med,
      deaths_averted_discounted_l95,
      deaths_averted_discounted_u95
    ),
    `Unncomplicated cases averted per 100 000 fully vaccinated children` = sprintf(
      "%d (%d–%d)",
      cases_averted_per_100000_FVC3_discounted_med,
      cases_averted_per_100000_FVC3_discounted_l95,
      cases_averted_per_100000_FVC3_discounted_u95
    ),
    `Deaths averted per 100 000 fully vaccinated children` = sprintf(
      "%d (%d–%d)",
      deaths_averted_per_100000_FVC3_discounted_med,
      deaths_averted_per_100000_FVC3_discounted_l95,
      deaths_averted_per_100000_FVC3_discounted_u95
    ),

    `Proportion of unncomplicated cases averted in children younger than 10 years` = sprintf(
      "%.1f (%.1f–%.1f)",
      perc_cases_averted_discounted_med,
      perc_cases_averted_discounted_l95,
      perc_cases_averted_discounted_u95
    ),
    `Proportion of deaths averted in children younger than 10 years` = sprintf(
      "%.1f (%.1f–%.1f)",
      perc_deaths_averted_discounted_med,
      perc_deaths_averted_discounted_l95,
      perc_deaths_averted_discounted_u95
    )

  ) |>
  # select key variables
  select(
    all_of(group_vars),
    label, seasonality, RTSS,
    `Proportion of unncomplicated cases averted in children younger than 10 years`,
    `Proportion of deaths averted in children younger than 10 years`,
    `Unncomplicated cases averted in children younger than 10 years`,
    `Deaths averted in children younger than 10 years`,
    `Unncomplicated cases averted per 100 000 fully vaccinated children`,
    `Deaths averted per 100 000 fully vaccinated children`

  ) |>
  # pivot wider PfPR to columns
  mutate(
    pfpr = paste0(pfpr*100, "%")
  ) |>
  # pivot LONGER over outcome measures
  tidyr::pivot_longer(
    cols = c(
      `Proportion of unncomplicated cases averted in children younger than 10 years`,
      `Proportion of deaths averted in children younger than 10 years`,
      `Unncomplicated cases averted in children younger than 10 years`,
      `Deaths averted in children younger than 10 years`,
      `Unncomplicated cases averted per 100 000 fully vaccinated children`,
      `Deaths averted per 100 000 fully vaccinated children`
    ),
    names_to  = "outcome",
    values_to = "value"
  ) |>
  # now pivot WIDER so PfPR becomes columns
  tidyr::pivot_wider(
    names_from  = pfpr,
    values_from = value
  )

# save this to excel .xlsx format
writexl::write_xlsx(list(`Table S2-1` = ph_impact_u5_table,
                         `Table S2-2` = ph_incremental_u5_table,
                         `Table S2-3` = ph_incremental_u5_table_smc_and_cov,
                         `Table S2-4` = ph_impact_u10_table),
                    path = "./05.plots/publication-plots/tables/supplementary-tables.xlsx")

### Incremental table------------------------
ph_incremental_u10_table <-
  incremental_outputs |>
  filter(
    age_category == "<10", # 0-5 year olds main text
    coverage_assumption == "S2-moderate-dropout",
    SMC == 0
  ) |>
  mutate(
    # ROUND ALL PERCENTAGE SUMMARY COLS TO 1 DECIMAL
    across(
      matches("^(incremental_percent_).*_(med|l95|u95)$"),
      \(x) round(x, 1)
    )
  ) |>
  # now make the summary within a single cell median values (lo - high)
  # cases and deaths averted and % reductions different cols
  # build "median (l95–u95)" strings for each metric
  mutate(
    `Incremental proportion of unncomplicated cases averted in children younger than 10 years (5 to 7 doses)` = sprintf(
      "%.1f (%.1f–%.1f)",
      incremental_percent_cases_med,
      incremental_percent_cases_l95,
      incremental_percent_cases_u95
    ),
    `Incremental proportion of deaths averted in children younger than 10 years (5 to 7 doses)` = sprintf(
      "%.1f (%.1f–%.1f)",
      incremental_percent_deaths_med,
      incremental_percent_deaths_l95,
      incremental_percent_deaths_u95
    )
  ) |>
  # select key variables
  select(
    RTSS,
    coverage_assumption,
    seasonality,
    pfpr,
    SMC,
    age_category,
    `Incremental proportion of unncomplicated cases averted in children younger than 10 years (5 to 7 doses)`,
    `Incremental proportion of deaths averted in children younger than 10 years (5 to 7 doses)`
  ) |>
  # pivot wider PfPR to columns
  mutate(
    pfpr = paste0(pfpr*100, "%")
  ) |>
  # pivot LONGER over outcome measures
  tidyr::pivot_longer(
    cols = c(
      `Incremental proportion of unncomplicated cases averted in children younger than 10 years (5 to 7 doses)`,
      `Incremental proportion of deaths averted in children younger than 10 years (5 to 7 doses)`
    ),
    names_to  = "outcome",
    values_to = "value"
  ) |>
  # now pivot WIDER so PfPR becomes columns
  tidyr::pivot_wider(
    names_from  = pfpr,
    values_from = value
  )

# save this to excel .xlsx format
writexl::write_xlsx(list(`Table S2-1` = ph_impact_u5_table,
                         `Table S2-2` = ph_incremental_u5_table,
                         `Table S2-3` = ph_incremental_u5_table_smc_and_cov,
                         `Table S2-4` = ph_impact_u10_table,
                         `Table S2-5` = ph_incremental_u10_table),
                    path = "./05.plots/publication-plots/tables/supplementary-tables.xlsx")

ggplot(
  incremental_long |>
    filter(
      SMC == "No SMC",
      coverage_assumption == "S2-moderate-dropout",
    ),
  aes(
    x = as.factor(pfpr_percent),
    y = med,
    color = RTSS,
    fill  = RTSS,
    alpha = age_category
  )
) +
  geom_point() +
  # geom_smooth(
  #   se = FALSE,
  #   show.legend = FALSE
  # ) +
  facet_grid(
    . ~ seasonality  + outcome
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 60)
  ) +
  labs(
    x = "PfPR 2–10 (%)",
    y = "Incremental impact (%)",
    color   = " ",
    fill    = " ",
    alpha = " "
  ) +
  theme_minimal() +
  theme(
    legend.position   = "top",
    legend.box        = "horizontal",
    legend.title.align = 0.5
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.7
  ) +
  scale_color_manual(values = c(
    "#C75687",
    "#BCBD8B"
  )) +
  scale_alpha_manual(
    values = c(
      "0-5" = 0.5,
      "<10" = 1
    )
  )

ggsave("./05.plots/publication-plots/appendix/incrremental-impact--age-comp.png", width = 8, height = 3.5, dpi = 600)

##-AGE DISAGGREGATED IMPACTS-----------------------------------------------------------------------------
ph_impact_ab <-
  averted_outputs_ab |>
  ungroup() |>
  filter(
    coverage_assumption == "S2-moderate-dropout"
  )  |>
  group_by(
    RTSS, dosing_assumption, coverage_assumption,
    seasonality, pfpr, SMC, age_lower
    ) |>
  summarise(
    across(
      all_of(averted_outcome_vars),
      list(
        med   = ~ median(.x, na.rm = TRUE),
        l95   = ~ quantile(.x, 0.025, na.rm = TRUE),
        u95   = ~ quantile(.x, 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  mutate(label = paste0(dosing_assumption, " | ", RTSS ),
         seasonality = str_to_title(seasonality),
         RTSS = paste0(RTSS, " delivery"),
         pfpr_label = paste0("PfPR: ", pfpr*100, "%")
  )

ph_impact_ab$pfpr_label <- factor(
  ph_impact_ab$pfpr_label,
  levels = c("PfPR: 3%", "PfPR: 20%", "PfPR: 65%")
)

ggplot(ph_impact_ab |>
         filter(
           SMC == 0,
           pfpr %in% c(0.03, 0.2, 0.65)
           )
       ) +
  geom_col(
    aes(
      x = as.factor(age_lower),
      y = cases_averted_discounted_med,
      fill = dosing_assumption
    ),
    position = "dodge"
  ) +
  geom_errorbar(
    aes(
      x = as.factor(age_lower),
      ymin = cases_averted_discounted_l95,
      ymax = cases_averted_discounted_u95,
      group = label
    ),
    position = "dodge"
  ) +
  labs(
    x = "Age Group (years)",
    y = "Uncomplicated cases averted",
    fill = "Dosing assumption"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma) +
  facet_grid(pfpr_label + RTSS ~ seasonality , scales="free_y") +
  scale_fill_manual(
    values =  colors_labels
  ) +
  theme(legend.position = "top",
        axis.text = element_text(size=8),
        axis.title = element_text(size = 8))




ggplot(ph_impact_ab |>
         filter(
           SMC == 0,
           pfpr %in% c(0.03, 0.2, 0.65)
         )
) +
  geom_col(
    aes(
      x = as.factor(age_lower),
      y = deaths_averted_discounted_med,
      fill = dosing_assumption
    ),
    position = "dodge"
  ) +
  geom_errorbar(
    aes(
      x = as.factor(age_lower),
      ymin = deaths_averted_discounted_l95,
      ymax = deaths_averted_discounted_u95,
      group = label
    ),
    position = "dodge"
  ) +
  labs(
    x = "Age Group (years)",
    y = "Malaria deaths averted",
    fill = "Dosing assumption"
  ) +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma) +
  facet_grid(pfpr_label + RTSS ~ seasonality , scales="free_y") +
  scale_fill_manual(
    values =  colors_labels
  ) +
  theme(legend.position = "top",
        axis.text = element_text(size=8),
        axis.title = element_text(size = 8))
