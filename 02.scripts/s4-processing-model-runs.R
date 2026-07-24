# ------------------------------------------------------------------------------
# Purpose: Process raw model outputs into paper-ready epidemiologic and CE datasets.
# Inputs: Individual simulation outputs and scenario metadata tables.
# Outputs: Processed outputs and paper-results RDS artifacts used by result scripts.
# Dependencies: data_and_libraries.R, f4-processing-epi.R, f5-deaths-dalys.R.
# Run Stage: Pipeline step 4: post-processing and derived outcomes.
# ------------------------------------------------------------------------------

#-OUTLINE-------------------------------------------------------------------------------
# This script is for the main paper analyses
# It reads in the outputs from the model simulations, processes them to calculate outcomes
# by age group, and then calculates events averted relative to the baseline no vaccine scenario.
# It also calculates the incremental impact of additional doses 
# The script is structured as follows:
# 1. Read in and combine all model outputs
# 2. Process outputs to calculate age-specific outcomes (cases, deaths, DALYs)
# 3. Calculate events averted relative to baseline no vaccine scenario
# 4. Calculate incremental impact of additional doses
# 5. Calculate cost-effectiveness metrics (ICERs) for vaccine scenarios compared to baseline
# 6. Sub-analysis using 5-dose schedule as the baseline for comparison
# 7. Save processed outputs for use in paper figures and tables
# Note: This script assumes that the model simulation outputs have already been generated and saved in the specified directory.
#-------------------------------------------------------------------------------

#-SOURCE FUNCTIONS AND JOIN ALL OUTPUT--------------------------------------------------
source("./02.scripts/data_and_libraries.R")
source("./02.scripts/functions/f4-processing-epi.R")
source("./02.scripts/functions/f5-deaths-dalys.R")

# read all scenarios in
paramlist <- readRDS("./01.data/vaccine-scenarios/run_parameters.rds")
x <- c(1:nrow(readRDS("./01.data/vaccine-scenarios/runscenarios.rds"))) # number of runs

# define all combinations of scenarios and draws - only run draw 0
# index <- tibble(x = x)
# Filter for drawID == 0 and get the row numbers
rows_to_run <- which(paramlist$drawID == 0)
runname_to_run <- paramlist$run_name[rows_to_run]

# Create index tibble
index <- tibble(x = rows_to_run, rn = runname_to_run)

# Set up parallel processing
library(furrr)

# keep FALSE for batch runs where occasional scenario failures are acceptable
stop_on_processing_error <- FALSE

# Safety wrapper around processing function
safe_processing <- purrr::safely(
  processing,
  otherwise = NULL
)

# Set up parallel processing (default 50; override with S4_WORKERS)
workers_s4 <- suppressWarnings(as.integer(Sys.getenv("S4_WORKERS", "50")))
if (is.na(workers_s4) || workers_s4 < 1) workers_s4 <- 50
plan(multisession, workers = workers_s4)

# Run the full loop in parallel, safely
results <- future_map2(
  index$x,
  index$rn,
  ~ safe_processing(.x, .y),
  .progress = TRUE
)

# Shut down parallel plan
plan(sequential)

# List and stop on failed runs
failed_idx <- which(map_lgl(results, ~ !is.null(.x$error)))
failed_runs <- index[failed_idx, ]

if (nrow(failed_runs) > 0) {
  message("Processing completed with errors in ", nrow(failed_runs), " runs; continuing.")
  print(failed_runs)
  if (isTRUE(stop_on_processing_error)) {
    stop("Stopping because one or more processing jobs failed.")
  }
} else {
  message("Processing completed with no caught errors.")
}

# read and join function
# Function: read_and_join
# Purpose: Read one processed RDS file for row-binding into a combined table.
# Inputs: path to processed .rds file.
# Outputs: Data frame read from the file path.
# Assumptions: Assumes the supplied path points to an RDS with consistent schema.

read_and_join <- function(path) {
  if (grepl(".rds", path)) {
    out <- read_rds(path)
    return(out)
  }
}

# all files to combine
y <-
  index |>
  mutate(path = paste0("./03.outputs/processed-model-sim-outputs/", x, "-", rn, ".rds")) |>
  mutate(file_exists = file.exists(path)) |>
  filter(file_exists) |>
  select(path) |>
  as.list()

# combined outputs
raw_out <-
  data.frame(rbindlist(purrr::pmap(y, read_and_join), fill = TRUE))

#-ADD EPI INDICATORS--------------------------------------------------------------------
# Age groupings, deaths and DALYs
processed_out <-
  raw_out |>
  mortality_rate() |>
  outcome_uncertainty() |>
  daly_components() |>
  ungroup() |>
  mutate(
    age_category = case_when(
      # standalone under-5 group (use this for u5)
      age_lower == 0 & age_upper == 5 ~ "0-5",

      # under-10 group: ONLY the 1-year bins 0–1 ... 9–10 (exclude standalone 0–5)
      age_lower %in% 0:9 & age_upper == age_lower + 1 ~ "0-10",

      # 5+ group: everything starting at 5 (non-overlapping with standalone 0–5)
      age_lower >= 5 ~ "5+",
      TRUE ~ NA_character_
    )
  ) |>
  filter(age_category %in% c("0-5", "0-10", "5+")) |>
  select(-age, -age_lower, -age_upper) |>
  group_by(
    ID, runname, drawID, EIR, warmup, sim_length, population, pfpr,
    seasonality, speciesprop, treatment, SMC, RTSS, dosing_assumption,
    coverage_assumption, seasonal_ages, vaccine_model, age_category
  ) |>
  summarise(
    across(where(is.numeric), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  distinct() |>
  rename(n_smc_treated_discounted = n_smc_treeated_discounted)

# save data
saveRDS(processed_out, "./03.outputs/paper-results/processed_outputs.rds")

#-EVENTS AVERTED-------------------------------------------------------------------------

## -Define baseline no vaccination runs-------------
baseline_data <-
  processed_out |>
  filter(vaccine_model == "no-vaccine") |>
  rename(
    n_treated_baseline = n_treated,
    dose1_baseline = dose1,
    dose2_baseline = dose2,
    dose3_baseline = dose3,
    dose4_baseline = dose4,
    dose5_baseline = dose5,
    dose6_baseline = dose6,
    dose7_baseline = dose7,
    n_smc_treated_baseline = n_smc_treated,
    n_baseline = n,
    n_treated_discounted_baseline = n_treated_discounted,
    dose1_discounted_baseline = dose1_discounted,
    dose2_discounted_baseline = dose2_discounted,
    dose3_discounted_baseline = dose3_discounted,
    dose4_discounted_baseline = dose4_discounted,
    dose5_discounted_baseline = dose5_discounted,
    dose6_discounted_baseline = dose6_discounted,
    dose7_discounted_baseline = dose7_discounted,
    n_smc_treated_discounted_baseline = n_smc_treated_discounted,
    inc_baseline = inc,
    sev_baseline = sev,
    inc_discounted_baseline = inc_discounted,
    sev_discounted_baseline = sev_discounted,
    cases_baseline = cases,
    severe_cases_baseline = severe_cases,
    cases_discounted_baseline = cases_discounted,
    severe_cases_discounted_baseline = severe_cases_discounted,
    mortality_rate_baseline = mortality_rate,
    deaths_baseline = deaths,
    deaths_discounted_baseline = deaths_discounted,
    cases_lower_baseline = cases_lower,
    cases_upper_baseline = cases_upper,
    deaths_lower_baseline = deaths_lower,
    deaths_upper_baseline = deaths_upper,
    cases_discounted_lower_baseline = cases_discounted_lower,
    cases_discounted_upper_baseline = cases_discounted_upper,
    deaths_discounted_lower_baseline = deaths_discounted_lower,
    deaths_discounted_upper_baseline = deaths_discounted_upper,
    yll_baseline = yll,
    yll_lower_baseline = yll_lower,
    yll_upper_baseline = yll_upper,
    yld_baseline = yld,
    yld_lower_baseline = yld_lower,
    yld_upper_baseline = yld_upper,
    yll_discounted_baseline = yll_discounted,
    yll_discounted_lower_baseline = yll_discounted_lower,
    yll_discounted_upper_baseline = yll_discounted_upper,
    yld_discounted_baseline = yld_discounted,
    yld_discounted_lower_baseline = yld_discounted_lower,
    yld_discounted_upper_baseline = yld_discounted_upper,
    daly_baseline = daly,
    daly_lower_baseline = daly_lower,
    daly_upper_baseline = daly_upper,
    daly_discounted_baseline = daly_discounted,
    daly_discounted_lower_baseline = daly_discounted_lower,
    daly_discounted_upper_baseline = daly_discounted_upper
  ) |>
  distinct() |>
  select(
    -vaccine_model, -coverage_assumption, -RTSS, -dosing_assumption,
    -seasonal_ages, -ID, -runname
  ) |>
  distinct()

## -Calculate events averted relative to baseline no vaccine---------------------------------
averted_outputs <-
  processed_out |>
  filter(vaccine_model != "no-vaccine") |>
  left_join(
    baseline_data,
    by = c(
      "drawID", "EIR", "warmup", "sim_length", "population", "pfpr", "seasonality", "speciesprop", "treatment", "SMC",
      "age_category"
    )
  ) |>
  # Processing the outputs to calculate cases, deaths, and DALYs averted
  mutate(
    # model outputted cases averted
    cases_averted = cases_baseline - cases,
    cases_averted_discounted = cases_discounted_baseline - cases_discounted,
    # model outputted severe cases averted
    severe_cases_averted = severe_cases_baseline - severe_cases,
    severe_cases_averted_discounted = severe_cases_discounted_baseline - severe_cases_discounted,
    # model outputted deaths averted
    deaths_averted = deaths_baseline - deaths,
    deaths_averted_discounted = deaths_discounted_baseline - deaths_discounted,
    # model outputted DALYs averted
    dalys_averted = daly_baseline - daly,
    dalys_averted_discounted = daly_discounted_baseline - daly_discounted,
    # processed cases averted per 100,000 children
    cases_averted_per_100000 = ((cases_baseline / n_baseline) - (cases / n)) * 100000,
    cases_averted_per_100000_discounted = ((cases_discounted_baseline / n_baseline) - (cases_discounted / n)) * 100000,
    # processed severe cases averted per 100,000 children
    severe_cases_averted_per_100000 = ((severe_cases_baseline / n_baseline) - (severe_cases / n)) * 100000,
    severe_cases_averted_per_100000_discounted = ((severe_cases_discounted_baseline / n_baseline) - (severe_cases_discounted / n)) * 100000,
    # processed deaths averted per 100,000 children
    deaths_averted_per_100000 = ((deaths_baseline / n_baseline) - (deaths / n)) * 100000,
    deaths_averted_per_100000_discounted = ((deaths_discounted_baseline / n_baseline) - (deaths_discounted / n)) * 100000,
    # processed DALYs averted per 100,000 children
    dalys_averted_per_100000 = ((daly_baseline / n_baseline) - (daly / n)) * 100000,
    dalys_averted_per_100000_discounted = ((daly_discounted_baseline / n_baseline) - (daly_discounted / n)) * 100000,
    # processed cases averted per 100,000 FVC (3 doses)
    cases_averted_per_100000_FVC3 = 100000 * (cases_averted / dose3),
    cases_averted_per_100000_FVC3_discounted = 100000 * (cases_averted_discounted / dose3_discounted),
    # processed severe cases averted per 100,000 FVC (3 doses)
    severe_cases_averted_per_100000_FVC3 = 100000 * (severe_cases_averted / dose3),
    severe_cases_averted_per_100000_FVC3_discounted = 100000 * (severe_cases_averted_discounted / dose3_discounted),
    # processed deaths averted per 100,000 FVC (3 doses)
    deaths_averted_per_100000_FVC3 = 100000 * (deaths_averted / dose3),
    deaths_averted_per_100000_FVC3_discounted = 100000 * (deaths_averted_discounted / dose3_discounted),
    # processed DALYs averted per 100,000 FVC (3 doses)
    dalys_averted_per_100000_FVC3 = 100000 * (dalys_averted / dose3),
    dalys_averted_per_100000_FVC3_discounted = 100000 * (dalys_averted_discounted / dose3_discounted)
  ) |>
  # add averted per 100,000 FVC depending on dosing assumption
  mutate(
    cases_averted_per_100000_FVCdose_dependent = case_when(
      dosing_assumption == "3-doses" ~ cases_averted_per_100000_FVC3,
      dosing_assumption == "4-doses" ~ 100000 * (cases_averted / dose4),
      dosing_assumption == "5-doses" ~ 100000 * (cases_averted / dose5),
      dosing_assumption == "6-doses" ~ 100000 * (cases_averted / dose6),
      dosing_assumption == "7-doses" ~ 100000 * (cases_averted / dose7)
    ),
    cases_averted_per_100000_FVCdose_dependent_discounted = case_when(
      dosing_assumption == "3-doses" ~ cases_averted_per_100000_FVC3_discounted,
      dosing_assumption == "4-doses" ~ 100000 * (cases_averted_discounted / dose4_discounted),
      dosing_assumption == "5-doses" ~ 100000 * (cases_averted_discounted / dose5_discounted),
      dosing_assumption == "6-doses" ~ 100000 * (cases_averted_discounted / dose6_discounted),
      dosing_assumption == "7-doses" ~ 100000 * (cases_averted_discounted / dose7_discounted)
    ),
    severe_cases_averted_per_100000_FVCdose_dependent = case_when(
      dosing_assumption == "3-doses" ~ severe_cases_averted_per_100000_FVC3,
      dosing_assumption == "4-doses" ~ 100000 * (severe_cases_averted / dose4),
      dosing_assumption == "5-doses" ~ 100000 * (severe_cases_averted / dose5),
      dosing_assumption == "6-doses" ~ 100000 * (severe_cases_averted / dose6),
      dosing_assumption == "7-doses" ~ 100000 * (severe_cases_averted / dose7)
    ),
    severe_cases_averted_per_100000_FVCdose_dependent_discounted = case_when(
      dosing_assumption == "3-doses" ~ severe_cases_averted_per_100000_FVC3_discounted,
      dosing_assumption == "4-doses" ~ 100000 * (severe_cases_averted_discounted / dose4_discounted),
      dosing_assumption == "5-doses" ~ 100000 * (severe_cases_averted_discounted / dose5_discounted),
      dosing_assumption == "6-doses" ~ 100000 * (severe_cases_averted_discounted / dose6_discounted),
      dosing_assumption == "7-doses" ~ 100000 * (severe_cases_averted_discounted / dose7_discounted)
    ),
    deaths_averted_per_100000_FVCdose_dependent = case_when(
      dosing_assumption == "3-doses" ~ deaths_averted_per_100000_FVC3,
      dosing_assumption == "4-doses" ~ 100000 * (deaths_averted / dose4),
      dosing_assumption == "5-doses" ~ 100000 * (deaths_averted / dose5),
      dosing_assumption == "6-doses" ~ 100000 * (deaths_averted / dose6),
      dosing_assumption == "7-doses" ~ 100000 * (deaths_averted / dose7)
    ),
    deaths_averted_per_100000_FVCdose_dependent_discounted = case_when(
      dosing_assumption == "3-doses" ~ deaths_averted_per_100000_FVC3_discounted,
      dosing_assumption == "4-doses" ~ 100000 * (deaths_averted_discounted / dose4_discounted),
      dosing_assumption == "5-doses" ~ 100000 * (deaths_averted_discounted / dose5_discounted),
      dosing_assumption == "6-doses" ~ 100000 * (deaths_averted_discounted / dose6_discounted),
      dosing_assumption == "7-doses" ~ 100000 * (deaths_averted_discounted / dose7_discounted)
    ),
    dalys_averted_per_100000_FVCdose_dependent = case_when(
      dosing_assumption == "3-doses" ~ dalys_averted_per_100000_FVC3,
      dosing_assumption == "4-doses" ~ 100000 * (dalys_averted / dose4),
      dosing_assumption == "5-doses" ~ 100000 * (dalys_averted / dose5),
      dosing_assumption == "6-doses" ~ 100000 * (dalys_averted / dose6),
      dosing_assumption == "7-doses" ~ 100000 * (dalys_averted / dose7)
    ),
    dalys_averted_per_100000_FVCdose_dependent_discounted = case_when(
      dosing_assumption == "3-doses" ~ dalys_averted_per_100000_FVC3_discounted,
      dosing_assumption == "4-doses" ~ 100000 * (dalys_averted_discounted / dose4_discounted),
      dosing_assumption == "5-doses" ~ 100000 * (dalys_averted_discounted / dose5_discounted),
      dosing_assumption == "6-doses" ~ 100000 * (dalys_averted_discounted / dose6_discounted),
      dosing_assumption == "7-doses" ~ 100000 * (dalys_averted_discounted / dose7_discounted)
    )
  ) |>
  # averted per 100,000 doses delivered
  mutate(
    cases_averted_per_100000_doses = 100000 * (cases_averted / (dose1 + dose2 + dose3 + dose4 + dose5 + dose6 + dose7)),
    cases_averted_per_100000_doses_discounted = 100000 * (cases_averted_discounted / (dose1_discounted + dose2_discounted + dose3_discounted + dose4_discounted + dose5_discounted + dose6_discounted + dose7_discounted)),
    severe_cases_averted_per_100000_doses = 100000 * (severe_cases_averted / (dose1 + dose2 + dose3 + dose4 + dose5 + dose6 + dose7)),
    severe_cases_averted_per_100000_doses_discounted = 100000 * (severe_cases_averted_discounted / (dose1_discounted + dose2_discounted + dose3_discounted + dose4_discounted + dose5_discounted + dose6_discounted + dose7_discounted)),
    deaths_averted_per_100000_doses = 100000 * (deaths_averted / (dose1 + dose2 + dose3 + dose4 + dose5 + dose6 + dose7)),
    deaths_averted_per_100000_doses_discounted = 100000 * (deaths_averted_discounted / (dose1_discounted + dose2_discounted + dose3_discounted + dose4_discounted + dose5_discounted + dose6_discounted + dose7_discounted)),
    dalys_averted_per_100000_doses = 100000 * (dalys_averted / (dose1 + dose2 + dose3 + dose4 + dose5 + dose6 + dose7)),
    dalys_averted_per_100000_doses_discounted = 100000 * (dalys_averted_discounted / (dose1_discounted + dose2_discounted + dose3_discounted + dose4_discounted + dose5_discounted + dose6_discounted + dose7_discounted))
  )

## save
saveRDS(averted_outputs, "./03.outputs/paper-results/averted_outputs.rds")

#-INCREMENTAL IMPACT OF ADDITIONAL DOSES------------------------------------------------
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
    incremental_cases_FVC3 = cases_averted_per_100000_FVC3_to - cases_averted_per_100000_FVC3_from,
    incremental_cases_FVCdose = cases_averted_per_100000_FVCdose_dependent_to - cases_averted_per_100000_FVCdose_dependent_from,
    incremental_cases_total_doses = cases_averted_per_100000_doses_to - cases_averted_per_100000_doses_from,
    incremental_cases_population = cases_averted_per_100000_to - cases_averted_per_100000_from,
    incremental_percent_cases_FVC3 = (incremental_cases_FVC3 / cases_averted_per_100000_FVC3_from) * 100,
    incremental_percent_cases_FVCdose = (incremental_cases_FVCdose / cases_averted_per_100000_FVCdose_dependent_from) * 100,
    incremental_percent_cases_total_doses = (incremental_cases_total_doses / cases_averted_per_100000_doses_from) * 100,
    incremental_percent_cases_population = (incremental_cases_population / cases_averted_per_100000_from) * 100,

    # SEVERE CASES
    incremental_severe_cases_FVC3 = severe_cases_averted_per_100000_FVC3_to - severe_cases_averted_per_100000_FVC3_from,
    incremental_severe_cases_FVCdose = severe_cases_averted_per_100000_FVCdose_dependent_to - severe_cases_averted_per_100000_FVCdose_dependent_from,
    incremental_severe_cases_total_doses = severe_cases_averted_per_100000_doses_to - severe_cases_averted_per_100000_doses_from,
    incremental_severe_cases_population = severe_cases_averted_per_100000_to - severe_cases_averted_per_100000_from,
    incremental_percent_severe_cases_FVC3 = (incremental_severe_cases_FVC3 / severe_cases_averted_per_100000_FVC3_from) * 100,
    incremental_percent_severe_cases_FVCdose = (incremental_severe_cases_FVCdose / severe_cases_averted_per_100000_FVCdose_dependent_from) * 100,
    incremental_percent_severe_cases_total_doses = (incremental_severe_cases_total_doses / severe_cases_averted_per_100000_doses_from) * 100,
    incremental_percent_severe_cases_population = (incremental_severe_cases_population / severe_cases_averted_per_100000_from) * 100,

    # DEATHS
    incremental_deaths_FVC3 = deaths_averted_per_100000_FVC3_to - deaths_averted_per_100000_FVC3_from,
    incremental_deaths_FVCdose = deaths_averted_per_100000_FVCdose_dependent_to - deaths_averted_per_100000_FVCdose_dependent_from,
    incremental_deaths_total_doses = deaths_averted_per_100000_doses_to - deaths_averted_per_100000_doses_from,
    incremental_deaths_population = deaths_averted_per_100000_to - deaths_averted_per_100000_from,
    incremental_percent_deaths_FVC3 = (incremental_deaths_FVC3 / deaths_averted_per_100000_FVC3_from) * 100,
    incremental_percent_deaths_FVCdose = (incremental_deaths_FVCdose / deaths_averted_per_100000_FVCdose_dependent_from) * 100,
    incremental_percent_deaths_total_doses = (incremental_deaths_total_doses / deaths_averted_per_100000_doses_from) * 100,
    incremental_percent_deaths_population = (incremental_deaths_population / deaths_averted_per_100000_from) * 100,

    # DALYS
    incremental_dalys_FVC3 = dalys_averted_per_100000_FVC3_to - dalys_averted_per_100000_FVC3_from,
    incremental_dalys_FVCdose = dalys_averted_per_100000_FVCdose_dependent_to - dalys_averted_per_100000_FVCdose_dependent_from,
    incremental_dalys_total_doses = dalys_averted_per_100000_doses_to - dalys_averted_per_100000_doses_from,
    incremental_dalys_population = dalys_averted_per_100000_to - dalys_averted_per_100000_from,
    incremental_percent_dalys_FVC3 = (incremental_dalys_FVC3 / dalys_averted_per_100000_FVC3_from) * 100,
    incremental_percent_dalys_FVCdose = (incremental_dalys_FVCdose / dalys_averted_per_100000_FVCdose_dependent_from) * 100,
    incremental_percent_dalys_total_doses = (incremental_dalys_total_doses / dalys_averted_per_100000_doses_from) * 100,
    incremental_percent_dalys_population = (incremental_dalys_population / dalys_averted_per_100000_from) * 100
  ) |>
  distinct()

saveRDS(incremental_outputs, "./03.outputs/paper-results/incremental_outputs.rds")


#-COST EFFECTIVENESS CALCULATIONS---------------------------------------------------------

##-Define cost data---------------------------------------------
# vaccine cost data
rtss_cost_per_dose <- c(2, 4, 5) # RTS,S per dose assunmption
rtss_consumables_per_dose <- c(0.78, 1.19, 1.39) # https://doi.org/10.1136/bmjgh-2022-011316 "Wastage (vaccine)", "Injection syringe", "Reconstitution syringe", "Safety box (100-capacity)", "Wastage (injection device and safety boxes)", "Vaccine/injection device and safety boxes (buffer stock)", "Freight and insurance", "Handlig fee"

rtss_cost_per_dose <- rtss_cost_per_dose + rtss_consumables_per_dose

delivery_cost_seasonal <- 1.93 # 10.1136/bmjgh-2022-011316 cost of delivery per dose (economic) Table 3
delivery_cost_hybrid <- 1.22 # 10.1136/bmjgh-2022-011316 cost of delivery per dose (economic) Table 3

# SMC and treatment cost data
SMCcost <- 1.07 # per dose including delivery - Gilmartin et al
RDT <- 0.57 + (0.57 * 0.15) # RDT $0.57 unit cost + (unit cost * 15% delivery markup)
AL_adult <- 0.49 * 24 # clinical treatment cost ($0.47 * 24 doses)
AL_child <- 0.49 * 12 # clinical treatment cost ($0.47 * 12 doses)
outpatient <- 2.19 # (Median WHO Choice cost for SSA)
inpatient <- 26.06 # (Median WHO Choice cost for SSA, assuming average duration of stay of 3 days)

# clinical: RDT cost + Drug course cost + facility cost (outpatient)
TREATcost_adult <- RDT + AL_adult + outpatient
TREATcost_child <- RDT + AL_child + outpatient

# severe: RDT cost + Drug course cost + facility cost (inpatient)
SEVcost_adult <- RDT + AL_adult + inpatient
SEVcost_child <- RDT + AL_child + inpatient

## -Isolate age group outputs for child and adult costs-----------
# Isolate over 5 data needed for adult treatment cost processing and CE analysis
o5_output <-
  processed_out |>
  filter(age_category == "5+") |>
  select(
    ID, runname, drawID, EIR, warmup, sim_length, population, pfpr,
    seasonality, speciesprop, treatment, SMC, RTSS, dosing_assumption,
    coverage_assumption, seasonal_ages, vaccine_model,
    # rename the over 5 data into columns
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

# Isolate under 5 data needed for child treatment cost processing and CE analysis
u5_output <-
  processed_out |>
  filter(age_category == "0-5") |>
  select(
    ID, runname, drawID, EIR, warmup, sim_length, population, pfpr,
    seasonality, speciesprop, treatment, SMC, RTSS, dosing_assumption,
    coverage_assumption, seasonal_ages, vaccine_model,
    # rename the under 5 data into columns
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

# Join the two datasets together with isolated column values now for age groups
cost_output <-
  left_join(u5_output, o5_output,
    by = c(
      "ID", "runname", "drawID", "EIR", "warmup", "sim_length", "population", "pfpr",
      "seasonality", "speciesprop", "treatment", "SMC", "RTSS", "dosing_assumption",
      "coverage_assumption", "seasonal_ages", "vaccine_model"
    )
  )

##-Add costs into scenarios---------------------------------------
# add vaccine cost data to the dataset by vaccine schedule RTSS
vaccine_costs <-
  tibble(
    RTSS = c("Seasonal", "Hybrid", "none")
  ) |>
  crossing(
    rtss_cost_per_dose
  ) |>
  mutate(
    delivery_cost = case_when(
      RTSS == "Seasonal" ~ delivery_cost_seasonal,
      RTSS == "Hybrid" ~ delivery_cost_hybrid,
      RTSS == "none" ~ 0,
    )
  ) |>
  # make rtss_cost_per_dose 0
  mutate(
    rtss_cost_per_dose = ifelse(RTSS == "none", 0, rtss_cost_per_dose)
  ) |>
  distinct()

cost_output <-
  cost_output |>
  left_join(vaccine_costs, by = "RTSS")

# Calculate costs of each scenario

cost_output <-
  cost_output |>
  mutate(

    ## uncomplicated treatment costs
    cost_clinical =
    # adult treatment component
      (((o5_cases - o5_severe) * treatment * TREATcost_adult) +
        # child treatment component
        ((u5_cases - u5_severe) * treatment * TREATcost_child)) * 0.77, # add treatment costs assume 77% of treatment costs are from the public sector

    ## severe treatment costs
    cost_severe =
    # adult treatment component
      ((o5_severe * treatment * SEVcost_adult) +
        # child treatment component
        (u5_severe * treatment * SEVcost_child)) * 0.77, # add treatment costs assume 77% of treatment costs are from the public sector

    ## SMC costs
    cost_smc = n_smc_treated * SMCcost,

    ## Vaccine costs
    cost_vaccine =
      (dose1 + dose2 + dose3 + dose4 + dose5 + dose6 + dose7) *
        (rtss_cost_per_dose + delivery_cost),

    ## TOTAL costs
    cost_total = cost_clinical + cost_severe + cost_smc + cost_vaccine,

    ## cost just among children
    cost_total_u5 =
    # cost clinical
      ((u5_cases - u5_severe) * treatment * TREATcost_child) * 0.77 + # add treatment costs assume 77% of treatment costs are from the public sector
        # cost severe
        (u5_severe * treatment * SEVcost_child) * 0.77 + # add treatment costs assume 77% of treatment costs are from the public sector
        # cost SMC and cost VAX are all among children
        cost_smc + cost_vaccine,

    ## DISCOUNTED COSTS
    ## uncomplicated costs
    cost_clinical_discounted =
    # adult treatment component
      (((o5_cases_discounted - o5_severe_discounted) * treatment * TREATcost_adult) +
        # child treatment component
        ((u5_cases_discounted - u5_severe_discounted) * treatment * TREATcost_child)) * 0.77, # add treatment costs assume 77% of treatment costs are from the public sector

    ## severe treatment costs
    cost_severe_discounted =
    # adult treatment component
      ((o5_severe_discounted * treatment * SEVcost_adult) +
        # child treatment component
        (u5_severe_discounted * treatment * SEVcost_child)) * 0.77, # add treatment costs assume 77% of treatment costs are from the public sector

    ## SMC costs
    cost_smc_discounted = n_smc_treated_discounted * SMCcost,

    ## Vaccine costs
    cost_vaccine_discounted =
      (dose1_discounted + dose2_discounted + dose3_discounted + dose4_discounted + dose5_discounted + dose6_discounted + dose7_discounted) *
        (rtss_cost_per_dose + delivery_cost),

    ## TOTAL costs
    cost_total_discounted = cost_clinical_discounted + cost_severe_discounted + cost_smc_discounted + cost_vaccine_discounted,

    ## cost just among children
    cost_total_u5_discounted =
    # cost clinical
      ((u5_cases_discounted - u5_severe_discounted) * treatment * TREATcost_child) * 0.77 + # add treatment costs assume 77% of treatment costs are from the public sector
        # cost severe
        (u5_severe_discounted * treatment * SEVcost_child) * 0.77 + # add treatment costs assume 77% of treatment costs are from the public sector
        # cost SMC and cost VAX are all among children
        cost_smc_discounted + cost_vaccine_discounted
  )

##-Isolate baseline scenarios with costs------------------------------
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

##-Isolate vaccine scenarios and calculate ICERs-----------------------------
cost_output_vaccine <-
  cost_output |>
  filter(vaccine_model != "no-vaccine") |>
  left_join(
    cost_output_baseline,
    by = c("drawID", "EIR", "warmup", "sim_length", "population", "pfpr", "seasonality", "speciesprop", "treatment", "SMC")
  ) |>
  # calculate cost effectiveness
  mutate(
    cases_averted = (o5_cases_baseline + u5_cases_baseline) - (o5_cases + u5_cases),
    dalys_averted = (o5_dalys_baseline + u5_dalys_baseline) - (o5_dalys + u5_dalys),
    incremental_cost = (cost_total - cost_total_baseline),
    cost_saved = (cost_clinical_baseline + cost_severe_baseline) - (cost_clinical + cost_severe),
    cost_vax_diff = cost_vaccine - cost_vaccine_baseline,
    cost_effectiveness_daly = (cost_total - cost_total_baseline) / ((o5_dalys_baseline + u5_dalys_baseline) - (o5_dalys + u5_dalys)),
    cost_effectiveness_case = (cost_total - cost_total_baseline) / ((o5_cases_baseline + u5_cases_baseline) - (o5_cases + u5_cases)),
    cost_effectiveness_daly_u5 = (cost_total - cost_total_baseline) / (u5_dalys_baseline - u5_dalys),
    cost_effectiveness_case_u5 = (cost_total - cost_total_baseline) / (u5_cases_baseline - u5_cases),

    # discounted outcomes
    cases_averted_discounted = (o5_cases_discounted_baseline + u5_cases_discounted_baseline) - (o5_cases_discounted + u5_cases_discounted),
    dalys_averted_discounted = (o5_dalys_discounted_baseline + u5_dalys_discounted_baseline) - (o5_dalys_discounted + u5_dalys_discounted),
    incremental_cost_discounted = (cost_total_discounted - cost_total_discounted_baseline),
    cost_saved_discounted = (cost_clinical_discounted_baseline + cost_severe_discounted_baseline) - (cost_clinical_discounted + cost_severe_discounted),
    cost_vax_diff_discounted = cost_vaccine_discounted - cost_vaccine_discounted_baseline,
    cost_effectiveness_daly_discounted = (cost_total_discounted - cost_total_discounted_baseline) / ((o5_dalys_discounted_baseline + u5_dalys_discounted_baseline) - (o5_dalys_discounted + u5_dalys_discounted)),
    cost_effectiveness_case_discounted = (cost_total_discounted - cost_total_discounted_baseline) / ((o5_cases_discounted_baseline + u5_cases_discounted_baseline) - (o5_cases_discounted + u5_cases_discounted)),
    cost_effectiveness_daly_u5_discounted = (cost_total_discounted - cost_total_discounted_baseline) / (u5_dalys_discounted_baseline - u5_dalys_discounted),
    cost_effectiveness_case_u5_discounted = (cost_total_discounted - cost_total_discounted_baseline) / (u5_cases_discounted_baseline - u5_cases_discounted)
  )

# save the cost output
saveRDS(cost_output_vaccine, "./03.outputs/paper-results/cost_outputs.rds")

##-SUBANALYSIS WITH 5-DOSE AS THE BASELINE-----------------------------------

# isolate baseline scenarios with costs
cost_output_baseline_5dose <-
  cost_output |>
  filter(dosing_assumption == "5-doses") |>
  select(
    drawID, EIR, warmup, sim_length, population, pfpr,
    seasonality, speciesprop, treatment, SMC, RTSS,
    coverage_assumption, vaccine_model, rtss_cost_per_dose,
    seasonal_ages,
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

# isolate vaccine scenarios with costs
cost_output_vaccine_5dose_baseline <-
  cost_output |>
  filter(dosing_assumption == "7-doses") |>
  inner_join(
    cost_output_baseline_5dose,
    by = c(
      "drawID", "EIR", "warmup", "sim_length", "population", "pfpr", "seasonality", "speciesprop", "treatment", "SMC", "RTSS", "coverage_assumption",
      "vaccine_model", "rtss_cost_per_dose", "seasonal_ages"
    )
  ) |>
  # calculate cost effectiveness
  mutate(
    cases_averted = (o5_cases_baseline + u5_cases_baseline) - (o5_cases + u5_cases),
    dalys_averted = (o5_dalys_baseline + u5_dalys_baseline) - (o5_dalys + u5_dalys),
    incremental_cost = (cost_total - cost_total_baseline),
    cost_saved = (cost_clinical_baseline + cost_severe_baseline) - (cost_clinical + cost_severe),
    cost_vax_diff = cost_vaccine - cost_vaccine_baseline,
    cost_effectiveness_daly = (cost_total - cost_total_baseline) / ((o5_dalys_baseline + u5_dalys_baseline) - (o5_dalys + u5_dalys)),
    cost_effectiveness_case = (cost_total - cost_total_baseline) / ((o5_cases_baseline + u5_cases_baseline) - (o5_cases + u5_cases)),
    cost_effectiveness_daly_u5 = (cost_total - cost_total_baseline) / (u5_dalys_baseline - u5_dalys),
    cost_effectiveness_case_u5 = (cost_total - cost_total_baseline) / (u5_cases_baseline - u5_cases),

    # discounted outcomes
    cases_averted_discounted = (o5_cases_discounted_baseline + u5_cases_discounted_baseline) - (o5_cases_discounted + u5_cases_discounted),
    dalys_averted_discounted = (o5_dalys_discounted_baseline + u5_dalys_discounted_baseline) - (o5_dalys_discounted + u5_dalys_discounted),
    incremental_cost_discounted = (cost_total_discounted - cost_total_discounted_baseline),
    cost_saved_discounted = (cost_clinical_discounted_baseline + cost_severe_discounted_baseline) - (cost_clinical_discounted + cost_severe_discounted),
    cost_vax_diff_discounted = cost_vaccine_discounted - cost_vaccine_discounted_baseline,
    cost_effectiveness_daly_discounted = (cost_total_discounted - cost_total_discounted_baseline) / ((o5_dalys_discounted_baseline + u5_dalys_discounted_baseline) - (o5_dalys_discounted + u5_dalys_discounted)),
    cost_effectiveness_case_discounted = (cost_total_discounted - cost_total_discounted_baseline) / ((o5_cases_discounted_baseline + u5_cases_discounted_baseline) - (o5_cases_discounted + u5_cases_discounted)),
    cost_effectiveness_daly_u5_discounted = (cost_total_discounted - cost_total_discounted_baseline) / (u5_dalys_discounted_baseline - u5_dalys_discounted),
    cost_effectiveness_case_u5_discounted = (cost_total_discounted - cost_total_discounted_baseline) / (u5_cases_discounted_baseline - u5_cases_discounted)
  )


# save the cost output
saveRDS(cost_output_vaccine_5dose_baseline, "./03.outputs/paper-results/cost_outputs_5dose_baseline.rds")

#-AGE DISSAGGREGATED IMPACT---------------------------------------------------------------------------
## -Add epi indicators and define age groups--------------------
processed_out_ab <-
  raw_out |>
  mortality_rate() |>
  outcome_uncertainty() |>
  daly_components() |>
  # remove pre-combined age outputs and have just single age bands
  ungroup() |>
  mutate(age_category = case_when(
    age_lower == 0 & age_upper == 1 ~ "ab",
    age_lower == 1 & age_upper == 2 ~ "ab",
    age_lower == 2 & age_upper == 3 ~ "ab",
    age_lower == 3 & age_upper == 4 ~ "ab",
    age_lower == 4 & age_upper == 5 ~ "ab",
    age_lower == 5 & age_upper == 6 ~ "ab",
    age_lower == 6 & age_upper == 7 ~ "ab",
    age_lower == 7 & age_upper == 8 ~ "ab",
    age_lower == 8 & age_upper == 9 ~ "ab",
    age_lower == 9 & age_upper == 10 ~ "ab",
    age_lower == 10 & age_upper == 11 ~ "ab",
    age_lower == 11 & age_upper == 12 ~ "ab",
    age_lower == 12 & age_upper == 13 ~ "ab",
    age_lower == 13 & age_upper == 14 ~ "ab",
    age_lower == 14 & age_upper == 15 ~ "ab",
    age_lower == 15 & age_upper == 16 ~ "ab",
    age_lower == 16 & age_upper == 17 ~ "ab",
    age_lower == 17 & age_upper == 18 ~ "ab",
    age_lower == 18 & age_upper == 19 ~ "ab",
    age_lower == 19 & age_upper == 20 ~ "ab",
  )) |>
  filter(age_category %in% c("ab")) |>
  select(-age) |>
  group_by(
    ID, runname, drawID, EIR, warmup, sim_length, population, pfpr,
    seasonality, speciesprop, treatment, SMC, RTSS, dosing_assumption,
    coverage_assumption, seasonal_ages, vaccine_model, age_category, age_lower, age_upper
  ) |>
  # sum all other variables over ages in the <10 category
  distinct() |>
  rename(n_smc_treated_discounted = n_smc_treeated_discounted)

##-Identify baseline runs-------------------------------------
baseline_data_ab <-
  processed_out_ab |>
  ungroup() |>
  # filter(age_category == "0-5") |>
  filter(vaccine_model == "no-vaccine") |>
  rename(
    n_treated_baseline = n_treated,
    dose1_baseline = dose1,
    dose2_baseline = dose2,
    dose3_baseline = dose3,
    dose4_baseline = dose4,
    dose5_baseline = dose5,
    dose6_baseline = dose6,
    dose7_baseline = dose7,
    n_smc_treated_baseline = n_smc_treated,
    n_baseline = n,
    n_treated_discounted_baseline = n_treated_discounted,
    dose1_discounted_baseline = dose1_discounted,
    dose2_discounted_baseline = dose2_discounted,
    dose3_discounted_baseline = dose3_discounted,
    dose4_discounted_baseline = dose4_discounted,
    dose5_discounted_baseline = dose5_discounted,
    dose6_discounted_baseline = dose6_discounted,
    dose7_discounted_baseline = dose7_discounted,
    n_smc_treated_discounted_baseline = n_smc_treated_discounted,
    inc_baseline = inc,
    sev_baseline = sev,
    inc_discounted_baseline = inc_discounted,
    sev_discounted_baseline = sev_discounted,
    cases_baseline = cases,
    severe_cases_baseline = severe_cases,
    cases_discounted_baseline = cases_discounted,
    severe_cases_discounted_baseline = severe_cases_discounted,
    mortality_rate_baseline = mortality_rate,
    deaths_baseline = deaths,
    deaths_discounted_baseline = deaths_discounted,
    cases_lower_baseline = cases_lower,
    cases_upper_baseline = cases_upper,
    deaths_lower_baseline = deaths_lower,
    deaths_upper_baseline = deaths_upper,
    cases_discounted_lower_baseline = cases_discounted_lower,
    cases_discounted_upper_baseline = cases_discounted_upper,
    deaths_discounted_lower_baseline = deaths_discounted_lower,
    deaths_discounted_upper_baseline = deaths_discounted_upper,
    yll_baseline = yll,
    yll_lower_baseline = yll_lower,
    yll_upper_baseline = yll_upper,
    yld_baseline = yld,
    yld_lower_baseline = yld_lower,
    yld_upper_baseline = yld_upper,
    yll_discounted_baseline = yll_discounted,
    yll_discounted_lower_baseline = yll_discounted_lower,
    yll_discounted_upper_baseline = yll_discounted_upper,
    yld_discounted_baseline = yld_discounted,
    yld_discounted_lower_baseline = yld_discounted_lower,
    yld_discounted_upper_baseline = yld_discounted_upper,
    daly_baseline = daly,
    daly_lower_baseline = daly_lower,
    daly_upper_baseline = daly_upper,
    daly_discounted_baseline = daly_discounted,
    daly_discounted_lower_baseline = daly_discounted_lower,
    daly_discounted_upper_baseline = daly_discounted_upper
  ) |>
  distinct() |>
  select(
    -vaccine_model, -coverage_assumption, -RTSS, -dosing_assumption,
    -seasonal_ages, -ID, -runname
  ) |>
  distinct()

## -Calculate events averted------------------------------------------
averted_outputs_ab <-
  processed_out_ab |>
  ungroup() |>
  filter(vaccine_model != "no-vaccine") |>
  left_join(
    baseline_data_ab,
    by = c(
      "drawID", "EIR", "warmup", "sim_length", "population", "pfpr", "seasonality", "speciesprop", "treatment", "SMC",
      "age_category", "age_lower", "age_upper"
    )
  ) |>
  # Processing the outputs to calculate cases, deaths, and DALYs averted
  mutate(
    # model outputted cases averted
    cases_averted = cases_baseline - cases,
    cases_averted_discounted = cases_discounted_baseline - cases_discounted,
    # model outputted severe cases averted
    severe_cases_averted = severe_cases_baseline - severe_cases,
    severe_cases_averted_discounted = severe_cases_discounted_baseline - severe_cases_discounted,
    # model outputted deaths averted
    deaths_averted = deaths_baseline - deaths,
    deaths_averted_discounted = deaths_discounted_baseline - deaths_discounted,
    # model outputted DALYs averted
    dalys_averted = daly_baseline - daly,
    dalys_averted_discounted = daly_discounted_baseline - daly_discounted,
    # processed cases averted per 100,000 children
    cases_averted_per_100000 = ((cases_baseline / n_baseline) - (cases / n)) * 100000,
    cases_averted_per_100000_discounted = ((cases_discounted_baseline / n_baseline) - (cases_discounted / n)) * 100000,
    # processed severe cases averted per 100,000 children
    severe_cases_averted_per_100000 = ((severe_cases_baseline / n_baseline) - (severe_cases / n)) * 100000,
    severe_cases_averted_per_100000_discounted = ((severe_cases_discounted_baseline / n_baseline) - (severe_cases_discounted / n)) * 100000,
    # processed deaths averted per 100,000 children
    deaths_averted_per_100000 = ((deaths_baseline / n_baseline) - (deaths / n)) * 100000,
    deaths_averted_per_100000_discounted = ((deaths_discounted_baseline / n_baseline) - (deaths_discounted / n)) * 100000,
    # processed DALYs averted per 100,000 children
    dalys_averted_per_100000 = ((daly_baseline / n_baseline) - (daly / n)) * 100000,
    dalys_averted_per_100000_discounted = ((daly_discounted_baseline / n_baseline) - (daly_discounted / n)) * 100000
  ) |>
  ungroup()

# save
saveRDS(averted_outputs_ab, "./03.outputs/paper-results/averted_outputs_age-based.rds")

## -Incremental impact in an age-based setting-------------------------
incremental_outputs_ab <-
  averted_outputs_ab |>
  filter(dosing_assumption %in% c("5-doses", "7-doses")) |>
  mutate(
    dosing_num = case_when(
      dosing_assumption == "5-doses" ~ 5,
      dosing_assumption == "7-doses" ~ 7
    )
  ) |>
  inner_join(
    averted_outputs_ab |>
      mutate(
        dosing_num = case_when(
          dosing_assumption == "5-doses" ~ 5,
          dosing_assumption == "7-doses" ~ 7
        )
      ),
    by = c(
      "RTSS", "seasonality", "pfpr", "SMC", "vaccine_model", "drawID", "seasonal_ages", "coverage_assumption",
      "age_category", "age_lower", "age_upper"
    ),
    suffix = c("_from", "_to")
  ) |>
  filter(dosing_num_to > dosing_num_from) |> # only compare to higher doses
  mutate(
    # CASES
    incremental_cases = cases_averted_to - cases_averted_from,
    incremental_percent_cases = (incremental_cases / cases_averted_from) * 100,

    # SEVERE CASES
    incremental_severe_cases = severe_cases_averted_to - severe_cases_averted_from,
    incremental_percent_severe_cases = (incremental_severe_cases / severe_cases_averted_from) * 100,

    # DEATHS
    incremental_deaths = deaths_averted_to - deaths_averted_from,
    incremental_percent_deaths = (incremental_deaths / deaths_averted_from) * 100,

    # DALYS
    incremental_dalys = dalys_averted_to - dalys_averted_from,
    incremental_percent_dalys = (incremental_dalys / dalys_averted_from) * 100
  ) |>
  distinct()

# save
saveRDS(incremental_outputs_ab, "./03.outputs/paper-results/incremental_outputs_age-based.rds")
