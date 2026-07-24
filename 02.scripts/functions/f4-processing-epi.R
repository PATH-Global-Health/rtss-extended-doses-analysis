# ------------------------------------------------------------------------------
# Purpose: Transform raw simulation outputs into processed epidemiologic summaries.
# Inputs: Raw scenario output files plus run identifiers and metadata.
# Outputs: Processed per-scenario RDS outputs for paper result pipelines.
# Dependencies: tidyverse/data.table utilities and model output schema assumptions.
# Run Stage: Function library called by s4-processing-model-runs.R.
# ------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Processing functions
#
# 1. Reduction of raw model outputs into smaller dataframe
# 2. Adding vaccine doses in better format
# 3. Discounting cases and doses
# 4. Summarizing data over the first 15 years of the simulation
# 5. Moving from wide to long age groups
# 6. Calculating incidence and severe incidence
# 7. Saving pre-processed outputs for use in paper analyses
# 
#-------------------------------------------------------------------------------

# input: index of malariasimulation HPC run
# process: reads in raw HPC output, reduces to monthly vars of interest adds doses, condenses output over simulation length
# output: data frame with one row per age group

# Function: processing
# Purpose: Process one raw simulation output into aggregated epidemiologic outputs.
# Inputs: x (scenario index), rn (run name).
# Outputs: Writes processed scenario-level output RDS.
# Assumptions: Assumes raw output schema from runsim().

processing <- function(x, rn){

  # read in specified rds file
   output <- readRDS(paste0("./03.outputs/individual-model-sim-outputs/", format(x, scientific = FALSE), "-", rn, ".rds"))

  # extract type of RTS,S  intervention
  RTSS = output$RTSS[1]
  SMC = output$SMC[1]
  dosing_assumption = output$dosing_assumption[1]

  # add vaccine doses in better format
  output <-
      output |>
      rowwise() |>
      mutate(dose1 = 0,
             dose2 = 0,
             dose3 = 0,
             dose4 = 0,
             dose5 = 0,
             dose6 = 0,
             dose7 = 0) |>
      ungroup()

  if(SMC == 0){
    output <-
      output |>
      rowwise() |>
      mutate(n_smc_treated = 0,
             n_smc_drug_efficacy_failures = 0,
             n_smc_successfully_treated = 0
      ) |>
      ungroup()
  }

  if(RTSS == "Age-based" | RTSS == "Hybrid"){
    output <-
      output |>
      rowwise() |>
      mutate(dose1 = n_pev_epi_dose_1,
             dose2 = n_pev_epi_dose_2,
             dose3 = n_pev_epi_dose_3) |>
      ungroup()
  }

  if(RTSS == "Seasonal" ){
    output <-
      output |>
      rowwise() |>
      mutate(dose1 = n_pev_mass_dose_1,
             dose2 = n_pev_mass_dose_2,
             dose3 = n_pev_mass_dose_3) |>
      ungroup()
  }

  # additional doses
  if(dosing_assumption == "4-doses"){
    if(RTSS == "Age-based" | RTSS == "Hybrid"){
      output <-
        output |>
        rowwise() |>
        mutate(dose4 = n_pev_epi_booster_1) |>
        ungroup()
    }

    if(RTSS == "Seasonal"){
      output <-
        output |>
        rowwise() |>
        mutate(dose4 = n_pev_mass_booster_1) |>
        ungroup()
    }
  }

  if(dosing_assumption == "5-doses"){
    if(RTSS == "Age-based" | RTSS == "Hybrid"){
      output <-
        output |>
        rowwise() |>
        mutate(
          dose4 = n_pev_epi_booster_1,
          dose5 = n_pev_epi_booster_2
          ) |>
        ungroup()
    }

    if(RTSS == "Seasonal"){
      output <-
        output |>
        rowwise() |>
        mutate(
          dose4 = n_pev_mass_booster_1,
          dose5 = n_pev_mass_booster_2
          ) |>
        ungroup()
    }
  }

  if(dosing_assumption == "6-doses"){
    if(RTSS == "Age-based" | RTSS == "Hybrid"){
      output <-
        output |>
        rowwise() |>
        mutate(
          dose4 = n_pev_epi_booster_1,
          dose5 = n_pev_epi_booster_2,
          dose6 = n_pev_epi_booster_3
          ) |>
        ungroup()
    }

    if(RTSS == "Seasonal"){
      output <-
        output |>
        rowwise() |>
        mutate(
          dose4 = n_pev_mass_booster_1,
          dose5 = n_pev_mass_booster_2,
          dose6 = n_pev_mass_booster_3
          ) |>
        ungroup()
    }
  }

  if(dosing_assumption == "7-doses"){
    if(RTSS == "Age-based" | RTSS == "Hybrid"){
      output <-
        output |>
        rowwise() |>
        mutate(
          dose4 = n_pev_epi_booster_1,
          dose5 = n_pev_epi_booster_2,
          dose6 = n_pev_epi_booster_3,
          dose7 = n_pev_epi_booster_4
          ) |>
        ungroup()
    }

    if(RTSS == "Seasonal"){
      output <-
        output |>
        rowwise() |>
        mutate(
          dose4 = n_pev_mass_booster_1,
          dose5 = n_pev_mass_booster_2,
          dose6 = n_pev_mass_booster_3,
          dose7 = n_pev_mass_booster_4
          ) |>
        ungroup()
    }
  }

  # summarize data over the first 15 years 
  sim_length <- output$sim_length[1] / 365

  ## DISCOUNTING
  discount_function <- function(x, rate, year) {
    x * (1/(1+rate)^(year-1))
  }

  # Discount cases and severe cases, then take sum over 15 years
  discounted_output <-
    output |>
    ungroup() |>
    # filter to first 15 years
    filter(
      year <= sim_length
      ) |>
    # apply discounting to incicdence variables and vaccine dose variables
    mutate_at(
      vars(contains("inc"), dose1:dose7),
        ~ discount_function(.x, rate = 0.03, year = year)
      ) |>
    # apply discounting to SMC treated
    mutate_at(
      vars(contains("smc_treated")),
      ~ discount_function(.x, rate = 0.03, year = year)
      ) |>
    # apply discounting to number treated
    mutate_at(
      vars(contains("n_treated")),
      ~ discount_function(.x, rate = 0.03, year = year)
      ) |>
    # take the average value of population data
    mutate_at(
      vars(contains("n_age")),
      mean,
      na.rm = TRUE
      ) |>
    # sum the incidence and vaccine dose variables and smc treated
    mutate_at(
      vars(contains("inc"), dose1:dose7, contains("smc_treated"),
           n_treated),
      sum, na.rm = TRUE
      ) |>
    # round values to whole numbers
    mutate_at(
      vars(contains("inc"), dose1:dose7, contains("smc_treated"),
            contains("n_age"), n_treated),
      round, 0
      ) |>
    # remove year and month values
    select(
     -month, -year, -contains("pev"),
     -n_successfully_treated, -n_infections, -n_detect_lm_730_3650, -p_detect_lm_730_3650,
     -n_drug_efficacy_failures, -n_detect_pcr_730_3650, -timestep, -starts_with("p_"),
     -n_smc_drug_efficacy_failures, -n_smc_successfully_treated
      ) |>
    # keep distinct rows
    distinct() |>

    # moving from wide to long age groups
    pivot_longer(
      cols = c(
        contains("n_age"), contains("n_inc_clinical"), contains("n_inc_severe")
        ),
      names_to = c(
        'age'
        ),
      values_to = c(
        'value'
        )
      ) |>
    mutate(
      n = ifelse(grepl('n_age', age), value, NA),                     # creating var for age group
      inc_clinical = ifelse(grepl('n_inc_clinical', age), value, NA), # creating var for inc_clinical
      inc_severe = ifelse(grepl('n_inc_severe', age), value, NA),     # creating var for inc_severe
      age = gsub('n_inc_clinical_', '', age),                         # combining age vars
      age = gsub('n_inc_severe_', '', age),
      age = gsub('n_age_', '', age)
      ) |>
    group_by(age) |>
    select(-value) |>
    # consolidate
    mutate_at(
      vars(n:inc_severe),
      sum,
      na.rm = TRUE
      ) |>
    distinct() |>
    ungroup()

  # specifiy discounted variables
  discounted_output <-
    discounted_output  |>
    rename(inc_clinical_discounted = inc_clinical,
           inc_severe_discounted = inc_severe,
           dose1_discounted = dose1,
           dose2_discounted = dose2,
           dose3_discounted = dose3,
           dose4_discounted = dose4,
           dose5_discounted = dose5,
           dose6_discounted = dose6,
           dose7_discounted = dose7,
           n_treated_discounted = n_treated,
           n_smc_treeated_discounted = n_smc_treated)

  # Outputs without discounting
  output <-
    output |>
    # first 15 years
    filter(
      year <= sim_length
      ) |>
    # mean of n in each age group
    mutate_at(
      vars(contains("n_age")),
      mean,
      na.rm = TRUE
    ) |>
    # sum the incidence and vaccine dose variables and smc treated
    mutate_at(
      vars(contains("inc"), dose1:dose7, contains("smc_treated"),
           n_treated),
      sum, na.rm = TRUE
    ) |>
    # round values to whole numbers
    mutate_at(
      vars(contains("inc"), dose1:dose7, contains("smc_treated"),
           contains("n_age"), n_treated),
      round, 0
    ) |>
    # remove year and month values
    select(
      -month, -year, -contains("pev"),
      -n_successfully_treated, -n_infections, -n_detect_lm_730_3650, -p_detect_lm_730_3650,
      -n_drug_efficacy_failures, -n_detect_pcr_730_3650, -timestep, -starts_with("p_"),
      -n_smc_drug_efficacy_failures, -n_smc_successfully_treated
    ) |>
    # keep distinct rows
    distinct() |>

    # moving from wide to long age groups
    pivot_longer(
      cols = c(
        contains("n_age"), contains("n_inc_clinical"), contains("n_inc_severe")
      ),
      names_to = c(
        'age'
      ),
      values_to = c(
        'value'
      )
    ) |>
    mutate(
      n = ifelse(grepl('n_age', age), value, NA),                     # creating var for age group
      inc_clinical = ifelse(grepl('n_inc_clinical', age), value, NA), # creating var for inc_clinical
      inc_severe = ifelse(grepl('n_inc_severe', age), value, NA),     # creating var for inc_severe
      age = gsub('n_inc_clinical_', '', age),                         # combining age vars
      age = gsub('n_inc_severe_', '', age),
      age = gsub('n_age_', '', age)
    ) |>
    group_by(age) |>
    select(-value) |>
    # consolidate
    mutate_at(
      vars(n:inc_severe),
      sum,
      na.rm = TRUE
    ) |>
    distinct() |>
    ungroup()

  output <- left_join(output, discounted_output)

  output <- output |>
    separate(col = age, into = c("age_lower", "age_upper"), sep="_", remove = F) |>
    mutate(age_lower = as.numeric(age_lower)/365,
           age_upper = as.numeric(age_upper)/365,
           inc = inc_clinical / n,
           sev = inc_severe / n,
           cases = inc_clinical,
           severe_cases = inc_severe,
           cases_discounted = inc_clinical_discounted,
           severe_cases_discounted = inc_severe_discounted,
           inc_discounted = inc_clinical_discounted / n,
           sev_discounted = inc_severe_discounted / n) |>
    select(-inc_clinical, -inc_severe, -inc_clinical_discounted, -inc_severe_discounted)




  saveRDS(output, paste0("./03.outputs/processed-model-sim-outputs/", x, "-", rn,  ".rds"))

}


