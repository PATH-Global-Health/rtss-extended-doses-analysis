# ------------------------------------------------------------------------------
# Purpose: Add mortality and DALY calculations to processed epidemiologic outputs.
# Inputs: Processed scenario-level epidemiologic data frame(s).
# Outputs: Data frames augmented with deaths, YLL, YLD, and DALY metrics.
# Dependencies: Reference assumptions and helper calculations defined in this file.
# Run Stage: Function library called by s4-processing-model-runs.R.
# ------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Deaths and DALYs addition to outputs
# 1. Read in processed outputs from model simulations
# 2. Add in calculations for deaths and DALYs using standard methods and assumptions
# 3. Save outputs with additional variables for use in paper analyses
#-------------------------------------------------------------------------------

# input: data frame
# process: reads in raw HPC output, adds R21 doses, condenses output over simulation length
# output: data frame with additional variables for deaths and DALYs


# DALYs = Years of life lost (YLL) + Years of life with disease (YLD)
# YLL = Deaths * remaining years of life
# YLD = cases and severe cases * disability weighting  * episode_length
# CE = $ per event (case, death DALY) averted

# reference code: https://github.com/mrc-ide/gf/blob/69910e798a2ddce240c238d291bc36ea40661b90/R/epi.R
# weights from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4772264/ {Gunda _et al_, 2016}


# mortality --------------------------------------------------------------------

# Function: mortality_rate
# Purpose: Add mortality-rate based death estimates to processed outputs.
# Inputs: x (data frame), scaler and mortality parameters.
# Outputs: Data frame with death estimates and uncertainty bounds.
# Assumptions: Assumes expected epidemiologic columns are present.

mortality_rate <- function(x,
                           scaler = 0.215, # severe case to death scaler
                           treatment_scaler = 0.578) { # treatment modifier (range 0.157-0.719); https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1003359
  x |>

    # mortality rate
    dplyr::mutate(mortality_rate = (1 - (treatment_scaler * .data$treatment)) * scaler * .data$sev) |>

    dplyr::mutate(deaths = .data$mortality_rate * .data$n) |> # deaths

    # apply to discounted severe cases
    dplyr::mutate(deaths_discounted = (1 - (treatment_scaler * .data$treatment)) * scaler * .data$severe_cases_discounted)

}


# case and death uncertainty ---------------------------------------------------

# Function: outcome_uncertainty
# Purpose: Attach uncertainty intervals around modeled outcomes.
# Inputs: x (data frame) plus lower/upper multipliers.
# Outputs: Data frame with lower/upper uncertainty columns.
# Assumptions: Assumes baseline outcome columns are already computed.

outcome_uncertainty <- function(x,
                                cases_cv = 0.227,   # case uncertainty SD scaler
                                deaths_cv = 0.265){ # death uncertainty SD scaler
  x |>
    dplyr::mutate(cases_lower = round(pmax(0, stats::qnorm(0.025, .data$cases, .data$cases * cases_cv))),
                  cases_upper = round(stats::qnorm(0.975, .data$cases, .data$cases * cases_cv)),
                  deaths_lower = round(pmax(0, stats::qnorm(0.025, .data$deaths, .data$deaths * deaths_cv))),
                  deaths_upper = round(stats::qnorm(0.975, .data$deaths, .data$deaths * deaths_cv)),
                  # discounted outcomes
                  cases_discounted_lower = round(pmax(0, stats::qnorm(0.025, .data$cases_discounted,
                                                                      .data$cases_discounted * cases_cv))),
                  cases_discounted_upper = round(stats::qnorm(0.975, .data$cases_discounted,
                                                              .data$cases_discounted * cases_cv)),
                  deaths_discounted_lower = round(pmax(0, stats::qnorm(0.025, .data$deaths_discounted,
                                                                       .data$deaths_discounted * deaths_cv))),
                  deaths_discounted_upper = round(stats::qnorm(0.975, .data$deaths_discounted,
                                                               .data$deaths_discounted * deaths_cv)))
}


# DALYs ------------------------------------------------------------------------
# weights are an approximation of YLD estimates from https://ghdx.healthdata.org/record/ihme-data/gbd-2017-disability-weights; disabilities from comorbid conditions such as motor impairment and anemia are excluded.
# Life expectancy at birth (years) [Data table]. Geneva, 2020 https://www.who.int/data/gho/data/indicators/indicator-details/GHO/life-expectancy-at-birth-(years) (accessed Nov, 2025).
# Function: daly_components
# Purpose: Compute YLD, YLL, and DALY components and discounted values.
# Inputs: x (data frame) with cases and deaths.
# Outputs: Data frame with DALY component columns added.
# Assumptions: Assumes mortality and outcome fields already exist.

daly_components <- function(x,
                            lifespan = 63.6,         # average life expectancy
                            episode_length = 0.01375, # average length of clinical episode
                            severe_episode_length = 0.04795, # average length of severe episode
                            treatment_coverage = 0.45, # treatment coverage
                            treatment_scaler = 0.578, # treatment scaler
                            mild_dw = 0.006, # disability weight mild malaria
                            moderate_dw = 0.051, # disability moderate severe malaria
                            severe_dw = 0.133) { # disability weight severe malaria
  output <- x |>
    dplyr::mutate(yll = .data$deaths * (lifespan - ((.data$age_lower + .data$age_upper) / 2)),
                  yll_lower = .data$deaths_lower * (lifespan - ((.data$age_lower + .data$age_upper) / 2)),
                  yll_upper = .data$deaths_upper * (lifespan - ((.data$age_lower + .data$age_upper) / 2)),

                  yll = ifelse(yll < 0, 0, yll),                    # should be no negative yll from older age groups
                  yll_lower = ifelse(yll_lower < 0, 0, yll_lower),  # should be no negative yll from older age groups
                  yll_upper =  ifelse(yll_upper < 0, 0, yll_upper), # should be no negative yll from older age groups

                  yld = dplyr::case_when(.data$age_upper < 5 ~ .data$cases * episode_length * moderate_dw + .data$severe_cases * severe_episode_length * severe_dw * (treatment_scaler * treatment_coverage + (1 - treatment_coverage)),

                                         .data$age_upper >= 5 ~ .data$cases * episode_length * mild_dw + .data$severe_cases * severe_episode_length * severe_dw * (treatment_scaler * treatment_coverage + (1 - treatment_coverage))),

                  yld_lower = dplyr::case_when(.data$age_upper < 5 ~ .data$cases_lower * episode_length * moderate_dw + .data$severe_cases * severe_episode_length * severe_dw * (treatment_scaler * treatment_coverage + (1 - treatment_coverage)),

                                               .data$age_upper >= 5 ~ .data$cases_lower * episode_length * mild_dw + .data$severe_cases * severe_episode_length * severe_dw * (treatment_scaler * treatment_coverage + (1 - treatment_coverage))),

                  yld_upper = dplyr::case_when(.data$age_upper < 5 ~ .data$cases_upper * episode_length * moderate_dw + .data$severe_cases * severe_episode_length * severe_dw * (treatment_scaler * treatment_coverage + (1 - treatment_coverage)),

                                               .data$age_upper >= 5 ~ .data$cases_upper * episode_length * mild_dw + .data$severe_cases * severe_episode_length * severe_dw * (treatment_scaler * treatment_coverage + (1 - treatment_coverage))),

                  # discounted outcomes

                  yll_discounted = .data$deaths_discounted * (lifespan - ((.data$age_lower + .data$age_upper) / 2)),
                  yll_discounted_lower = .data$deaths_discounted_lower * (lifespan - ((.data$age_lower + .data$age_upper) / 2)),
                  yll_discounted_upper = .data$deaths_discounted_upper * (lifespan - ((.data$age_lower + .data$age_upper) / 2)),

                  yll_discounted = ifelse(yll_discounted < 0, 0, yll_discounted),                    # should be no negative yll from older age groups
                  yll_discounted_lower = ifelse(yll_discounted_lower < 0, 0, yll_discounted_lower),  # should be no negative yll from older age groups
                  yll_discounted_upper =  ifelse(yll_discounted_upper < 0, 0, yll_discounted_upper), # should be no negative yll from older age groups

                  yld_discounted = dplyr::case_when(.data$age_upper < 5 ~ .data$cases_discounted * episode_length * moderate_dw + .data$severe_cases_discounted * severe_episode_length * severe_dw * (treatment_scaler * treatment_coverage + (1 - treatment_coverage)),

                                                    .data$age_upper >= 5 ~ .data$cases_discounted * episode_length * mild_dw + .data$severe_cases_discounted * severe_episode_length * severe_dw * (treatment_scaler * treatment_coverage + (1 - treatment_coverage))),

                  yld_discounted_lower = dplyr::case_when(.data$age_upper < 5 ~ .data$cases_discounted_lower * episode_length * moderate_dw + .data$severe_cases_discounted * severe_episode_length * severe_dw * (treatment_scaler * treatment_coverage + (1 - treatment_coverage)),

                                                          .data$age_upper >= 5 ~ .data$cases_discounted_lower * episode_length * mild_dw + .data$severe_cases_discounted * severe_episode_length * severe_dw * (treatment_scaler * treatment_coverage + (1 - treatment_coverage))),

                  yld_discounted_upper = dplyr::case_when(.data$age_upper < 5 ~ .data$cases_discounted_upper * episode_length * moderate_dw + .data$severe_cases_discounted * severe_episode_length * severe_dw * (treatment_scaler * treatment_coverage + (1 - treatment_coverage)),

                                                          .data$age_upper >= 5 ~ .data$cases_discounted_upper * episode_length * mild_dw + .data$severe_cases_discounted * severe_episode_length * severe_dw * (treatment_scaler * treatment_coverage + (1 - treatment_coverage)))
    ) |>

    dplyr::mutate(daly = yll + yld,
                  daly_upper = yll_lower + yld_lower,
                  daly_lower = yll_upper + yld_upper,
                  daly_discounted = yll_discounted + yld_discounted,
                  daly_discounted_upper = yll_discounted_lower + yld_discounted_lower,
                  daly_discounted_lower = yll_discounted_upper + yld_discounted_upper)

  return(output)

}

