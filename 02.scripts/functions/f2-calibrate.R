# ------------------------------------------------------------------------------
# Purpose: Calibrate starting EIR values to match target PfPR for baseline scenarios.
# Inputs: Baseline parameter table and scenario/draw indices.
# Outputs: Calibration outputs in 03.outputs/pf-eir-match for each indexed run.
# Dependencies: cali and malariasimulation dependencies loaded upstream.
# Run Stage: Function library called by s2-baseline-calibration.R.
# Notes: Comment-only documentation update; logic and outputs unchanged.
# ------------------------------------------------------------------------------


# calibration ref: https://mrc-ide.github.io/cali/articles/Basic_calibration.html

# Function to run calibration for a given scenario and parameter draw
# Function: PRmatch
# Purpose: Calibrate starting EIR to match target PfPR for one scenario/draw.
# Inputs: x (scenario index), y (parameter draw).
# Outputs: Returns and saves PfPR-EIR calibration output for one run.
# Assumptions: Assumes baseline parameter file exists and calibration package is available.

PRmatch <- function(x, y){ # x = scenario # , y = parameter draw #

  #-read in selected scenario--------------------------------
  data <- readRDS("./01.data/baseline/baseline_parameters.rds")[x,]

  #-unlist parameter set from baseline scenarios------------
  p <- unlist(data$params, recursive = F)

  #-define target: PfPR2-10 value--------------------------
  target <- data$pfpr

  # function for calibration-------------------------------
  # average PfPR - years 4-6 of the simulation
  # 6 years is enough to match by PfPR
  # 21 years is needed to match by c inc
  year <- 365
  p$timesteps <- 9 * year   # simulation run time = 9 years

  summary_mean_pfpr_2_10_6y9y <- function(x){

    x$year <- ceiling(x$timestep / year)
    x <- x |> filter(year >= 7)
    prev_2_10 <- mean(x$n_detect_pcr_730_3650 / x$n_age_730_3650)
    return(prev_2_10)

  }

  #-run calibration model--------------------------------
  set.seed(123)
  out <- cali::calibrate(
    parameters = p,
    target = target,
    summary_function = summary_mean_pfpr_2_10_6y9y,
    eq_prevalence = target[1],
    eq_ft = data$treatment,
    eir_limits = c(0.001, 1500),
    max_attempts = 40,
    human_population = c(5000, 50000, 100000)
    )

  #-store init_EIR results-----------------------------
  #.rds file to be read in later
  PR <- data.frame(scenarioID = x, drawID = y)
  PR$starting_EIR <- out
  PR$ID <- data$ID

  saveRDS(PR, paste0('./03.outputs/pf-eir-match/PRmatch_draws_', data$ID, '.rds'))

}

