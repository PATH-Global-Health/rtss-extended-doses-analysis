# ------------------------------------------------------------------------------
# Purpose: Run one configured simulation scenario and persist model outputs.
# Inputs: Scenario index with run parameters and calibration match table.
# Outputs: One scenario-level simulation output RDS in individual-model-sim-outputs.
# Dependencies: malariasimulation runtime plus prepared parameter and calibration files.
# Run Stage: Function library called by s3-model-runs.R.
# ------------------------------------------------------------------------------

# Set interventions and run malariasimulation ----------------------------------

# Function: runsim
# Purpose: Run one malaria simulation scenario and write model output.
# Inputs: x (scenario index).
# Outputs: Writes one scenario output RDS with metadata columns.
# Assumptions: Assumes calibration outputs and run parameters are available.

runsim <- function(x){ # x = scenario #

  year <- 365
  month <- year / 12

  # read in selected scenario
  data <- readRDS("./01.data/vaccine-scenarios/run_parameters.rds")[x,]
  match <- readRDS("./03.outputs/pf-eir-match/EIRestimates.rds") |> select(-scenarioID)

  # EIR / prev match from "PfPR_EIR_match.R"
  data <- data |> left_join(match, by = c("drawID", "ID"))

  # EIR equilibrium ----------
  params <- set_equilibrium(unlist(data$params, recursive = F), as.numeric(data$starting_EIR))

  # run simulation ----------
  # Fixed seed is intentional for reproducible and comparable outputs under draw-0 design.
  set.seed(123)

  output <- run_simulation(
    timesteps = data$warmup + data$sim_length,
    # correlations = correlations,
    parameters = params) |>

    # add vars to output
    mutate(ID = data$ID,
           runname = data$run_name,
           drawID = data$drawID,
           EIR = data$starting_EIR,
           warmup = data$warmup,
           sim_length = data$sim_length,
           population = data$population,
           pfpr = data$pfpr,
           timestep = timestep - data$warmup,
           seasonality = data$seas_name,
           speciesprop = paste(data$speciesprop, sep = ",", collapse = ""),
           treatment = data$treatment,
           SMC = data$SMC,
           RTSS = data$RTSS,
           dosing_assumption = data$dosing_assumption,
           coverage_assumption = data$coverage_assumption,
           seasonal_ages = data$seasonal_ages,
           vaccine_model = data$vaccine_model) |>
    ungroup() |>
    filter(timestep > 0) |> # remove warmup period

    # statistics by month
    mutate(year = ceiling(timestep/year),
           month = ceiling(timestep/month)) |>

    # keep only necessary variables
    dplyr::select(ID,
                  runname,
                  drawID,
                  EIR,
                  warmup,
                  sim_length,
                  population,
                  pfpr,
                  timestep,
                  seasonality,
                  speciesprop,
                  treatment,
                  SMC,
                  RTSS,
                  dosing_assumption,
                  coverage_assumption,
                  seasonal_ages,
                  vaccine_model,
                  year,
                  month,
                  starts_with("n_inc_severe"), starts_with("p_inc_severe"),
                  starts_with("n_pev"),
                  starts_with("n_inc"), starts_with("p_inc"),
                  starts_with("n_detect"), starts_with("p_detect"),
                  starts_with("n_"), -n_bitten, n_treated, n_infections)


  # save output ----------
  saveRDS(output, paste0("./03.outputs/individual-model-sim-outputs/", x, "-", data$run_name,  ".rds"))
}
