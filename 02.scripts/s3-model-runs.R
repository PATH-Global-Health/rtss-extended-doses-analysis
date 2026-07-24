# ------------------------------------------------------------------------------
# Purpose: Construct intervention scenarios and run model simulations for paper analyses.
# Inputs: Scenario constants, uncertainty draw files, and calibration outputs.
# Outputs: Scenario tables and individual model simulation output files.
# Dependencies: data_and_libraries.R plus functions f1 and f3.
# Run Stage: Pipeline step 3: scenario generation and simulation execution.
# ------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 3-65% pfpr | 2 seasonality profiles | SMC or no SMC
# Hybrid with 5 or 7 dose schedules
# Seasonal RTS,S with 5 or 7 dose schedules
#
# 2 x vaccine efficacy models
# Combination synergy model
# Variety of coverage assumptions
#
# Run for 15 years
#-------------------------------------------------------------------------------

#-libraries---------------------------------------------------------------------
source("./02.scripts/data_and_libraries.R")

#-set up------------------------------------------------------------------------
# make all combinations of baseline scenarios

# defined parameters universal
year <- 365             #year
population <- 1000000   #population size increase because of severe outcomes

# run time
warmup <- 3 * year       # keep the same as in the calibration
sim_length <- 15 * year   # value > 0

# number of parameter draws
# draws <- rbind(draws, data.frame(drawID = 0))
drawID <- readRDS("./01.data/inputs/uncertainty-params/parameter_draws.rds")$draw # sampled in s1_parameter_draws.R
draws  <- tibble(drawID)

# drawID = 0 is intentional: uncertainty is driven by intervention efficacy models,
# while transmission model uncertainty is held fixed due to compute constraints.
draws$drawID <- 0
#-transmission set-up-----------------------------------------------------------
# parasite prevalence 2-10 year olds
pfpr <- c(0.01, 0.03, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.45, 0.55, 0.65) #c(0.1, 0.2, 0.35)

# seasonal profiles: c(g0, g[1], g[2], g[3], h[1], h[2], h[3])
# drawn from mlgts: https://github.com/mrc-ide/mlgts/tree/master/data
# g0 = a0, a = g, b = h
seas_name <- "highly seasonal"
seasonality <- list(c(0.284596,-0.317878,-0.0017527,0.116455,-0.331361,0.293128,-0.0617547))
s1 <- tibble(seasonality, seas_name)

seas_name <- "seasonal"
seasonality <- list(c(0.285505,-0.325352,-0.0109352,0.0779865,-0.132815,0.104675,-0.013919))
s2 <- tibble(seasonality, seas_name)

stable <- bind_rows(s1, s2)

#-Vectors-----------------------------------------------------------------------
# list(arab_params, fun_params, gamb_params)
speciesprop <- data.frame(speciesprop = rbind(list(c(0.25, 0.25, 0.5))),
                          row.names = NULL)
#-Interventions-----------------------------------------------------------------
# treatment coverage
treatment <- c(0.45)

# SMC: 0, 0.75
SMC <- c(0, 0.75)

# RTS,S: none, EPI, SV, hybrid
RTSS <- c("none", "Seasonal", "Hybrid") #"Age-based"

# RTS,S coverage assumptions
# i.	Scenario 1 = 90% coverage no drop out
# ii.	 Scenario 2 = 90% coverage d3, d3-d4 drop out 10% and 10% drop out total over the remaining doses
# iii.	Scenario 3 = 90% coverage d3, 5% drop out each additional dose
# iiii.  Scenario 4 = 90% coverage d3, 10% drop out each additional dose

# Potential for more granular analysis on drop out rates as required
coverage_assumption <- c(
  "S1-optimal-campaigns",
  "S2-enhanced-communication",
  # "S3-sustained-dropout-5perc",
  "S3-sustained-dropout-10perc"
)

dosing_assumption <- c( "5-doses", "7-doses") # remove 4 and 6 as not of interest in WHO comparison "4-doses", "6-doses"

interventions <- crossing(treatment, SMC, RTSS, coverage_assumption, dosing_assumption)

# create combination of all runs
combo <-
  crossing(population, pfpr, stable, warmup, sim_length,
           speciesprop, interventions, draws) |>
  mutate(ID = paste(pfpr, seas_name, drawID, SMC, sep = "_"))

# remove non-applicable scenarios - where vaccine is none drop all the coverage assumptions
combo <-
  combo |>
  mutate(coverage_assumption = case_when(RTSS == "none" ~ "no-vaccine",
                                         TRUE ~ coverage_assumption),
         dosing_assumption = case_when(RTSS == "none" ~ "no-vaccine",
                                       TRUE ~ dosing_assumption)) |>
  distinct()

# add in ID for vaccine model type
vaccine_model <-
  data.frame(
    RTSS = c(
      "Seasonal", "Seasonal", "Seasonal",
      "Hybrid", "Hybrid", "Hybrid",
      "none"
    ),
    vaccine_model = c(
      "model-1", "model-2", "model-3",
      "model-1", "model-2", "model-3",
      "no-vaccine"
    )
  )

# join together
combo <-
  combo |>
  left_join(
    vaccine_model,
    relationship = "many-to-many"
  )

# remove non-applicable scenarios
combo <-
  combo |>
  # remove model 3 when smc is not delivered
  filter(!(SMC == 0 & vaccine_model == "model-3"))

# seasonal age model
seasonal_ages <-
  data.frame(
    RTSS = c("Seasonal", "Seasonal", "Hybrid", "none"),
    seasonal_ages = c("5-months-1-year", "5-months-3-years", "age-based-primary", "no-vaccine")
  )

# # replacing with 5 to 17 months to have a back up
# seasonal_ages <-
#   data.frame(
#     RTSS = c("Seasonal", "Hybrid", "none"),
#     seasonal_ages = c("5-months-17-months", "age-based-primary", "no-vaccine")
#   )

# join together
combo <-
  combo |>
  left_join(seasonal_ages)

# make a new run name long name to save by
combo <-
  combo |>
  mutate(
    run_name = paste(pfpr, seas_name, drawID, RTSS,
                     dosing_assumption, coverage_assumption, seasonal_ages,
                     "smc", SMC, vaccine_model, "with-age-cap", sep = "-")
  )

# put variables into the same order as function arguments
combo <-
  combo |>
  select(population,        # simulation population
         seasonality,       # seasonal profile
         seas_name,         # name of seasonal profile
         pfpr,              # corresponding PfPR
         warmup,            # warm-up period
         sim_length,        # length of simulation run
         speciesprop,       # proportion of each vector species
         treatment,         # treatment coverage
         SMC,               # SMC coverage
         RTSS,              # RTS,S strategy
         coverage_assumption, # RTS,S coverage assumption
         seasonal_ages,      # Seasonal age values
         dosing_assumption, # RTSS dose number
         ID,                # name of output file
         drawID,            # parameter draw no.
         vaccine_model,     # vaccine model
         run_name           # long ID name
  ) |> as.data.frame()

saveRDS(combo, "./01.data/vaccine-scenarios/runscenarios.rds")

#-Generate the parameter lists for each scenario-------------------------------------------
source("./02.scripts/functions/f1-generate-parameters.R")

generate_params("./01.data/vaccine-scenarios/runscenarios.rds",   # file path to pull
                "./01.data/vaccine-scenarios/run_parameters.rds") # file path to push

# glance
paramlist <- readRDS("./01.data/vaccine-scenarios/run_parameters.rds")

# Function: check_paramlist_params_column
# Purpose: Check that parameter list-column rows have the expected nested structure.
# Inputs: paramlist data frame and verbose flag.
# Outputs: Integer indices of malformed rows (empty if none).
# Assumptions: Assumes params column contains list objects from generate_params().

check_paramlist_params_column <- function(paramlist, verbose = TRUE) {
  bad_rows <- which(sapply(paramlist$params, function(p) {
    is.list(p) && length(p) == 1 && is.list(p[[1]]) && !is.null(names(p[[1]]))
  }))

  if (verbose) {
    if (length(bad_rows) == 0) {
      message("✅ All rows in `params` column are correctly structured.")
    } else {
      message("❌ Found ", length(bad_rows), " improperly nested `params` rows.")
      message("Affected row indices: ", paste(bad_rows, collapse = ", "))

      # Show the structure of the first bad example
      cat("\n🔍 Structure of first problematic row:\n")
      print(str(paramlist$params[[bad_rows[1]]], max.level = 2))
    }
  }

  return(bad_rows)
}

bad_rows <- check_paramlist_params_column(paramlist)

#-Run sims--------------------------------------------------------------------------------
x <- c(1:nrow(readRDS("./01.data/vaccine-scenarios/runscenarios.rds"))) # number of runs

# define all combinations of scenarios and draws - only run draw 0
# index <- tibble(x = x)
# Filter for drawID == 0 and get the row numbers
# This enforces the study uncertainty design (intervention efficacy uncertainty only).
rows_to_run <- which(paramlist$drawID == 0)

# Create index tibble
index <- tibble(x = rows_to_run)

source("./02.scripts/functions/f3-runsims.R")

# # Example run
# runsim(x=index$x[777])

# run all combos---------------------------------------------
library(furrr)

# Set up parallel processing (default 40; override with S3_WORKERS_MAIN)
workers_s3_main <- suppressWarnings(as.integer(Sys.getenv("S3_WORKERS_MAIN", "40")))
if (is.na(workers_s3_main) || workers_s3_main < 1) workers_s3_main <- 40
plan(multisession, workers = workers_s3_main)

# Run PRmatch in parallel
results <- future_map(index$x, runsim, .progress = TRUE)

# Shut down parallel workers
plan(sequential)

#-find which ones didn't run-----------------------------------------------------------
# Paths
output_dir <- "./03.outputs/individual-model-sim-outputs/"
runscenarios_path <- "./01.data/vaccine-scenarios/runscenarios.rds"

# Expected run indices (already filtered for drawID == 0)
expected_runs <- index$x

# Get all saved filenames
saved_files <- list.files(output_dir, full.names = FALSE)

# Extract scenario indices from filenames (assuming filenames like "1234-runname.rds")
saved_indices <- saved_files |>
  str_extract("^\\d+") |>
  as.integer()

# Identify which expected indices are missing
missing_runs <- setdiff(expected_runs, saved_indices)

# Summary
cat("Total expected runs:", length(expected_runs), "\n")
cat("Completed runs:", length(saved_indices), "\n")
cat("Missing runs:", length(missing_runs), "\n")

#-run missing scenarios----------------------------------------------------------------
missed_index <- tibble(x = missing_runs)

# Set up parallel processing for missing scenarios (default 50; override with S3_WORKERS_MISSING)
workers_s3_missing <- suppressWarnings(as.integer(Sys.getenv("S3_WORKERS_MISSING", "50")))
if (is.na(workers_s3_missing) || workers_s3_missing < 1) workers_s3_missing <- 50
plan(multisession, workers = workers_s3_missing)

# Run PRmatch in parallel
results <- future_map(missed_index$x, runsim, .progress = TRUE)

# Shut down parallel workers
plan(sequential)

