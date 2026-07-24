# ------------------------------------------------------------------------------
# Purpose: Generate baseline scenarios and calibrate PfPR-to-EIR matching outputs.
# Inputs: Baseline draw inputs and scenario configuration constants.
# Outputs: Calibration artifacts in 03.outputs/pf-eir-match and baseline parameter files.
# Dependencies: data_and_libraries.R, f1-generate-parameters.R, f2-calibrate.R.
# Run Stage: Pipeline step 2: baseline calibration.
# ------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# 3-65% pfpr | 2 seasonality profiles | SMC or no SMC
#-------------------------------------------------------------------------------

#-libraries---------------------------------------------------------------------
source("./02.scripts/data_and_libraries.R")

#-set up------------------------------------------------------------------------
# make all combinations of baseline scenarios

# defined parameters universal
year <- 365            #year
population <- 100000   #population size

# run time
warmup <- 6 * year       # multiple of 3
sim_length <- 3 * year   # value > 0

# number of parameter draws
drawID <- readRDS("./01.data/inputs/uncertainty-params/parameter_draws.rds")$draw # sampled in s1_parameter_draws.R
# drawID_R21 <- c(1:50)
draws <- tibble(drawID)
draws <- rbind(draws, data.frame(drawID = 0))

#-transmission set-up-----------------------------------------------------------
# parasite prevalence 2-10 year olds
pfpr <- c(0.01, 0.03, 0.05, 0.1, 0.15, 0.2,
          0.25, 0.35, 0.45, 0.55, 0.65)

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
RTSS <- c("none")

# RTS,S coverage
RTSScov <- c(0)

# adding a fifth, sixth, seventh RTS,S dose: 0, 1
RTSS_4 <- c(0)
RTSS_5 <- c(0)
RTSS_6 <- c(0)
RTSS_7 <- c(0)

interventions <- crossing(treatment, SMC, RTSS, RTSScov, RTSS_4,
                          RTSS_5, RTSS_6, RTSS_7)

vaccine_model <- "baseline"
#-create combination of all runs------------------------------------------------
combo <-
  crossing(population, pfpr, stable, warmup, sim_length,
           speciesprop, interventions, draws, vaccine_model) |>
  mutate(ID = paste(pfpr, seas_name, drawID, SMC, sep = "_"))

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
         RTSScov,           # RTS,S coverage
         RTSS_4,            # status of 4th dose for SV or hybrid strategies
         RTSS_5,            # status of 5th dose for SV or hybrid strategies
         RTSS_6,            # status of 6th dose for SV or hybrid strategies
         RTSS_7,            # status of 7th dose for SV or hybrid strategies
         ID,                # name of output file
         drawID,            # parameter draw no.
         vaccine_model      # vaccine model
  ) |> as.data.frame()

saveRDS(combo, "./01.data/baseline/baselinescenarios.rds")

#-generate parameter list-------------------------------------------------------
#in malariasimulation format and output as baseline_parameters.rds
source("./02.scripts/functions/f1-generate-parameters.R")

generate_params("./01.data/baseline/baselinescenarios.rds",    # file path to pull
                "./01.data/baseline/baseline_parameters.rds")  # file path to push

# glance
paramlist <- readRDS("./01.data/baseline/baseline_parameters.rds")

#-Run cali-- -------------------------------------------------------------------
x <- c(1:nrow(combo)) # baseline scenarios

# define all combinations of scenarios and draws
index <- tibble(x = x, y = combo$drawID, ID = combo$ID)

source("./02.scripts/functions/f2-calibrate.R")

# # Example run
# PRmatch(x=index$x[1452], y = index$y[1452])

# run all combos---------------------------------------------
library(furrr)

# Set up parallel processing (default 44; override with S2_WORKERS)
workers_s2 <- suppressWarnings(as.integer(Sys.getenv("S2_WORKERS", "44")))
if (is.na(workers_s2) || workers_s2 < 1) workers_s2 <- 44
plan(multisession, workers = workers_s2)

index0 <- index |> filter(y == 0)

# Run PRmatch in parallel
results <- future_map2(index0$x, index0$y, PRmatch, .progress = TRUE)

# Shut down parallel workers
plan(sequential)

#-Results-----------------------------------------------------------------------
# read in results
files <-
  list.files("./03.outputs/pf-eir-match",
             pattern = "PRmatch_draws_",
             full.names = TRUE)
dat_list <- lapply(files, function (x) readRDS(x))

# concatenate
match <-  do.call("rbind", dat_list) |> as_tibble()

summary(match$starting_EIR)

# take a look at failed jobs
anti_join(combo, match, by = "ID") # Nothing failed

# save EIR estimates
saveRDS(match, "./03.outputs/pf-eir-match/EIRestimates.rds")
