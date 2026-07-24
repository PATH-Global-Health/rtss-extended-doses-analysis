# ------------------------------------------------------------------------------
# Purpose: Create uncertainty draw inputs used by downstream scenario generation.
# Inputs: Intervention uncertainty CSV files in 01.data/inputs/uncertainty-params.
# Outputs: Saved draw subsets for RTS,S and SMC uncertainty inputs.
# Dependencies: data_and_libraries.R packages and uncertainty input files.
# Run Stage: Pipeline step 1: uncertainty input preparation.
#-------------------------------------------------------------------------------

# libraries
source("./02.scripts/data_and_libraries.R")

#-malaria simulation draws------------------------------------------------------

# choose a random sample of 50 draws to use
# in malariasimulation::set_parameter_draw()
set.seed(123)
rsample <- sample(x = c(1:1000), size = 50, replace = FALSE, prob = NULL)
rsample <- tibble(draw = rsample)

saveRDS(rsample, "./01.data/inputs/uncertainty-params/parameter_draws.rds")

#-RTS,S parameter draw----------------------------------------------------------
rtss_draw <-
  read_csv("./01.data/inputs/uncertainty-params/reduced_vacc_draws.csv") |>
  bind_cols(rsample) |>
  select(-param_draw)

saveRDS(rtss_draw, "./01.data/inputs/uncertainty-params/rtss_draws.rds")


#-SMC parameter draw------------------------------------------------------------
smc_draw <-
  read_csv("./01.data/inputs/uncertainty-params/reduced_spaq_draws.csv") |>
  bind_cols(rsample) |>
  select(-param_draw)

saveRDS(smc_draw, "./01.data/inputs/uncertainty-params/smc_draws.rds")
