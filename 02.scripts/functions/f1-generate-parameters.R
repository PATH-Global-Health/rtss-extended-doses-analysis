# ------------------------------------------------------------------------------
# Purpose: Build simulation parameter payloads from scenario inputs.
# Inputs: Scenario table RDS with intervention, transmission, and draw settings.
# Outputs: Parameter table RDS used by downstream calibration and simulation steps.
# Dependencies: data_and_libraries.R packages; scenario input files in 01.data.
# Run Stage: Function library used by baseline and model-run scripts.
# Notes: Comment-only documentation update; logic and outputs unchanged.
# ------------------------------------------------------------------------------

# RTS,S coverage assumptions
# i.	Scenario 1 = 90% coverage no drop out
# ii.	 Scenario 2 = 90% coverage d3, d3-d4 drop out 10% and 10% drop out total over the remaining doses
# iii.	Scenario 3 = 90% coverage d3, 5% drop out each additional dose
# iiii.  Scenario 4 = 90% coverage d3, 10% drop out each additional dose

# Potential for more granular analysis on drop out rates as required
coverage_assumption <- c(
  "S1-optimal-campaigns",
  "S2-enhanced-communication",
  "S3-sustained-dropout-5perc",
  "S3-sustained-dropout-10perc"
)


#-helper function to set rtss coverage----------------------------------------------------------------------------------------------------------------------
# Function: get_rtss_retention
# Purpose: Map coverage assumption to per-dose RTS,S retention multipliers.
# Inputs: assumption (character), dose_number (numeric).
# Outputs: Numeric retention fraction for the requested dose.
# Assumptions: Assumes project-standard coverage assumption labels.

get_rtss_retention <- function(assumption, dose_number) {
  if (assumption == "S1-optimal-campaigns") return(1)
  if (assumption == "S3-sustained-dropout-5perc") return(0.95)
  if (assumption == "S3-sustained-dropout-10perc") return(0.9)

  if (assumption == "S2-enhanced-communication") {
    if (dose_number == "3-doses") return(numeric(0))
    if (dose_number == "4-doses") return(c(0.9))
    if (dose_number == "5-doses") return(c(0.9, 0.9))
    if (dose_number == "6-doses") return(c(0.9, 0.95, 0.95))
    if (dose_number == "7-doses") return(c(0.9, 0.97, 0.97, 0.96))
  }

  stop(paste("Unknown assumption or dose number:", assumption, dose_number))
}

#-Generate parameters for malariasimulation runs translate input dataframe to parameter list-----------------------------------------------------------------
# Function: generate_params
# Purpose: Generate scenario-specific parameter lists for simulation runs.
# Inputs: inputpath (RDS scenarios), outputpath (target RDS).
# Outputs: Writes parameter table with list-column of model parameters.
# Assumptions: Assumes scenario schema used by this repository.

generate_params <- function(inputpath,   # path to input scenarios
                            outputpath){ # path where output file will be stored

  # read in dataframe of all scenario combinations
  scenarios <- readRDS(inputpath)

  #-generate parameters------------------------------------------------------------
  generate_params2 <- function(x){ # x = scenario number

    #-pull one scenario at a time-----------------------
    data <- scenarios[x, ]

    #-assign values from inputs-------------------------
    population = data$population
    seasonality = data$seasonality
    seas_name = data$seas_name
    pfpr = data$pfpr
    warmup = data$warmup
    sim_length = data$sim_length
    speciesprop = data$speciesprop
    treatment = data$treatment
    SMC = data$SMC
    RTSS = data$RTSS
    dosing_assumption <- data$dosing_assumption
    coverage_assumption <- data$coverage_assumption
    ID = data$ID
    drawID = data$drawID
    vaccine_model = data$vaccine_model
    seasonal_ages = data$seasonal_ages
    run_name = data$run_name

    #-timesteps values---------------------------------
    year <- 365
    month <- year / 12

    #-starting parameters-----------------------------
    params <-
      get_parameters(
        list(
          human_population = population,
          model_seasonality = TRUE,
          # rainfall fourier parameters
          g0 = unlist(seasonality)[1],
          g = unlist(seasonality)[2:4],
          h = unlist(seasonality)[5:7],
          individual_mosquitoes = FALSE
        )
      )

    #-parameter draws-------------------------------
    if(drawID >0){
      params <- set_parameter_draw(
        parameters = params,
        draw = drawID)
    }

    # RTSS vax and booster profiles-----------------
    RTSS_profile <- malariasimulation::rtss_profile
    RTSS_booster_profile <- malariasimulation::rtss_booster_profile


    if(drawID > 0 ){
      # pull values from RTSS parameter draws
      RTSS_params <-
        readRDS("./01.data/inputs/uncertainty-params/rtss_draws.rds")|>
        filter(draw == drawID)

      # adjust for parameter draw
      RTSS_profile <- create_pev_profile(
        vmax = RTSS_params$pev_Vmax,
        alpha = RTSS_params$pev_alpha,
        beta = RTSS_params$pev_beta,
        cs = RTSS_profile$cs,
        rho = RTSS_profile$rho,
        ds = RTSS_profile$ds,
        dl = RTSS_profile$dl
      )

      # adjsut for parameter draw
      RTSS_booster_profile <-
        create_pev_profile(
          vmax = RTSS_params$pev_Vmax,
          alpha = RTSS_params$pev_alpha,
          beta = RTSS_params$pev_beta,
          cs = RTSS_booster_profile$cs,
          rho = RTSS_booster_profile$rho,
          ds = RTSS_booster_profile$ds,
          dl = RTSS_booster_profile$dl
        )
    }

    # modify for the different vaccine models if improved efficacy needed
    if(vaccine_model %in% c("model-2", "model-3")){
      RTSS_booster_profile$cs <- c(6.37008, 0.35)
    }


    #-SMC parameter draws------------------------
    if(drawID >0){
      smc_draws <-
        readRDS("./01.data/inputs/uncertainty-params/smc_draws.rds") |>
        filter(draw == drawID)
    }

    if(drawID == 0){
      smc_draws <-
        data.frame(scale = 39.34,
                   shape = 3.40,
                   draw = 0)
    }

    #-modify in model_3--------------------------

    # if including 'synergy'
    # synergy ----------
    if(vaccine_model == "model-3"){

      # Thompson et al. 2022
      # vaccine
      RTSS_profile$beta  <- 70.9
      RTSS_profile$alpha <- 0.868
      RTSS_profile$vmax  <- 0.843
      RTSS_booster_profile$beta  <- 70.9
      RTSS_booster_profile$alpha <- 0.868
      RTSS_booster_profile$vmax  <- 0.843

      # smc
      smc_draws$scale <- c(45.76)
      smc_draws$shape <- c(2.87)

    }

    #-outcome definitions-------------------------
    # incidence for every 1 year age group up to 20, 10-year age groups from 20 to 100 years and a 0-5 category

    # 0–5 group (standalone)
    special_group_min <- 0
    special_group_max <- 5

    # 1-year bands from 0–20
    min_1yr <- seq(0, 19, 1)
    max_1yr <- seq(1, 20, 1)

    # 10-year bands from 20–100
    min_10yr <- seq(20, 90, 10)
    max_10yr <- seq(30, 100, 10)

    # Combine all
    min_ages <- c(special_group_min, min_1yr, min_10yr)
    max_ages <- c(special_group_max, max_1yr, max_10yr)

    # Convert to days
    params$clinical_incidence_rendering_min_ages <- min_ages * year
    params$clinical_incidence_rendering_max_ages <- max_ages * year
    params$severe_incidence_rendering_min_ages   <- min_ages * year
    params$severe_incidence_rendering_max_ages   <- max_ages * year

    # prevalence 2-10 year olds
    params$prevalence_rendering_min_ages = 2 * year
    params$prevalence_rendering_max_ages = 10 * year

    #-demography---------------------------------
    africa_demog <- read.csv("./01.data/inputs/demography/ssa_demography_2021.csv")
    ages <- round(africa_demog$age_upper * year) # top of age bracket
    deathrates <- africa_demog$mortality_rate / 365   # age-specific death rates

    params <- set_demography(
      params,
      agegroups = ages,
      timesteps = 0,
      deathrates = matrix(deathrates, nrow = 1))

    #-vectors---------------------------------
    params <- set_species(
      parameters = params,
      species = list(arab_params, fun_params, gamb_params),
      proportions = unlist(speciesprop))

    # proportion of bites taken in bed for each species
    # find values in S.I. of 10.1038/s41467-018-07357-w Table 3
    params$phi_bednets <- c(0.9, 0.9, 0.89) # Hogan et al. 2020
    # proportion of bites taken indoors for each species
    params$phi_indoors <- c(0.96, 0.98, 0.97) # Hogan et al. 2020

    #-treatment------------------------------
    if (treatment > 0) {
      params <- set_drugs(
        parameters = params,
        list(AL_params, SP_AQ_params))

      # AL default, SP: https://doi.org/10.1016/S2214-109X(22)00416-8 supplement
      params$drug_prophylaxis_scale <- c(10.6, smc_draws$scale)
      params$drug_prophylaxis_shape <- c(11.3, smc_draws$shape)

      params <- set_clinical_treatment(
        parameters = params,
        drug = 1,
        timesteps = c(1),
        coverages = c(treatment)
      )  }

    #-SMC-----------------------------------
    smc_timesteps <- 0

    if (SMC > 0) {
      peak <- peak_season_offset(params)

      if(seas_name == "highly seasonal"){
        # 4 doses, centered around the peak
        first <- round(c(peak + c(-1, 0, 1, 2) * month), 0)
        firststeps <- sort(rep(first, (warmup + sim_length)/year))
        yearsteps <- rep(c(0, seq(year, (warmup + sim_length) - year, year)), length(first))
        timesteps <- yearsteps + firststeps

      }

      if(seas_name == "seasonal"){
        # 5 doses, centered around peak
        first <- round(c(peak + c(-2, -1, 0, 1, 2) * month), 0)
        firststeps <- sort(rep(first, (warmup + sim_length)/year))
        yearsteps <- rep(c(0, seq(year, (warmup + sim_length) - year, year)), length(first))
        timesteps <- yearsteps + firststeps

      }

      params <- set_drugs(
        parameters = params,
        list(AL_params, SP_AQ_params)
      )

      # AL default, SP: https://doi.org/10.1016/S2214-109X(22)00416-8 supplement
      params$drug_prophylaxis_scale <- c(10.6, smc_draws$scale)
      params$drug_prophylaxis_shape <- c(11.3, smc_draws$shape)

      params <- set_smc(
        parameters = params,
        drug = 2,
        timesteps = sort(timesteps),
        coverages = rep(SMC, length(timesteps)),
        min_ages = rep(round(0.25 * year),length(timesteps)),
        max_ages = rep(round(5 * year),length(timesteps))
      )

      # var for outputting to check SMC timings are correct
      smc_timesteps <- params$smc_timesteps - warmup
    }


    # #-RTS,S EPI ---------------------------------------------------------
    # # Age based 6,7,9 and 15 months later for a booster if RTSS_4 == 1
    # # otherwise just a three dose schedule
    # if (RTSS == "Age-based") {
    #   params$pev_doses <- round(c(0, 1 * month, 3 * month))
    #
    #   # boosters variable - need values but where 3 dose set coverage to 0
    #   boosters <-  round(c(15 * month))
    #   booster_cov <- if(RTSS_4 == 0) matrix(0) else matrix(.80)
    #   booster_prof <- rep(list(RTSS_booster_profile), length(boosters))
    #
    #   params <- set_pev_epi(
    #     parameters = params,
    #     profile = RTSS_profile,
    #     coverages = RTSScov,
    #     timesteps = warmup,
    #     age = round(6 * month), # 6, 7, 9 months
    #     min_wait = 0,
    #     booster_spacing = boosters,
    #     booster_coverage = booster_cov,
    #     booster_profile = booster_prof,
    #     seasonal_boosters = FALSE
    #   )
    #
    #   }


    #-RTS,S SV -------------------------------------------------------------------------------
    rtss_mass_timesteps <- 0

    if (RTSS == "Seasonal") {
      # dose values
      n_doses_total <- as.integer(stringr::str_extract(dosing_assumption, "\\d+"))
      n_boosters <- n_doses_total - 3

      #offsets
      peak <- peak_season_offset(params)
      if(seas_name == "highly seasonal"){
        first <- round(warmup + (peak - month * 3.5), 0)}
      if(seas_name == "seasonal"){
        first <- round(warmup + (peak - month * 5.5), 0)}
      timesteps <- c(first, first+seq(year, sim_length, year))
      params$pev_doses <- round(c(0, 1 * month, 2 * month))

      # ages
      age_lower <- 5 * month
     if(seasonal_ages == "5-months-1-year"){
       age_upper = 1*year
     }
     if(seasonal_ages == "5-months-3-years"){
       age_upper = 3*year
     }
     if(seasonal_ages == "5-months-17-months"){
       age_upper = 17*month
     }

      # variable boosters and coverage and profile setting
      primary_cov <- 0.9#ifelse(coverage_assumption == "A1-100pct-all", 1, 0.9)       # 100% if set otherwise 90%
      boosters <- round(seq(12 * month, by = 12 * month, length.out = n_boosters)) # booster time steps annually dependent on dose number specified
      retention <- get_rtss_retention(coverage_assumption, dosing_assumption)      # get retention assumption ()

        if (coverage_assumption == "S2-enhanced-communication"){
          booster_cov <- matrix(rep(retention, times = length(timesteps)),
                 nrow = length(timesteps),
                 byrow = TRUE)
        } else {
          booster_cov <- matrix(retention, nrow=length(timesteps), ncol=length(boosters))
        }

        # set up for malaria sim format
      booster_prof <- rep(list(RTSS_booster_profile), length(boosters))            # efficacy profile

      params <- set_mass_pev(
        parameters = params,
        profile = RTSS_profile,
        timesteps = timesteps,
        coverages = rep(primary_cov, length(timesteps)),
        min_ages = round(age_lower),
        max_ages = round(age_upper),
        min_wait = 15*year,
        booster_spacing = boosters,
        booster_coverage = booster_cov,
        booster_profile = booster_prof,
        vaccine_max_age_cap = 1825 # Children stop getting boosters at 5 years old
      )

      # var for outputting to check RTS,S timings are correct
      rtss_mass_timesteps <- params$rtss_mass_timesteps - warmup

    }

    #-RTS,S hybrid -------------------------------------------------------------------------
    if (RTSS == "Hybrid") {

      # set dose spacing
      params$pev_doses <- round(c(0, 1 * month, 3 * month))

      # set timestep of first seasonal booster
      peak <- peak_season_offset(params)

      if(seas_name == "highly seasonal"){
        boost <- round((peak - month * 1.5), 0)}
      if(seas_name == "seasonal"){
        boost <- round((peak - month * 3.5), 0)}

      # variable boosters and coverage
      # Get total dose count and number of boosters
      n_doses_total <- as.integer(stringr::str_extract(dosing_assumption, "\\d+"))
      n_boosters <- n_doses_total - 3

      # Booster timing: annual from boost offset
      boosters <- boost + seq(0, by = year, length.out = n_boosters)

      # Booster coverage based on assumption
      retention <- get_rtss_retention(coverage_assumption, dosing_assumption)

      booster_cov <- matrix(retention, nrow=1, ncol=length(boosters)) # set up for malaria sim format
      booster_prof <- rep(list(RTSS_booster_profile), length(boosters))            # efficacy profile

      # primary coverage
      primary_cov <- 0.9#ifelse(coverage_assumption == "A1-100pct-all", 1, 0.9)       # 100% if set otherwise 90%

      # set parameters
      params <- set_pev_epi(
        parameters = params,
        profile = RTSS_profile,
        coverages = primary_cov,
        timesteps = warmup,
        age = round(5 * month),  # 5, 6, 7 months
        min_wait = round(6 * month), #6 month min wait between getting seasonal doses
        booster_spacing = boosters,
        booster_coverage = booster_cov,
        booster_profile = booster_prof,
        seasonal_boosters = TRUE,
        vaccine_max_age_cap = 1825 # Children stop getting boosters at 5 years old
        )

    }

    #-save as data.frame------------------------------------------------------------------------
    data$params <- list(params)
    data$scenarioID <- x

    # print count
    print(x)

    return(data)

  }

  #-loop through function to generate parameters one by one-------------------------------------
  output <- map_dfr(1:nrow(scenarios), generate_params2)


  #-save output----------------------------------------------------------------------------------
  saveRDS(output, outputpath)

}
