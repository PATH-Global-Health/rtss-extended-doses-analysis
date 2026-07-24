# ------------------------------------------------------------------------------
# Purpose: Generate efficacy-profile figures for SP+AQ and seasonal RTS,S model assumptions.
# Inputs: Uncertainty input model objects and efficacy parameter draws.
# Outputs: Supporting efficacy profile figures for manuscript interpretation.
# Dependencies: Specialized plotting/model packages loaded in this script.
# Run Stage: Supporting analysis: intervention efficacy profiles.
# Notes: Requires model.4.RDS, a large fitted model object not stored in GitHub.
# ------------------------------------------------------------------------------

#-Libraries---------------------------------------------------------------------
library(compiler)
library(leaflet)
library(sf)
library(rnaturalearth)
library(webshot2)
library(htmlwidgets)
library(tidyverse)
library(lubridate)
# library(PATHtools)
library(umbrella)
library(terra)
library(openxlsx)
library(purrr)
library(patchwork)

#-SP+AQ BASELINE----------------------------------------------------------------

## load data----------------------------
model_4_path <- "./01.data/inputs/uncertainty-params/model.4.RDS"

if (!file.exists(model_4_path)) {
  stop(
    paste(
      "Missing required fitted model object:",
      model_4_path,
      "This large file is not included in the GitHub repository.",
      "It is only needed to regenerate intervention-efficacy.png.",
      "Place model.4.RDS at this path before running this supporting figure script."
    ),
    call. = FALSE
  )
}

model.4 <-
  readRDS(model_4_path)

post_draws_model <- rstan::extract(model.4)
output_mcmc <- cbind(post_draws_model[["lambda"]], post_draws_model[["k"]])
par_med <- c(median(post_draws_model[["lambda"]]), median(post_draws_model[["k"]]))
time_seq <- seq(0, 60, 1)

## protection over time function-------
# Function: model_PE
# Purpose: Compute protective efficacy at a time point from fitted parameters.
# Inputs: time (numeric), pars (lambda/k parameter vector).
# Outputs: Numeric efficacy estimate for each time value.
# Assumptions: Assumes parameterization from fitted efficacy model.

model_PE <- function(time, pars) {
  lambda <- pars[1]
  k <- pars[2]
  ## Drug response - protection against infection
  drug_eff <- exp(-(time / lambda)^k)
}

## sample and create output------------
# sampled across posterior draws to calculate CrIs
model_PE <- cmpfun(model_PE, options = list(optimize = 3))
median_pred <- sapply(time_seq, model_PE, par = par_med)

# 20,000 samples
N_sam <- 20000

# indexes
sam_seq <- round(seq(from = 1, to = nrow(output_mcmc), length = N_sam))

# empty matrix to store values in
sam_cs <- matrix(NA, nrow = N_sam, ncol = length(time_seq))

for (k in 1:N_sam) {
  sam_cs[k, ] <- sapply(time_seq, model_PE, par = output_mcmc[sam_seq[k], 1:2])
}

quant_cs <- matrix(NA, nrow = 3, ncol = length(time_seq))

for (j in 1:length(time_seq)) {
  quant_cs[, j] <- quantile(sam_cs[, j], prob = c(0.025, 0.5, 0.975))
}

quant_cs5 <- matrix(NA, nrow = 3, ncol = length(time_seq))

for (j in 1:length(time_seq)) {
  quant_cs5[, j] <- quantile(sam_cs[, j], prob = c(0.25, 0.5, 0.75))
}

quant_cs9 <- matrix(NA, nrow = 3, ncol = length(time_seq))

for (j in 1:length(time_seq)) {
  quant_cs9[, j] <- quantile(sam_cs[, j], prob = c(0.001, 0.5, 1))
}

## output-----------------------------
smc_bl_efficacy <- data.frame(
  time = time_seq,
  median = quant_cs[2, ],
  lower_50 = quant_cs5[1, ],
  upper_50 = quant_cs5[3, ],
  lower = quant_cs[1, ],
  upper = quant_cs[3, ],
  name = "SP+AQ",
  model = "Model 1"
)

# replicate for model 2 as no changes here
smc_m2_efficacy <- data.frame(
  time = time_seq,
  median = quant_cs[2, ],
  lower_50 = quant_cs5[1, ],
  upper_50 = quant_cs5[3, ],
  lower = quant_cs[1, ],
  upper = quant_cs[3, ],
  name = "SP+AQ",
  model = "Model 2"
)


#-BASELINE RTSS SEASONAL VACCINE PROFILE----------------------------------------

## Functions Anitbodies-----------------------------

# draw from logitnormal distribution.
#  Input arguments include the mean and standard deviation
#  of the RAW normal variate that will be converted
#  to the logit scale. Therefore the mean and
#  standard deviation of the final random variable will
#   not equal these values.
# Function: rlogitnorm
# Purpose: Sample from a logit-normal distribution for bounded efficacy inputs.
# Inputs: n, mean_raw, sd_raw.
# Outputs: Numeric vector on [0,1] scale.
# Assumptions: Assumes raw parameters are on logit scale.

rlogitnorm <- function(n, mean_raw, sd_raw) {
  x <- exp(rnorm(n, mean_raw, sd_raw))
  x / (1 + x)
}

# draw from lognormal distribution.
# Input arguments include the desired mean and standard
#  deviation of the final lognormal random variable;
#  the mean and standard deviation of the raw normal
#   variate are calculated from these values.
# Function: rlnorm2
# Purpose: Generate log-normal draws using mean and sd on natural scale.
# Inputs: n, mean, sd.
# Outputs: Numeric vector of log-normal samples.
# Assumptions: Assumes positive-valued target variable.

rlnorm2 <- function(n, mean, sd) {
  meanlog <- log(mean^2 / sqrt(mean^2 + sd^2)) # mean of the lognormal distribution
  sdlog <- sqrt(log(1 + (sd / mean)^2)) # standard deviation of the lognormal distribution
  rlnorm(n, meanlog, sdlog) # draw from the lognormal distribution
}

# draw single antibody titre from the model
# Function: antibody_titre_sv
# Purpose: Simulate seasonal-vaccination antibody titre trajectories.
# Inputs: Dose-level distribution parameters and timing settings.
# Outputs: Data frame/vector of simulated titres by time.
# Assumptions: Assumes model-1 seasonal vaccination structure and timing.

antibody_titre_sv <- function(d1_mu = 45, d1_sigma = 16, d2_mu = 591, d2_sigma = 245,
                              rho_mu = 2.37832, rho_sigma = 1.00813, ab_mu = 621, ab_sigma = 0.35,
                              rho_mu_boost = 1.03431, rho_sigma_boost = 1.02735,
                              t_boost = 1, t_boost1 = 1, t_boost2 = 1, t_boost3 = 1,
                              ab_mu_boost = 277, ab_sigma_boost = 0.35, years = 6) {
  # if 7 doses
  if (!is.null(t_boost2) & !is.null(t_boost3)) {
    # draw parameters from distributions
    rho <- rlogitnorm(1, mean_raw = rho_mu, sd_raw = rho_sigma)
    rho_boost <- rlogitnorm(1, mean_raw = rho_mu_boost, sd_raw = rho_sigma_boost)
    d1 <- rlnorm2(1, mean = d1_mu, sd = d1_sigma) # half-life of short-lived component of antibody response
    d2 <- rlnorm2(1, mean = d2_mu, sd = d2_sigma) # half-life of long-lived component of antibody response
    ab0 <- exp(rnorm(1, log(ab_mu) - ab_sigma^2 / 2, sd = ab_sigma)) # ab0=CS_peak
    ab0_boost <- exp(rnorm(1, log(ab_mu_boost) - ab_sigma_boost^2 / 2, sd = ab_sigma_boost))
    ab0_boost1 <- exp(rnorm(1, log(ab_mu_boost) - ab_sigma_boost^2 / 2, sd = ab_sigma_boost))
    ab0_boost2 <- exp(rnorm(1, log(ab_mu_boost) - ab_sigma_boost^2 / 2, sd = ab_sigma_boost))
    ab0_boost3 <- exp(rnorm(1, log(ab_mu_boost) - ab_sigma_boost^2 / 2, sd = ab_sigma_boost))

    # Simulate antibody titre over time
    r1 <- log(2) / d1
    r2 <- log(2) / d2

    ab <- vector(mode = "numeric", length = 365 * years)

    tvec1 <- 1:t_boost
    ab[1:t_boost] <- ab0 * (rho * exp(-r1 * tvec1) + (1 - rho) * exp(-r2 * tvec1))

    tvec2 <- seq((t_boost + 1), (365 * years))
    ab[(t_boost + 1):(365 * years)] <- ab0_boost * (rho_boost * exp(-r1 * (tvec2 - t_boost)) + (1 - rho_boost) * exp(-r2 * (tvec2 - t_boost)))

    tvec3 <- seq((t_boost1 + 1), (365 * years))
    ab[(t_boost1 + 1):(365 * years)] <- ab0_boost1 * (rho_boost * exp(-r1 * (tvec3 - t_boost1)) + (1 - rho_boost) * exp(-r2 * (tvec3 - t_boost1)))

    tvec4 <- seq((t_boost2 + 1), (365 * years))
    ab[(t_boost2 + 1):(365 * years)] <- ab0_boost2 * (rho_boost * exp(-r1 * (tvec4 - t_boost2)) + (1 - rho_boost) * exp(-r2 * (tvec4 - t_boost2)))

    tvec5 <- seq((t_boost3 + 1), (365 * years))
    ab[(t_boost3 + 1):(365 * years)] <- ab0_boost3 * (rho_boost * exp(-r1 * (tvec5 - t_boost3)) + (1 - rho_boost) * exp(-r2 * (tvec5 - t_boost3)))

    ret <- list(ab = ab)
    return(ret)
  } else {
    # draw parameters from distributions
    rho <- rlogitnorm(1, mean_raw = rho_mu, sd_raw = rho_sigma)
    rho_boost <- rlogitnorm(1, mean_raw = rho_mu_boost, sd_raw = rho_sigma_boost)
    d1 <- rlnorm2(1, mean = d1_mu, sd = d1_sigma) # half-life of short-lived component of antibody response
    d2 <- rlnorm2(1, mean = d2_mu, sd = d2_sigma) # half-life of long-lived component of antibody response
    ab0 <- exp(rnorm(1, log(ab_mu) - ab_sigma^2 / 2, sd = ab_sigma)) # ab0=CS_peak
    ab0_boost <- exp(rnorm(1, log(ab_mu_boost) - ab_sigma_boost^2 / 2, sd = ab_sigma_boost))
    ab0_boost1 <- exp(rnorm(1, log(ab_mu_boost) - ab_sigma_boost^2 / 2, sd = ab_sigma_boost))


    # Simulate antibody titre over time
    r1 <- log(2) / d1
    r2 <- log(2) / d2

    ab <- vector(mode = "numeric", length = 365 * years)

    tvec1 <- 1:t_boost
    ab[1:t_boost] <- ab0 * (rho * exp(-r1 * tvec1) + (1 - rho) * exp(-r2 * tvec1))

    tvec2 <- seq((t_boost + 1), (365 * years))
    ab[(t_boost + 1):(365 * years)] <- ab0_boost * (rho_boost * exp(-r1 * (tvec2 - t_boost)) + (1 - rho_boost) * exp(-r2 * (tvec2 - t_boost)))

    tvec3 <- seq((t_boost1 + 1), (365 * years))
    ab[(t_boost1 + 1):(365 * years)] <- ab0_boost1 * (rho_boost * exp(-r1 * (tvec3 - t_boost1)) + (1 - rho_boost) * exp(-r2 * (tvec3 - t_boost1)))


    ret <- list(ab = ab)
    return(ret)
  }
}

## Functions efficacy---------------------------------

# Simulate efficacy profile over time from the antibody titre
# Function: efficacy_profile
# Purpose: Translate titre values into efficacy trajectory over time.
# Inputs: alpha, beta, Vmax, titre, draw.
# Outputs: Efficacy profile outputs for plotting/summary.
# Assumptions: Assumes parameter draws follow fitted intervention model.

efficacy_profile <- function(alpha = 0.74, beta = 99.2, Vmax = 0.93, titre = 1, draw) {
  # Convert antibody titre ab to vaccine efficacy (Hill function)
  efficacy <- Vmax * (1 - 1 / (1 + (titre / beta)^alpha))

  return(efficacy)
}

## Initial stuff------------------------
years <- 6
tvec <- seq(0, 365 * years, length.out = 365 * years) # time vector
reps <- 5000 # number of runs

## Parameter values-----------------------------
# Anitbody model
# Fitted parameter values: from White et al (2015) Lancet ID
ab_mu <- 621 # median of the geometric means of observed antibody titres
ab_sigma <- 0.35 # observational variance of antibody titre (log-normal)
ab_sigma_boost <- 0.35
d1_mu <- 45
d1_sigma <- 16
d2_mu <- 591
d2_sigma <- 245
rho_mu <- 2.37832
rho_sigma <- 1.00813
rho_mu_boost <- 1.03431
rho_sigma_boost <- 1.02735

# Phase 3 profile Hill Function
alpha <- 0.74 # Shape parameter of dose-response curve
beta <- 99.2 # Scale parameter of dose-response curve
Vmax <- 0.93 # maximum efficacy against infection

## 5 doses total------------------------------
# Initialise matrices to store simulations of antibody
#  titre and vaccine efficacy over time.
ab_matrix_1 <- matrix(NA_real_, nrow = length(tvec), ncol = reps)
c1 <- matrix(NA_real_, nrow = length(tvec), ncol = 2)

### Simulate antibodies------------------------------
# seasonal vaccine profile 5 doses
t_boost <- 365 # time of booster dose in days after third vaccine dose
t_boost1 <- 365 * 2 # time of the second booster dose in days after the third vaccine dose
ab_mu_boost <- 277 # model 1 value

for (i in 1:reps) {
  sim_1 <- antibody_titre_sv(d1_mu, d1_sigma, d2_mu, d2_sigma,
    rho_mu, rho_sigma, ab_mu, ab_sigma,
    rho_mu_boost, rho_sigma_boost,
    t_boost, t_boost1,
    t_boost2 = NULL, t_boost3 = NULL,
    ab_mu_boost, ab_sigma_boost,
    years
  ) # run simulation
  ab_matrix_1[, i] <- sim_1$ab # store values
}

medianAB_1 <- apply(ab_matrix_1, 1, median)
quant_1 <- apply(ab_matrix_1, 1, quantile)

for (j in 1:length(tvec)) {
  c1[j, ] <- quantile(ab_matrix_1[j, ], c(0.025, 0.975))
}

write(medianAB_1, ncolumns = 1, file = "./02.scripts/supporting-figure-creation/titre_1.txt")
write.table(ab_matrix_1[, 1:100], file = "./02.scripts/supporting-figure-creation/ab_runs_1.txt", col.names = FALSE, row.names = FALSE)

### Simulate efficacy---------------------------------
# simulate the vaccine efficacy component of the model

# Read in outputs from antibody titre simulation
titre_1 <- read.table("./02.scripts/supporting-figure-creation/titre_1.txt", header = F)
ab_runs_1 <- read.table("./02.scripts/supporting-figure-creation/ab_runs_1.txt", header = F)

# Initialise matrices to store simulations of vaccine efficacy over time
# Initialise matrices to store simulations of vaccine efficacy over time
efficacy_matrix_1 <- matrix(NA_real_, nrow = length(tvec), ncol = 1)
efficacy_matrix_runs_1 <- matrix(NA_real_, nrow = length(tvec), ncol = 100)
c1 <- matrix(NA_real_, nrow = length(tvec), ncol = 2)

# 5 dose sRTS,S efficacy
sim1 <- efficacy_profile(alpha, beta, Vmax, titre_1[, 1, drop = TRUE])

for (i in 1:100) {
  efficacy_matrix_runs_1[, i] <- efficacy_profile(
    alpha, beta, Vmax,
    ab_runs_1[, i, drop = TRUE]
  )
}

c1 <- t(apply(
  efficacy_matrix_runs_1,
  1,
  stats::quantile,
  probs = c(0.025, 0.975),
  na.rm = TRUE
))

quant_1 <- apply(
  efficacy_matrix_runs_1,
  1,
  stats::quantile,
  na.rm = TRUE
)

srtss_efficacy_5dose_bl <- data.frame(
  time = tvec,
  median = sim1, # <- just the vector
  lower_50 = quant_1[2, ], # 25th percentile  (from apply(..., quantile))
  upper_50 = quant_1[4, ], # 75th percentile
  lower = c1[, 1], # 2.5th percentile (95% CI lower)
  upper = c1[, 2], # 97.5th percentile (95% CI upper)
  name = "RTSS 5 dose",
  model = "Model 1"
)

## 7 doses total-----------------------------
# Initialise matrices to store simulations of antibody
#  titre and vaccine efficacy over time.
ab_matrix_2 <- matrix(NA_real_, nrow = length(tvec), ncol = reps)
c2 <- matrix(NA_real_, nrow = length(tvec), ncol = 2)

### Simulate antibodies------------------------------

# seasonal vaccine profile 7 doses
t_boost <- 365 # time of booster dose in days after third vaccine dose
t_boost1 <- 365 * 2 # time of the second booster dose in days after the third vaccine dose
t_boost2 <- 365 * 3 # time of the third booster dose in days after the third vaccine dose
t_boost3 <- 365 * 4 # time of the fourth booster dose in days after the third vaccine dose
ab_mu_boost <- 277
# ab_mu_boost <- 621

for (i in 1:reps) {
  sim_2 <- antibody_titre_sv(
    d1_mu, d1_sigma, d2_mu, d2_sigma,
    rho_mu, rho_sigma, ab_mu, ab_sigma,
    rho_mu_boost, rho_sigma_boost,
    t_boost, t_boost1, t_boost2, t_boost3,
    ab_mu_boost, ab_sigma_boost,
    years
  ) # run simulation
  ab_matrix_2[, i] <- sim_2$ab # store values
}

medianAB_2 <- apply(ab_matrix_2, 1, median)
quant_2 <- apply(ab_matrix_2, 1, quantile)

for (j in 1:length(tvec)) {
  # c1[j,]<-quantile(ab_matrix_1[j,],c(0.025,0.975))
  c2[j, ] <- quantile(ab_matrix_2[j, ], c(0.025, 0.975))
  # c3[j,]<-quantile(ab_matrix_3[j,],c(0.025,0.975))
}

write(medianAB_2, ncolumns = 1, file = "./02.scripts/supporting-figure-creation/titre_2.txt")
write.table(ab_matrix_2[, 1:100], file = "./02.scripts/supporting-figure-creation/ab_runs_2.txt", col.names = FALSE, row.names = FALSE)

### Simulate efficacy---------------------------------
# simulate the vaccine efficacy component of the model
# Read in outputs from antibody titre simulation
# Read in outputs from antibody titre simulation
titre_2 <- read.table("./02.scripts/supporting-figure-creation/titre_2.txt", header = F)
ab_runs_2 <- read.table("./02.scripts/supporting-figure-creation/ab_runs_2.txt", header = F)

# Initialise matrices to store simulations of vaccine efficacy over time
# Initialise matrices to store simulations of vaccine efficacy over time
efficacy_matrix_2 <- matrix(NA_real_, nrow = length(tvec), ncol = 1)
efficacy_matrix_runs_2 <- matrix(NA_real_, nrow = length(tvec), ncol = 100)
c2 <- matrix(NA_real_, nrow = length(tvec), ncol = 2)

# 5 dose sRTS,S efficacy
sim2 <- efficacy_profile(alpha, beta, Vmax, titre_2[, 1, drop = TRUE])

for (i in 1:100) {
  efficacy_matrix_runs_2[, i] <- efficacy_profile(
    alpha, beta, Vmax,
    ab_runs_2[, i, drop = TRUE]
  )
}

c2 <- t(apply(
  efficacy_matrix_runs_2,
  1,
  stats::quantile,
  probs = c(0.025, 0.975),
  na.rm = TRUE
))

quant_2 <- apply(
  efficacy_matrix_runs_2,
  1,
  stats::quantile,
  na.rm = TRUE
)

srtss_efficacy_7dose_bl <- data.frame(
  time = tvec,
  median = sim2, # <- just the vector
  lower_50 = quant_2[2, ], # 25th percentile  (from apply(..., quantile))
  upper_50 = quant_2[4, ], # 75th percentile
  lower = c2[, 1], # 2.5th percentile (95% CI lower)
  upper = c2[, 2], # 97.5th percentile (95% CI upper)
  name = "RTSS 7 dose",
  model = "Model 1"
)

#-MODEL 2 RTSS EFFICACY---------------------------------------------------------
# Returning peak dose efficacy after dose 3 to peak levels

## New parameter----------------------------
ab_mu_boost_m2 <- 621

## - 5 doses Model 2-------------------------

### -Anitbody simulations----------
# Initialise matrices to store simulations of antibody
#  titre and vaccine efficacy over time.
ab_matrix_3 <- matrix(NA_real_, nrow = length(tvec), ncol = reps)
c3 <- matrix(NA_real_, nrow = length(tvec), ncol = 2)

for (i in 1:reps) {
  sim_3 <- antibody_titre_sv(d1_mu, d1_sigma, d2_mu, d2_sigma,
    rho_mu, rho_sigma, ab_mu, ab_sigma,
    rho_mu_boost, rho_sigma_boost,
    t_boost, t_boost1,
    t_boost2 = NULL, t_boost3 = NULL,
    ab_mu_boost_m2, ab_sigma_boost,
    years
  ) # run simulation
  ab_matrix_3[, i] <- sim_3$ab # store values
}

medianAB_3 <- apply(ab_matrix_3, 1, median)
quant_3 <- apply(ab_matrix_3, 1, quantile)

for (j in 1:length(tvec)) {
  c3[j, ] <- quantile(ab_matrix_3[j, ], c(0.025, 0.975))
}

write(medianAB_3, ncolumns = 1, file = "./02.scripts/supporting-figure-creation/titre_3.txt")
write.table(ab_matrix_3[, 1:100], file = "./02.scripts/supporting-figure-creation/ab_runs_3.txt", col.names = FALSE, row.names = FALSE)

### -Efficacy simulations----------
titre_3 <- read.table("./02.scripts/supporting-figure-creation/titre_3.txt", header = F)
ab_runs_3 <- read.table("./02.scripts/supporting-figure-creation/ab_runs_3.txt", header = F)

# Initialise matrices to store simulations of vaccine efficacy over time
# Initialise matrices to store simulations of vaccine efficacy over time
efficacy_matrix_3 <- matrix(NA_real_, nrow = length(tvec), ncol = 1)
efficacy_matrix_runs_3 <- matrix(NA_real_, nrow = length(tvec), ncol = 100)
c3 <- matrix(NA_real_, nrow = length(tvec), ncol = 2)

# 5 dose sRTS,S efficacy
sim3 <- efficacy_profile(alpha, beta, Vmax, titre_3[, 1, drop = TRUE])

for (i in 1:100) {
  efficacy_matrix_runs_3[, i] <- efficacy_profile(
    alpha, beta, Vmax,
    ab_runs_3[, i, drop = TRUE]
  )
}

c3 <- t(apply(
  efficacy_matrix_runs_3,
  1,
  stats::quantile,
  probs = c(0.025, 0.975),
  na.rm = TRUE
))

quant_3 <- apply(
  efficacy_matrix_runs_3,
  1,
  stats::quantile,
  na.rm = TRUE
)

srtss_efficacy_5dose_m2 <- data.frame(
  time = tvec,
  median = sim3, # <- just the vector
  lower_50 = quant_3[2, ], # 25th percentile  (from apply(..., quantile))
  upper_50 = quant_3[4, ], # 75th percentile
  lower = c3[, 1], # 2.5th percentile (95% CI lower)
  upper = c3[, 2], # 97.5th percentile (95% CI upper)
  name = "RTSS 5 dose",
  model = "Model 2"
)

## - 5 doses Model 2-------------------------

### -Anitbody simulations----------
# Initialise matrices to store simulations of antibody
#  titre and vaccine efficacy over time.
ab_matrix_4 <- matrix(NA_real_, nrow = length(tvec), ncol = reps)
c4 <- matrix(NA_real_, nrow = length(tvec), ncol = 2)

for (i in 1:reps) {
  sim_4 <- antibody_titre_sv(
    d1_mu, d1_sigma, d2_mu, d2_sigma,
    rho_mu, rho_sigma, ab_mu, ab_sigma,
    rho_mu_boost, rho_sigma_boost,
    t_boost, t_boost1, t_boost2, t_boost3,
    ab_mu_boost_m2, ab_sigma_boost,
    years
  ) # run simulation
  ab_matrix_4[, i] <- sim_4$ab # store values
}

medianAB_4 <- apply(ab_matrix_4, 1, median)
quant_4 <- apply(ab_matrix_4, 1, quantile)

for (j in 1:length(tvec)) {
  c4[j, ] <- quantile(ab_matrix_4[j, ], c(0.025, 0.975))
}

write(medianAB_4, ncolumns = 1, file = "./02.scripts/supporting-figure-creation/titre_4.txt")
write.table(ab_matrix_4[, 1:100], file = "./02.scripts/supporting-figure-creation/ab_runs_4.txt", col.names = FALSE, row.names = FALSE)

### -Efficacy simulations----------
titre_4 <- read.table("./02.scripts/supporting-figure-creation/titre_4.txt", header = F)
ab_runs_4 <- read.table("./02.scripts/supporting-figure-creation/ab_runs_4.txt", header = F)

# Initialise matrices to store simulations of vaccine efficacy over time
# Initialise matrices to store simulations of vaccine efficacy over time
efficacy_matrix_4 <- matrix(NA_real_, nrow = length(tvec), ncol = 1)
efficacy_matrix_runs_4 <- matrix(NA_real_, nrow = length(tvec), ncol = 100)
c4 <- matrix(NA_real_, nrow = length(tvec), ncol = 2)

# 5 dose sRTS,S efficacy
sim4 <- efficacy_profile(alpha, beta, Vmax, titre_4[, 1, drop = TRUE])

for (i in 1:100) {
  efficacy_matrix_runs_4[, i] <- efficacy_profile(
    alpha, beta, Vmax,
    ab_runs_4[, i, drop = TRUE]
  )
}

c4 <- t(apply(
  efficacy_matrix_runs_4,
  1,
  stats::quantile,
  probs = c(0.025, 0.975),
  na.rm = TRUE
))

quant_4 <- apply(
  efficacy_matrix_runs_4,
  1,
  stats::quantile,
  na.rm = TRUE
)

srtss_efficacy_7dose_m2 <- data.frame(
  time = tvec,
  median = sim4, # <- just the vector
  lower_50 = quant_4[2, ], # 25th percentile  (from apply(..., quantile))
  upper_50 = quant_4[4, ], # 75th percentile
  lower = c4[, 1], # 2.5th percentile (95% CI lower)
  upper = c4[, 2], # 97.5th percentile (95% CI upper)
  name = "RTSS 7 dose",
  model = "Model 2"
)


#-MODEL 3 ENHANCEMENTS---------------------------------------------------------------

## -RTSS------------------------------------
# updated Hill Function parameters
alpha_m3 <- 0.86805 # Shape parameter of dose-response curve
beta_m3 <- 70.9131 # Scale parameter of dose-response curve
Vmax_m3 <- 0.842998 # maximum efficacy against infection

### 5-dose
srtss_efficacy_5dose_m3 <- data.frame(
  time = tvec,
  median = sim3, # <- just the vector
  lower_50 = quant_3[2, ], # 25th percentile  (from apply(..., quantile))
  upper_50 = quant_3[4, ], # 75th percentile
  lower = c3[, 1], # 2.5th percentile (95% CI lower)
  upper = c3[, 2], # 97.5th percentile (95% CI upper)
  name = "RTSS 5 dose",
  model = "Model 3",
  combined = efficacy_profile(alpha_m3, beta_m3, Vmax_m3, titre_3[, 1, drop = TRUE]) # updated synergy efficacy
)

### 7-dose
srtss_efficacy_7dose_m3 <- data.frame(
  time = tvec,
  median = sim4, # <- just the vector
  lower_50 = quant_4[2, ], # 25th percentile  (from apply(..., quantile))
  upper_50 = quant_4[4, ], # 75th percentile
  lower = c4[, 1], # 2.5th percentile (95% CI lower)
  upper = c4[, 2], # 97.5th percentile (95% CI upper)
  name = "RTSS 7 dose",
  model = "Model 3",
  combined = efficacy_profile(alpha_m3, beta_m3, Vmax_m3, titre_4[, 1, drop = TRUE]) # updated synergy efficacy
)

## -SMC SP+AQ-------------------------------

smc_synergy <- data.frame(eff = exp(-(time_seq / 45.76)^2.87))

smc_m3_efficacy <- data.frame(
  time = time_seq,
  median = quant_cs[2, ],
  lower_50 = quant_cs5[1, ],
  upper_50 = quant_cs5[3, ],
  lower = quant_cs[1, ],
  upper = quant_cs[3, ],
  name = "SP+AQ",
  model = "Model 3",
  combined = smc_synergy$eff
)


#-PLOT--------------------------------------------------------------------------
# single dataframe
efficacy_models <-
  bind_rows(
    smc_bl_efficacy,
    smc_m2_efficacy,
    srtss_efficacy_5dose_bl,
    srtss_efficacy_7dose_bl,
    srtss_efficacy_5dose_m2,
    srtss_efficacy_7dose_m2,
    srtss_efficacy_5dose_m3,
    srtss_efficacy_7dose_m3,
    smc_m3_efficacy
  )

efficacy_models <- efficacy_models |>
  mutate(
    name  = factor(name, levels = c("RTSS 5 dose", "RTSS 7 dose", "SP+AQ")),
    model = factor(model, levels = c("Model 1", "Model 2", "Model 3"))
  )

intervention_colors <- c(
  "RTSS 5 dose" = "#94623d",
  "RTSS 7 dose" = "#5489d0",
  "SP+AQ" = "#3f3f3f"
)

# plot
intervention_efficacy_plot <- ggplot(efficacy_models, aes(x = time)) +
  # 95% CrI
  geom_ribbon(
    aes(ymin = lower * 100, ymax = upper * 100, fill = name),
    alpha = 0.15,
    colour = NA
  ) +
  # 50% CrI
  geom_ribbon(
    aes(ymin = lower_50 * 100, ymax = upper_50 * 100, fill = name),
    alpha = 0.35,
    colour = NA
  ) +
  # Median profile
  geom_line(
    aes(y = median * 100, color = name),
    linewidth = 0.6
  ) +
  # Model 3 “combined” enhancement – only in Model 3 panels
  geom_line(
    data = efficacy_models |> dplyr::filter(model == "Model 3"),
    aes(y = combined * 100),
    inherit.aes = TRUE,
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    "Time (days)"
  ) +
  scale_y_continuous(
    "Intervention efficacy (%)",
    breaks = seq(0, 100, 20),
    limits = c(0, 100)
  ) +
  scale_fill_manual(values = intervention_colors) +
  scale_color_manual(values = intervention_colors) +
  facet_grid(model ~ name, scales = "free_x") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "none",
    legend.box        = "horizontal",
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold"),
    strip.background  = element_blank(),
    axis.title.x      = element_text(margin = margin(t = 8)),
    axis.title.y      = element_text(margin = margin(r = 8))
  )

ggsave("./05.plots/publication-plots/main-text/intervention-efficacy.png", plot = intervention_efficacy_plot, width = 8, height = 6, dpi = 600)
ggsave("./05.plots/publication-plots/main-text/intervention-efficacy.pdf", plot = intervention_efficacy_plot, width = 8, height = 6, dpi = 600)
