# ------------------------------------------------------------------------------
# Purpose: Simulate and visualize vaccine efficacy model trajectories across scenario settings.
# Inputs: Parameter scenario table and simulation controls.
# Outputs: Supporting efficacy trajectory plots for manuscript appendices.
# Dependencies: ggplot2/dplyr/tidyr/purrr used within simulation function.
# Run Stage: Supporting analysis: efficacy model scenario exploration.
# Notes: Comment-only documentation update; logic and outputs unchanged.
# ------------------------------------------------------------------------------

# Function: simulate_vaccine_efficacy
# Purpose: Simulate and optionally plot vaccine efficacy trajectories for scenarios.
# Inputs: Parameter table and simulation controls (years, reps, facets).
# Outputs: Simulation summary data and/or ggplot object.
# Assumptions: Assumes scenario parameter table matches expected columns.

simulate_vaccine_efficacy <- function(
    params_df,                        # Data frame of parameter scenarios
    years = 6,                        # Simulation duration in years
    reps = 5000,                      # Number of repetitions
    return_plot = TRUE,              # Whether to return the plot
    facet_var = c("n_doses", "scenario"),   # Faceting structure
    color_map = NULL,                 # Named vector: scenario -> color
    title_var = NULL
) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)

  # Distributions
  rlogitnorm <- function(n, mean_raw, sd_raw) {
    x <- exp(rnorm(n, mean_raw, sd_raw))
    x / (1 + x)
  }

  rlnorm2 <- function(n, mean, sd) {
    meanlog <- log(mean^2 / sqrt(mean^2 + sd^2))
    sdlog <- sqrt(log(1 + (sd / mean)^2))
    rlnorm(n, meanlog, sdlog)
  }

  # Antibody titre simulation
  antibody_titre_sv <- function(tvec, t_boost_times, ab_mu_vals, params) {
    ab_mu_vals <- as.numeric(ab_mu_vals)  # <-- fix

    rho <- rlogitnorm(1, params$rho_mu, params$rho_sigma)
    rho_boost <- rlogitnorm(1, params$rho_mu_boost, params$rho_sigma_boost)
    d1 <- rlnorm2(1, params$d1_mu, params$d1_sigma)
    d2 <- rlnorm2(1, params$d2_mu, params$d2_sigma)
    r1 <- log(2) / d1
    r2 <- log(2) / d2

    ab <- numeric(length(tvec))
    # First dose (primary)
    boost_vals <- numeric(length(t_boost_times) + 1)
    boost_vals[1] <- exp(rnorm(1, log(params$ab_mu) - params$ab_sigma^2 / 2, params$ab_sigma))

    # Booster doses (if any)
    if (length(t_boost_times) > 0) {
      for (j in seq_along(t_boost_times)) {
        boost_vals[j + 1] <- exp(rnorm(1, log(as.numeric(ab_mu_vals)) - params$ab_sigma_boost^2 / 2, params$ab_sigma_boost))
      }
    }

    time_points <- c(1, t_boost_times + 1)
    for (i in seq_along(time_points)) {
      if (i == length(time_points)) {
        idx <- which(tvec >= time_points[i])
      } else {
        idx <- which(tvec >= time_points[i] & tvec < time_points[i + 1])
      }
      t_offset <- tvec[idx] - tvec[time_points[i]]
      ab[idx] <- boost_vals[i] * (rho_boost * exp(-r1 * t_offset) + (1 - rho_boost) * exp(-r2 * t_offset))
    }
    return(ab)
  }

  # Efficacy from antibody titre
  efficacy_profile <- function(alpha, beta, Vmax, titre) {
    Vmax * (1 - 1 / (1 + (titre / beta)^alpha))
  }

  # Time vector
  tvec <- seq(0, 365 * years, length.out = 365 * years)

  # Loop through each scenario
  results <- purrr::map_dfr(1:nrow(params_df), function(i) {
    row <- params_df[i, ]

    # Ensure n_doses is present
    if (!"n_doses" %in% names(row)) row$n_doses <- 4

    # Automatically space boosters based on booster_spacing
    t_boost_times <- if (row$n_doses > 0) {
      seq(from = row$booster_spacing, by = row$booster_spacing, length.out = row$n_doses)
    } else {
      numeric(0)
    }

    ab_mu_vals <- row$ab_mu_boost

    param_list <- list(
      d1_mu = row$d1_mu, d1_sigma = row$d1_sigma,
      d2_mu = row$d2_mu, d2_sigma = row$d2_sigma,
      rho_mu = row$rho_mu, rho_sigma = row$rho_sigma,
      rho_mu_boost = row$rho_mu_boost, rho_sigma_boost = row$rho_sigma_boost,
      ab_mu = row$ab_mu, ab_sigma = row$ab_sigma,
      ab_sigma_boost = row$ab_sigma_boost
    )

    # Simulate antibody titres
    ab_matrix <- replicate(reps, antibody_titre_sv(tvec, t_boost_times, ab_mu_vals, param_list))
    median_ab <- apply(ab_matrix, 1, median)
    ci_ab <- t(apply(ab_matrix, 1, quantile, probs = c(0.025, 0.25, 0.75, 0.975)))

    # Simulate efficacy
    efficacy_matrix <- apply(
      ab_matrix[, 1:100, drop = FALSE],
      2,
      function(x) efficacy_profile(row$alpha, row$beta, row$Vmax, x)
    )

    median_eff <- apply(efficacy_matrix, 1, median)
    ci_eff <- t(apply(efficacy_matrix, 1, quantile, probs = c(0.025, 0.25, 0.75, 0.975)))
    median_eff <- apply(efficacy_matrix, 1, median)
    ci_eff <- t(apply(efficacy_matrix, 1, quantile, probs = c(0.025, 0.25, 0.75, 0.975)))

    tibble(
      time = tvec,
      median = median_eff,
      lower = ci_eff[, 1],
      upper = ci_eff[, 4],
      lower_50 = ci_eff[, 2],
      upper_50 = ci_eff[, 3],
      scenario = row$scenario,
      n_doses = row$n_doses
    )
  })

  dose_labels <- setNames(
    paste0("Additional doses: ", 0:4),
    as.character(0:4)
  )

  # Plotting
  if (return_plot) {
    plot <- ggplot(results, aes(x = time, y = median, color = scenario, fill = scenario)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA, show.legend = F) +
      geom_ribbon(aes(ymin = lower_50, ymax = upper_50), alpha = 0.5, color = NA, show.legend = F) +
      geom_line(show.legend = F) +
      scale_x_continuous(
        breaks = seq(1, 365 * years, by = 365),
        labels = 1:years,
        limits = c(1, 365 * years)
      ) +
      scale_y_continuous(
        breaks = seq(0, 1, 0.2),
        labels = scales::percent_format(accuracy = 1),
        limits = c(0, 1)
      ) +
      labs(x = "Time since dose three (years)",
           y = "Vaccine Efficacy (%)",
           title = title_var) +
      theme_minimal(base_size = 14) +
      facet_grid(
        rows = vars(n_doses),
        cols = vars(scenario),
        labeller = labeller(n_doses = as_labeller(dose_labels))
      )

    if (!is.null(color_map)) {
      plot <- plot +
        scale_color_manual(values = color_map) +
        scale_fill_manual(values = color_map)
    }

    return(plot)
  } else {
    return(results)
  }
}


param_inputs <- expand.grid(
  scenario = c("model-1", "model-2", "model-3"),
  n_doses = 0:4
) %>%
  dplyr::mutate(
    ab_mu = 621,
    ab_sigma = 0.35,
    ab_sigma_boost = 0.35,
    d1_mu = 45,
    d1_sigma = 16,
    d2_mu = 591,
    d2_sigma = 245,
    rho_mu = 2.37832,
    rho_sigma = 1.00813,
    rho_mu_boost = 1.03431,
    rho_sigma_boost = 1.02735,
    ab_mu_boost = dplyr::case_when(
      scenario %in% c("model-2", "model-3") ~ 621,
      TRUE ~ 277
    ),
    alpha = dplyr::case_when(
      scenario == "model-3" ~ 0.86805,
      TRUE ~ 0.74
    ),
    beta = dplyr::case_when(
      scenario == "model-3" ~ 70.9131,
      TRUE ~ 99.2
    ),
    Vmax = dplyr::case_when(
      scenario == "model-3" ~ 0.842998,
      TRUE ~ 0.93
    ),
    booster_spacing = 365  # default: one year between boosters
  )

simulate_vaccine_efficacy(
  param_inputs,
  years = 6,
  reps = 5000,
  color_map = c("model-1" = "#243E36", "model-2" = "#7CA982", "model-3" = "#BFE37C"),
  title_var = "Simulated Vaccine Efficacy Profile - Seasonal Vaccination"
)

ggsave("./05.plots/fig-3-seasonal-vaccination-models-efficacy.png", width=8.5, height=8.5, dpi=600)


# Age efficacy based model - 15 month gap between 3 and 4
param_inputs2 <- expand.grid(
  scenario = c("model-1"),
  n_doses = 0:1
) %>%
  dplyr::mutate(
    ab_mu = 621,
    ab_sigma = 0.35,
    ab_sigma_boost = 0.35,
    d1_mu = 45,
    d1_sigma = 16,
    d2_mu = 591,
    d2_sigma = 245,
    rho_mu = 2.37832,
    rho_sigma = 1.00813,
    rho_mu_boost = 1.03431,
    rho_sigma_boost = 1.02735,
    ab_mu_boost = dplyr::case_when(
      scenario %in% c("model-2", "model-3") ~ 621,
      TRUE ~ 277
    ),
    alpha = dplyr::case_when(
      scenario == "model-3" ~ 0.86805,
      TRUE ~ 0.74
    ),
    beta = dplyr::case_when(
      scenario == "model-3" ~ 70.9131,
      TRUE ~ 99.2
    ),
    Vmax = dplyr::case_when(
      scenario == "model-3" ~ 0.842998,
      TRUE ~ 0.93
    ),
    booster_spacing = 450
  )

simulate_vaccine_efficacy(
  param_inputs2,
  years = 6,
  reps = 5000,
  color_map = c("model-1" = "#243E36", "model-2" = "#7CA982", "model-3" = "#BFE37C"),
  title_var = "Simulated Vaccine Efficacy Profile - Age-based Vaccination"
)

ggsave("./05.plots/fig-4-age-based-vaccination-models-efficacy.png", width=8.5, height=6, dpi=600)
