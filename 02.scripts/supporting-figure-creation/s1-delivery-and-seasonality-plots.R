# ------------------------------------------------------------------------------
# Purpose: Create supporting plots for delivery timing, seasonality, and age-pattern diagnostics.
# Inputs: Simulation outputs and run parameter tables from model-run steps.
# Outputs: Supporting figures for manuscript appendix/contextual interpretation.
# Dependencies: data_and_libraries.R plus processed individual simulation outputs.
# Run Stage: Supporting analysis: delivery and seasonality diagnostics.
# Notes: Comment-only documentation update; logic and output paths unchanged.
# ------------------------------------------------------------------------------

source("./02.scripts/data_and_libraries.R")

# Function: plot_rtss_doses_pull_data
# Purpose: Read one simulation output and compute RTS,S dose/SMC summaries for plotting.
# Inputs: x (scenario index), rn (run name).
# Outputs: Returns processed data frame or invisible NULL when file is unavailable.
# Assumptions: Assumes expected columns in individual simulation output files.

plot_rtss_doses_pull_data <- function(x, rn) {

  # Construct file path
  file_path <- paste0("./03.outputs/individual-model-sim-outputs/", format(x, scientific = FALSE), "-", rn, ".rds")

  # Try to read in the file safely
  output <- tryCatch(
    {
      readRDS(file_path)
    },
    error = function(e) {
      message("File not found or error reading: ", file_path)
      return(NULL)  # return NULL if file doesn't exist
    }
  )

  # Exit early if file doesn't exist
  if (is.null(output)) return(invisible(NULL))

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
      mutate(n_smc_treated = 0) |>
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

  # summarize data over the first 15 years (mult of 3 for ITNs)
  sim_length <- output$sim_length[1] / 365

  # Filter timesteps up to the desired simulation length in years
  output_filtered <- output %>%
    filter(year <= sim_length)

  # Pivot the dose columns into long format
  dose_data <- output_filtered %>%
    select(
      ID,
      runname,
      pfpr,
      seasonality,
      SMC,
      RTSS,
      dosing_assumption,
      coverage_assumption,
      seasonal_ages,
      vaccine_model,
      month,
      vaccine_model,
      month,
      starts_with("dose")
    ) %>%
    pivot_longer(
      cols = starts_with("dose"),
      names_to = "dose_number",
      values_to = "n_vaccinated"
    ) |>
    group_by(
      ID,
      runname,
      pfpr,
      seasonality,
      SMC,
      RTSS,
      dosing_assumption,
      coverage_assumption,
      seasonal_ages,
      vaccine_model,
      month,
      dose_number
    ) |>
    summarise(n_vaccinated = sum(n_vaccinated, na.rm = TRUE), .groups = 'drop')

  # Convert dose_number to a factor with proper ordering
  dose_data$dose_number <- factor(dose_data$dose_number, levels = c("dose7", "dose6", "dose5", "dose4", "dose3", "dose2", "dose1"))

  return(dose_data)

}



#-seasonality-------------------------------------------------------------------
# Establish the length of time, in daily time steps, over which to simulate:
year <- 365
years <- 8
month <- 30
sim_length <- years * year

# Set the size of the human population and an initial entomological inoculation rate (EIR)
human_population <- 100000
starting_EIR <- 50

# highly seasonal
simparams_hs <- get_parameters(
  list(
    human_population = human_population,
    clinical_incidence_rendering_min_ages = 0,
    clinical_incidence_rendering_max_ages = 5 * year,
    model_seasonality = TRUE,
    g0 = 0.284596,
    g = c(-0.317878,-0.0017527,0.116455),
    h = c(-0.331361,0.293128,-0.0617547)
  )
)

simparams_hs <- set_equilibrium(simparams_hs, starting_EIR)

# seasonal
simparams_s <- get_parameters(
  list(
    human_population = human_population,
    model_seasonality = TRUE,
    clinical_incidence_rendering_min_ages = 0,
    clinical_incidence_rendering_max_ages = 5 * year,
    g0 = 0.285505,
    g = c(-0.325352,-0.0109352,0.0779865),
    h = c(-0.132815,0.104675,-0.013919)
  )
)

simparams_s <- set_equilibrium(simparams_s, starting_EIR)

# run model
highly_seas_out <- run_simulation(sim_length, simparams_hs)
seas_out <- run_simulation(sim_length, simparams_s)

# add ID vars to these runs
highly_seas_out$seas <- "Highly Seasonal"
seas_out$seas <- "Seasonal"

# combine mosq values == rainfall
# total mosquito population through time:
highly_seas_out$mosq_total =
  highly_seas_out$Sm_gamb_count +
  highly_seas_out$Im_gamb_count +
  highly_seas_out$Pm_gamb_count

seas_out$mosq_total =
  seas_out$Sm_gamb_count +
  seas_out$Im_gamb_count +
  seas_out$Pm_gamb_count

highly_seas_out$under_5_inc =
  highly_seas_out$n_inc_clinical_0_1825 / highly_seas_out$n_age_0_1825

seas_out$under_5_inc =
  seas_out$n_inc_clinical_0_1825 / seas_out$n_age_0_1825

out_out <-
  bind_rows(highly_seas_out, seas_out)

# get the peak in seasona
peak_hs <- peak_season_offset(simparams_hs)
peak_s <- peak_season_offset(simparams_s)

# Schedule drug administration times (in daily time steps) before, during and after the seasonal peak:
admin_days_hs <- c( -30, 0, 30, 60)
admin_days_s <- c(-60, -30, 0, 30, 60)

# Use the peak-offset, number of simulation years and number of drug admin. days to calculate the days
# on which to administer drugs in each year
smc_days_hs <- rep((365 * seq(1, years-1, by = 1)), each = length(admin_days_hs)) + peak_hs + rep(admin_days_hs, 2)
smc_days_s <- rep((365 * seq(1, years-1, by = 1)), each = length(admin_days_s)) + peak_s + rep(admin_days_s, 2)

# Create a dataframe for the SMC days with appropriate facet labels
smc_df <- bind_rows(
  data.frame(timestep = smc_days_hs, seas = "Highly Seasonal"),
  data.frame(timestep = smc_days_s, seas = "Seasonal")
)

# Plot with vertical lines
ggplot(out_out, aes(x = timestep, y = under_5_inc)) +
  geom_area(alpha=0.2) +
  facet_wrap(~seas, ncol=1) +
  geom_vline(data = smc_df, aes(xintercept = timestep), linetype = "dashed", color = "#56B4E9", linewidth=0.9) +
  theme_minimal(14) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(title = "Seasonality and SMC delivery schedules",
       x = "Time (days)",
       y = "Clinical Incidence (under 5's)")

ggplot(out_out, aes(x = timestep, y = under_5_inc)) +
  geom_area(fill = "grey20", alpha = 0.9) +  # Solid dark area
  geom_vline(data = smc_df, aes(xintercept = timestep),
             linetype = "dashed", color = "#008C9B", linewidth = 0.5, inherit.aes = FALSE) +
  facet_wrap(~seas, ncol = 2, scales="free") +
  labs(
    x = "Time (days)",
    y = "Clinical Incidence (under 5's)"
  ) +
  # theme_classic(14)+
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # panel.grid = element_blank(),
    # axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    # axis.line.x = element_blank(),
    # axis.title = element_text(size = 8),
    strip.background = element_blank()#,
    # strip.text = element_text(size = 8)
  )

ggsave("./05.plots/publication-plots/appendix/smc-delivery-and-seasonality.png", width=6, height=3, dpi=600)

# Plot vaccine timing for seasonal doses
first_hs <- round(1 + (peak_hs - month * 3.5), 0)
first_s  <- round(1 + (peak_s - month * 5.5), 0)

doses_hs <- c(first_hs, first_hs+30, first_hs+60, first_hs+425, first_hs+790)
doses_s <- c(first_s, first_s+30, first_s+60, first_s+425, first_s+790)

vax_df <- bind_rows(
  data.frame(timestep = doses_hs, seas = "Highly Seasonal"),
  data.frame(timestep = doses_s,  seas = "Seasonal")
)

ggplot(out_out, aes(x = timestep, y = under_5_inc)) +
  geom_area(alpha=0.2) +
  facet_wrap(~seas, ncol=1) +
  # Add vaccine dose lines
  geom_vline(data = vax_df, aes(xintercept = timestep),
             linetype = "dotted", color = "#009E73", linewidth = 0.9) +
  theme_minimal(14) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank()) +
  labs(title = "Seasonality and Vaccine Dose Timings",
       x = "Time (days)",
       y = "Clinical Incidence (under 5's)")

ggsave("./05.plots/fig-2-vacc-delivery-and-seasonality.png", width=6, height=7, dpi=600)

#-Paper figures of vaccine actual delivery and seasonality----------------------
# select for seasonal and highly seasonal setting
# optimal coverage
# 7 doses delivered
paramlist <- readRDS("./01.data/vaccine-scenarios/run_parameters.rds")
x <- c(1:nrow(readRDS("./01.data/vaccine-scenarios/runscenarios.rds"))) # number of runs

# define all combinations of scenarios and draws - only run draw 0
# index <- tibble(x = x)
# Filter for drawID == 0 and get the row numbers
rows_to_run <- which(paramlist$drawID == 0)
runname_to_run <- paramlist$run_name[rows_to_run]

# Create index tibble
index <- tibble(x = rows_to_run, rn = runname_to_run) |>
  filter(rn %in% c(
    "0.65-highly seasonal-0-Seasonal-7-doses-S1-optimal-campaigns-5-months-1-year-smc-0.75-model-3-with-age-cap",
    "0.65-seasonal-0-Seasonal-7-doses-S1-optimal-campaigns-5-months-1-year-smc-0.75-model-3-with-age-cap",
    "0.65-highly seasonal-0-Hybrid-7-doses-S1-optimal-campaigns-age-based-primary-smc-0.75-model-3-with-age-cap",
    "0.65-seasonal-0-Hybrid-7-doses-S1-optimal-campaigns-age-based-primary-smc-0.75-model-3-with-age-cap"
  ))

library(furrr)

month <- year / 12

combo <- future_map2_dfr(index$x, index$rn, plot_rtss_doses_pull_data, .progress = TRUE)
combo$timestep <- (combo$month * month) - month
combo$seasonality <- str_to_title(combo$seasonality)

combo$dose_number <- factor(combo$dose_number, levels = c("dose7", "dose6", "dose5", "dose4", "dose3", "dose2", "dose1"))

ggplot(combo, aes(x = timestep, y = n_vaccinated, fill = dose_number)) +
  geom_bar(stat = "identity") +
  facet_wrap(~runname, ncol=1) +
  scale_fill_brewer(
    palette = "Set3",
    name = "Dose",
    breaks = c("dose1", "dose2", "dose3", "dose4", "dose5", "dose6", "dose7"),
    guide = guide_legend(nrow = 1)
  ) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    title = "RTS,S Doses Administered Over Time",
    x = "Timestep (months)",
    y = "Number of Doses Administered"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

out_out <-
  out_out |>
  crossing(tibble(RTSS = c("Hybrid", "Seasonal")))

out_out$seasonality <- str_to_title(out_out$seas)

# plot over the seasonality
ggplot(out_out, aes(x = timestep, y = under_5_inc*10000000)) +
  geom_area(alpha=0.2) +
  geom_bar(combo,
           mapping=aes(x = timestep, y = n_vaccinated, fill = dose_number),
           # col = "black",
           stat = "identity",
           width=month) +
  facet_grid(RTSS~seasonality, scales="free_y") +
  theme_minimal() +
  labs( x = "Time (days)",
       y="") +
  coord_cartesian(xlim=c(0,2000)) +
  scale_fill_brewer(palette = "Spectral",
                    name = "",
                    breaks = c("dose1", "dose2", "dose3", "dose4", "dose5", "dose6", "dose7"),
                    guide = guide_legend(nrow = 1),
                    direction = -1)+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "top")

ggsave("./05.plots/publication-plots/appendix/vacc-delivery-and-seasonality.png", width=8, height=6, dpi=600)
