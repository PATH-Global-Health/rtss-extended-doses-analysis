# ------------------------------------------------------------------------------
# Purpose: Load and centralize R package dependencies used across analysis scripts.
# Inputs: None at runtime; package installation handled externally.
# Outputs: Attached packages for downstream scripts in this repository.
# Dependencies: R package library, including malariasimulation and tidyverse stack.
# Run Stage: Shared dependency loader sourced by pipeline scripts.
# Notes: Comment-only documentation update; package calls unchanged.
# ------------------------------------------------------------------------------

# Libraries --------------------------------------------------------------------

# # Install modified version
# devtools::install_github("ht1212/malariasimulation-vacc-cap", ref = "add-vaccine-age-cap")
library(malariasimulation)
# OLD VERSION devtools::install_github("mrc-ide/malariasimulation@03cd3d5", force=T)
library(cali)              # devtools::install_github("https://github.com/mrc-ide/cali", force = T)
library(tidyverse)
library(data.table)

# plotting
library(patchwork)
library(MetBrewer)      # devtools::install_github("BlakeRMills/MetBrewer")

