# Extended RTS,S Doses Analysis

Code and curated outputs for:

**The public health impact and cost-effectiveness of extended RTS,S malaria vaccine schedules of up to 7 doses in seasonal settings: a mathematical modelling study**

## Repository Contents

- `01.data/`: small input files required by the analysis pipeline.
- `02.scripts/`: R scripts for scenario generation, model runs, post-processing, cost-effectiveness analysis, sensitivity analysis, and supporting figures.
- `03.outputs/paper-results/`: processed analysis outputs used to generate manuscript figures and tables.
- `05.plots/publication-plots/`: final manuscript and supplementary figures/tables, including the supplementary Excel tables.

## Run Order

Run scripts from the repository root in this order:

1. `02.scripts/s1-parameter_draws.R`
2. `02.scripts/s2-baseline-calibration.R`
3. `02.scripts/s3-model-runs.R`
4. `02.scripts/s4-processing-model-runs.R`
5. `02.scripts/s5-results-public-health-impact.R`
6. `02.scripts/s6-results-cost-effectiveness.R`
7. `02.scripts/s7-results-cost-sensitivity.R`
8. `02.scripts/s8-results-5dose-baseline.R`

The full simulation steps are computationally intensive. The curated files in `03.outputs/paper-results/` can be used to regenerate the main paper results without rerunning every individual simulation.

## Large Fitted Model Object

The file `01.data/inputs/uncertainty-params/model.4.RDS` is not included in this GitHub repository because it is approximately 466 MB and is used only by `02.scripts/supporting-figure-creation/s2-vaccine-efficacy-plot.R` to regenerate `05.plots/publication-plots/main-text/intervention-efficacy.png`.

The final generated figure is included in `05.plots/publication-plots/main-text/`. To rerun that supporting figure script, place `model.4.RDS` at:

```text
01.data/inputs/uncertainty-params/model.4.RDS
```

## Software

The analysis is written in R. The `renv.lock` file records the package versions used in this project with R 4.2.1. To restore the package environment, run:

```r
install.packages("renv")
renv::restore()
```

Core packages are loaded in `02.scripts/data_and_libraries.R`, including `malariasimulation`, `cali`, `tidyverse`, `data.table`, `patchwork`, and `MetBrewer`.
