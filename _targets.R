# Load packages required to define the pipeline:
library(targets)
library(here)
source(here("R/utils.r"))
source(here("R/pipelines.r"))
options(tidyverse.quiet = TRUE)

# global seed
.seed <- 11112022

# Set target options:
tar_option_set(
  packages = c(
    "tidyverse", "bayesDCA", "rmda",
    "dcurves", "OptimalCutpoints"
  )
)

# Replace the target list below with your own:
list(
  tar_target(
    name = results_01_subsection,
    command = run_bayes_vs_frequentist(
      thr_interval = 0.05,
      .seed = .seed
    )
  ),
  tar_target(
    name = results_02_subsection,
    command = run_simulation_study(
      thresholds = results_01_subsection$thresholds,
      .seed = .seed
    )
  )
)
