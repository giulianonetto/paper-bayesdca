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

simulation_output_dir <- str_path("output/simulation_study/tmp")
dir.create(
  simulation_output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# Replace the target list below with your own:
list(
  tar_target(
    name = results_01_subsection,
    command = run_bayes_vs_frequentist(
      thresholds = c(0.1, 0.2, 0.3),
      .seed = .seed
    )
  ),
  tar_target(
    name = results_02_subsection,
    command = run_simulation_study(
      n_sim = 40,
      thresholds = results_01_subsection$thresholds,
      n_pop = 1e4,
      output_dir = simulation_output_dir,
      overwrite = FALSE,
      .seed = .seed
    )
  )
)
