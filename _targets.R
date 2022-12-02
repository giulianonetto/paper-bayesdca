# Load packages required to define the pipeline:
library(targets)
library(here)
source(here("R/utils.r"))
source(here("R/pipeline_functions.r"))
options(tidyverse.quiet = TRUE)

# global seed
.seed <- 12022022

# Set target options:
tar_option_set(
  packages = c(
    "tidyverse", "bayesDCA", "rmda",
    "dcurves", "OptimalCutpoints", "furrr"
  )
)

simulation_dir <- str_path("output/simulation_study")

# Replace the target list below with your own:
list(
  tar_target(
    name = results_01_subsection,
    command = run_bayes_vs_frequentist(
      thresholds = seq(0, 0.9, 0.01),
      .seed = .seed
    )
  ),
  tar_target(
    name = results_02_subsection,
    command = run_simulation_study(
      n_sim = 5000,
      thresholds = c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
      n_pop = 2e6,
      outdir = simulation_dir,
      overwrite = FALSE,
      .workers = 32,
      .seed = .seed,
      .verbose = TRUE
    )
  ),
  tar_target(
    # plotting gets a step of its own to leave simulation alone
    name = results_02_subsection_plots,
    command = plot_simulation_results(
      simulation_results = results_02_subsection,
      outdir = simulation_dir,
      global_simulation_seed = .seed
    )
  ),
  tar_target(
    name = results_03_subsection,
    command = run_case_study(
      thresholds = seq(0, 0.5, 0.01),
      .seed = .seed
    )
  ),
  tar_target(
    name = results_03_subsection_plots,
    command = plot_case_study_results(
      fit = results_03_subsection,
      outdir = str_path("output/case-study")
    )
  )
)
