# Load packages required to define the pipeline:
library(targets)
library(here)
source(here("R/utils.r"))
source(here("R/pipelines.r"))
options(tidyverse.quiet = TRUE)

# global seed
.seed <- 12112022

# Set target options:
tar_option_set(
  packages = c(
    "tidyverse", "bayesDCA", "rmda",
    "dcurves", "OptimalCutpoints"
  )
)

simulation_dir <- str_path("output/simulation_study")
thresholds <- seq(0, 0.9, 0.05)

# Replace the target list below with your own:
list(
  tar_target(
    name = results_01_subsection,
    command = run_bayes_vs_frequentist(
      thresholds = thresholds,
      .seed = .seed
    )
  ),
  tar_target(
    name = results_02_subsection,
    command = run_simulation_study(
      n_sim = 50,
      thresholds = results_01_subsection$thresholds,
      n_pop = 1e6,
      outdir = simulation_dir,
      overwrite = FALSE,
      .seed = .seed
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
  )
)
