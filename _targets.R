# Load packages required to define the pipeline:
library(targets)
library(here)
source(here("R/utils.r"))
source(here("R/pipeline_functions.r"))
options(tidyverse.quiet = TRUE)

# global seed
.seed <- 20012024

# decision thresholds used in simulations and examples
simulation_thresholds <- c(
  1e-9, 0.001, 0.01,
  0.05, 0.1, 0.25,
  0.5, 0.75
)

# other simulation settings
n_sim <- 1e3
n_pop <- 1e6
workers <- parallel::detectCores() - 1

readr::write_lines(paste0("workers=", workers), "workers.txt")

# Set target options:
tar_option_set(
  packages = c(
    "tidyverse", "bayesDCA", "rmda",
    "dcurves", "OptimalCutpoints", "furrr", "predtools"
  )
)

# Replace the target list below with your own:
list(
  tar_target(
    name = gusto_trial_example,
    command = run_gusto_trial_example(
      outdir = str_path("output/gusto-trial-example"),
      thresholds = seq(0, 0.95, 0.01),
      .seed = .seed
    )
  ),
  tar_target(
    name = simulation_binary_outcomes,
    command = run_simulation_study(
      n_sim = n_sim,
      thresholds = simulation_thresholds,
      n_pop = n_pop,
      outdir = str_path("output/simulation-study-binary"),
      overwrite = FALSE,
      .workers = workers,
      .seed = .seed,
      .verbose = TRUE
    )
  ),
  tar_target(
    name = plot_simulation_binary_outcomes,
    command = plot_simulation_results(
      simulation_results = simulation_binary_outcomes,
      outdir = str_path("output/simulation-study-binary"),
      global_simulation_seed = .seed
    )
  ),
  tar_target(
    name = adnex_case_study,
    command = run_case_study(
      thresholds = seq(0, 0.5, 0.01),
      .seed = .seed
    )
  ),
  tar_target(
    name = plot_adnex_case_study,
    command = plot_case_study_results(
      fit = adnex_case_study,
      outdir = str_path("output/adnex-case-study")
    )
  ),
  tar_target(
    name = simulation_survival_outcomes,
    command = run_simulation_study_surv(
      n_sim = n_sim,
      thresholds = simulation_thresholds,
      n_pop = n_pop,
      pred_time = 12,
      outdir = str_path("output/simulation-study-survival"),
      overwrite = FALSE,
      .workers = workers,
      .seed = .seed,
      .verbose = TRUE
    )
  ),
  tar_target(
    name = plot_simulation_survival_outcomes,
    command = plot_simulation_results(
      simulation_results = simulation_survival_outcomes,
      outdir = str_path("output/simulation-study-survival"),
      surv = TRUE,
      global_simulation_seed = .seed
    )
  ),
  tar_target(
    name = get_final_figure_survival_simulation,
    command = merge_survival_simulation_plots(
      survival_simulation_plots = plot_simulation_survival_outcomes,
      outdir = str_path("output/simulation-study-survival")
    )
  ),
  tar_target(
    name = prior_predictive_checks,
    command = plot_informative_priors_ppc(
      thresholds = seq(0, 0.5, 0.01),
      outdir = str_path("output/informative-priors-ppc")
    )
  ),
  tar_target(
    name = evpi_simulation,
    command = run_evpi_simulation(
      n_sim = n_sim,
      outdir = str_path("output/evpi-simulation")
    )
  ),
  tar_target(
    name = sample_size_simulation,
    command = run_sample_size_simulation(
      n_sim = 300, overwrite = FALSE,
      outdir = str_path("output/sample-size-simulation")
    )
  ),
  tar_target(
    name = expected_regret_simulation,
    command = run_expected_regret_simulation(
      n_sim = 1000, overwrite = FALSE,
      outdir = str_path("output/expected-regret-simulation")
    )
  )
)
