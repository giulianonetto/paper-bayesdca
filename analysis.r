library(tidyverse)
library(here)
import::from(
    "R/utils.r",
    str_path, expit,
    plot_bdca_vs_rmda
)
import::from(rmda, rmda_data = dcaData)
import::from(dcurves, dcurves_dca = dca)
if (!require(bayesDCA)) {
    devtools::install_github("giulianonetto/bayesdca")
    require(bayesDCA)
}
theme_set(
    theme_bw(base_size = 14)
)

# Comparison with other packages ----
thr_interval <- 0.02
thresholds <- seq(0, 1, by = thr_interval)

# Simulation section ----
# Simulation section ----
ggplot2::theme_set(ggplot2::theme_bw(base_size = 14))
outdir <- str_path("output-tmp/simulation_study")
thresholds <- pmin(thresholds, get_max_threshold())

## simulate predictors, true p|x, and y|p
simulation <- simulate_dca_data(
    n_pop = 1e6,
    thresholds = thresholds,
    true_beta = c(-2.1, -log(1.5), log(1.5)),
    beta_hat = c(-3, -log(1.5) * 3, log(1.5) * 3),
    events = 100,
    .seed = .seed,
    .plot = FALSE
)



ggsave(
    str_path("{outdir}/simulation.png"),
    .plot,
    width = 7, height = 4, dpi = 600
)

auc_true_p <- pROC::auc(y, true_p)
auc_phat <- pROC::auc(y, p_hat)
