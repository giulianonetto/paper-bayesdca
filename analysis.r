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

# Comparison with other packages
thresholds <- seq(0, 1, by = 0.005)

## RMDA data

rmda_data <- rmda_data %>%
    mutate(
        model_predictions = expit(
            -10.5 + 0.22 * Age - 0.01 * Female
                + 0.91 * Smokes + 2.03 * Marker1
                - 1.56 * Marker2
        ),
    )

plot_bdca_vs_rmda(
    dataset = rmda_data,
    outcomes = "Cancer",
    predictor = "model_predictions",
    thresholds = thresholds,
    bootstraps = 2e3
)
ggsave(
    str_path("output/bayes_vs_frequentist/rmda_data.png"),
    width = 8, height = 4.5, dpi = 600
)

# GUSTO-I trial

load(str_path("data/gusto.rda"))

## Small sample size

set.seed(123)
test_ix <- sample(nrow(gusto), 500)
dev_data <- gusto[-test_ix, ]
test_data <- gusto[test_ix, ]


fit_gusto <- glm(
    day30 ~ age + sysbp + pulse + Killip,
    data = dev_data,
    family = binomial
)

test_data$phat1 <- predict(
    fit_gusto,
    newdata = test_data,
    type = "response"
)

test_data %>%
    plot_bdca_vs_rmda(
        outcomes = "day30",
        predictor = "phat1",
        thresholds = thresholds,
        bootstraps = 2e3
    ) +
    theme(legend.position = c(.7, .7))

ggsave(
    str_path("output/bayes_vs_frequentist/gusto.png"),
    last_plot() + theme(legend.position = c(.7, .7)) +
        labs(title = NULL),
    width = 8, height = 4.5, dpi = 600
)


## Reasonable sample size

set.seed(123)
test_ix <- sample(nrow(gusto), 1500)
dev_data <- gusto[-test_ix, ]
test_data_large <- gusto[test_ix, ]


fit_gusto <- glm(
    day30 ~ age + sysbp + pulse + Killip,
    data = dev_data,
    family = binomial
)

test_data_large$phat1 <- predict(
    fit_gusto,
    newdata = test_data_large,
    type = "response"
)

test_data_large %>%
    plot_bdca_vs_rmda(
        outcomes = "day30",
        predictor = "phat1",
        thresholds = thresholds,
        bootstraps = 2e3
    ) +
    theme(legend.position = c(.7, .7))

ggsave(
    str_path("output/bayes_vs_frequentist/gusto-1500.png"),
    width = 8, height = 4.5, dpi = 600
)

# Copare dca

misclass <- summary(
    OptimalCutpoints::optimal.cutpoints(
        X = "phat1",
        status = "day30",
        tag.healthy = 0,
        data = test_data_large,
        methods = c("MCT")
    )
)
t_misclass <- misclass$MCT$Global$optimal.cutoff$cutoff[1]

df <- data.frame(
    outcomes = test_data_large$day30,
    model = test_data_large$phat1,
    classifier = as.numeric(test_data_large$phat1 > t_misclass),
    noisy_model = expit(
        qlogis(test_data_large$phat1) + rnorm(nrow(test_data_large), sd = 1.5)
    )
)

bdca <- dca(df, cores = 4, thresholds = thresholds)

compare_dca(bdca, type = "best")
ggsave(
    str_path("output/bayes_vs_frequentist/gusto-1500-compare-best.png"),
    width = 10, height = 6, dpi = 600
)
compare_dca(bdca, type = "useful")
ggsave(
    str_path("output/bayes_vs_frequentist/gusto-1500-compare-useful.png"),
    width = 10, height = 6, dpi = 600
)
compare_dca(bdca, models_or_tests = c("model", "classifier"), type = "pairwise")
ggsave(
    str_path("output/bayes_vs_frequentist/gusto-1500-compare-pairwise-model-vs-classifier.png"),
    width = 10, height = 6, dpi = 600
)

compare_dca(bdca, models_or_tests = c("model", "noisy_model"), type = "pairwise")
ggsave(
    str_path("output/bayes_vs_frequentist/gusto-1500-compare-pairwise-model-vs-noisy-model.png"),
    width = 10, height = 6, dpi = 600
)

# TODO: fix compare_dca to allow multiple models_or_tests (it's capping at two)
# WRITE
