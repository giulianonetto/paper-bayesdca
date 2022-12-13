library(tidyverse)
theme_set(theme_bw(base_size = 14))
library(bayesDCA)
# prior strength
plot(ts, pmax(50 - 300 * ts + 300 * ts^2, 5))

load(
    here::here("data/gusto.rda")
)
source("R/utils.r")
.seed <- 12112022
set.seed(.seed)
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
thresholds <- seq(0, 1, 0.01)
plot_gusto <- plot_bdca_vs_rmda(
    dataset = test_data,
    outcomes = "day30",
    predictor = "phat1",
    thresholds = thresholds,
    constant_prior = TRUE,
    bootstraps = 2e3
) +
    theme(legend.position = c(.7, .7)) +
    labs(title = NULL)

plot_gusto2 <- plot_bdca_vs_rmda(
    dataset = test_data,
    outcomes = "day30",
    predictor = "phat1",
    thresholds = thresholds,
    constant_prior = FALSE,
    bootstraps = 2e3
) +
    theme(legend.position = c(.7, .7)) +
    labs(title = NULL)

df <- data.frame(
    outcomes = test_data$day30,
    preds = test_data$phat1
)

system.time({
    du0 <- dca(
        df,
        thresholds = thresholds,
        prior_only = TRUE,
        refresh = 1
    )
})
system.time({
    du1 <- dca(
        df,
        thresholds = thresholds,
        prior_only = FALSE,
        refresh = 1
    )
})
system.time({
    di0 <- dca(
        df,
        thresholds = thresholds,
        prior_only = TRUE,
        constant_prior = FALSE,
        prior_sample_size = 25,
        refresh = 1
    )
})
system.time({
    di1 <- dca(
        df,
        thresholds = thresholds,
        prior_only = FALSE,
        constant_prior = FALSE,
        prior_sample_size = 25,
        refresh = 1
    )
})
