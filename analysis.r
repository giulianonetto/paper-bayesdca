library(tidyverse)
theme_set(theme_bw(base_size = 14))
library(bayesDCA)

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
thresholds <- seq(0, 0.99, 0.02)
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
