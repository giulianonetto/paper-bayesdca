library(tidyverse)
x1 <- readRDS("setting1.rds")
system.time({
    fit <- x1$df_pop %>%
        sample_n(131) %>%
        mutate(outcomes = survival::Surv(obsTime, status)) %>%
        select(outcomes, p_hat) %>%
        bayesDCA::dca_surv(
            keep_fit = T,
            prediction_time = 12,
            thresholds = seq(0, .9, length = 51),
            keep_prior = TRUE,
            cutpoints = seq(0, 12, length = 10),
            prior_only = TRUE,
            iter = 2000,
            cores = 1,
            chains = 1,
            refresh = 0,
            prior_scaling_factor = 2,
            min_events = 0
        )
})
fit

plot(fit) +
    geom_point(
        data = x1$true_nb, aes(.thr, nb),
        inherit.aes = TRUE,
        color = "red",
        size = 1
    ) +
    coord_cartesian(ylim = c(-1, 1))

system.time({
    fitw <- x1$df_pop %>%
        sample_n(131) %>%
        mutate(outcomes = survival::Surv(obsTime, status)) %>%
        select(outcomes, p_hat) %>%
        bayesDCA::dca_surv_weibull(
            prediction_time = 12,
            iter = 3000,
            cores = 2,
            chains = 2,
            refresh = 1,
            thresholds = seq(0, .9, length = 51),
            mean_mu = 0,
            sd_mu = 5,
            mean_log_alpha = 1,
            sd_log_alpha = 1,
            prior_only = FALSE,
            keep_fit = TRUE
        )
})

plot(fitw) +
    geom_point(
        data = x1$true_nb, aes(.thr, nb),
        inherit.aes = TRUE,
        color = "red",
        size = 1
    ) +
    coord_cartesian(ylim = c(-1, 1))
s3 <- rstan::extract(fitw$fit, "St_positives")[[1]][, 1, ]
hist(s3[, 1], breaks = 100, col = "red", add = T)
fitw
plot(fitw) +
    geom_point(
        data = x1$true_nb, aes(.thr, nb),
        inherit.aes = TRUE,
        color = "red",
        size = 1
    ) +
    coord_cartesian(ylim = c(-1, 1))


summary(survival::survfit(outcomes ~ 1, data = fit$.data), times = 12)
