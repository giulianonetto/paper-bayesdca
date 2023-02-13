library(tidyverse)
x1 <- readRDS("setting1.rds")
system.time({
    fit <- x1$df_pop %>%
        mutate(outcomes = survival::Surv(obsTime, status)) %>%
        select(outcomes, p_hat) %>%
        sample_n(200) %>%
        bayesDCA::dca_surv(
            prediction_time = 12,
            thresholds = seq(0, .9, length = 51),
            keep_prior = TRUE,
            iter = 2000,
            cores = 2,
            chains = 2,
            refresh = 0,
            cutpoints = seq(0, 12, length = 10)
        )
})
fit
plot(fit) +
    geom_point(
        data = x1$true_nb, aes(.thr, nb),
        inherit.aes = TRUE,
        color = "red",
        size = 1
    )


summary(survival::survfit(outcomes ~ 1, data = fit$.data), times = 12)
