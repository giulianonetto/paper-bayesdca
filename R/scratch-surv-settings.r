library(tidyverse)
library(survival)
source("R/utils.r")
ts <- validate_thresholds(
    thresholds = seq(0, 1, .1)
)

{
    setting <- get_simulation_settings_surv()[[2]]
    d <- simulate_dca_population_surv(
        sim_setting = setting,
        n_pop = 2e5, # ceiling(10000 / setting$event_fraction),
        thresholds = ts,
        .seed = 122445,
        pred_time = 1,
        .verbose = TRUE
    )$df_pop
}
mean(d$survTime == 0)
mean(d$obsTime == 0)
hist(d$survTime, breaks = 200, xlim = c(0, 2.5))
mean(d$survTime < .5)
hist(d$censorTime, breaks = 200, xlim = c(0, 2.5))
hist(d$obsTime, breaks = 200, xlim = c(0, 2.5))
rms::survplot(rms::psm(Surv(obsTime, status) ~ 1, data = d), what = "hazard", xlim = c(1e-5, 2.5))
survfit(Surv(obsTime, status) ~ 1, data = d) %>%
    survminer::ggsurvplot(
        xlim = c(0, 2.2),
        break.time.by = .3,
        surv.median.line = "hv",
        risk.table = "abs_pct",
        cumevents = TRUE,
        cumcensor = TRUE,
        tables.height = 0.15,
    )
median(survfit(Surv(obsTime, status) ~ 1, data = d))
mean(d$status)
mean(d$survTime >= 1)
coxph(Surv(obsTime, status) ~ x1 + x2, data = d) %>%
    summary()


# new interval function

a <- a[a > 0]
.previous <- 0
new_cuts <- c(0)
for (i in seq_along(a)) {
    .current <- a[i]
    events <- sum(event_times > .previous & event_times <= .current)
    if (events >= 3) {
        new_cuts <- c(new_cuts, .current)
        .previous <- .current
    }
}
