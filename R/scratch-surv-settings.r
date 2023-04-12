library(tidyverse)
library(survival)

{
    source("R/utils.r")
    setting <- get_simulation_settings_surv()$sim9
    d <- simulate_dca_population_surv(
        sim_setting = setting,
        n_pop = 2e3, # ceiling(10000 / setting$event_fraction),
        thresholds = validate_thresholds(
            thresholds = c(0, .1)
        ),
        .seed = 6698,
        pred_time = 12,
        .verbose = TRUE
    )$df_pop
}

coxph(Surv(obsTime, status) ~ x1 + x2, data = d) %>%
    summary()

fitdistrplus::fitdist(d$survTime, "weibull")
fitdistrplus::fitdistcens(fitdistrplus::Surv2fitdistcens(time = d$obsTime, event = d$status, type = "right"), "weibull")
rms::survplot(rms::psm(Surv(obsTime, status) ~ 1, data = d), what = "hazard")
mean(d$survTime == 0)
mean(d$obsTime == 0)
hist(d$survTime, breaks = 200)
mean(d$survTime >= 12)
hist(d$censorTime, breaks = 200, add = T, col = "red")
hist(d$obsTime, breaks = 200)

survfit(Surv(obsTime, status) ~ 1, data = d) %>%
    survminer::ggsurvplot(
        xlim = c(0, 50),
        break.time.by = 10,
        surv.median.line = "hv",
        risk.table = "abs_pct",
        cumevents = TRUE,
        cumcensor = TRUE,
        tables.height = 0.15,
    )
median(survfit(Surv(obsTime, status) ~ 1, data = d))
mean(d$status)



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
