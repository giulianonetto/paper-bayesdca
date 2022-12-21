library(tidyverse)
library(simstudy)
library(survival)

gen_surv_data <- function(surv_formula,
                          surv_shape,
                          surv_scale,
                          censor_shape,
                          censor_scale,
                          sample_size) {
    def_x <- defData(
        varname = "x1",
        dist = "exponential",
        formula = 1
    ) %>%
        defData(
            varname = "x2",
            dist = "exponential",
            formula = 1
        )
    def_surv <- defSurv(
        varname = "survTime",
        formula = surv_formula,
        scale = surv_shape,
        shape = surv_scale
    ) %>%
        defSurv(
            varname = "censorTime",
            shape = censor_shape,
            scale = censor_scale
        )

    surv_data <- genData(
        n = sample_size,
        def_x
    ) %>%
        genSurv(
            def_surv,
            timeName = "obsTime",
            censorName = "censorTime",
            eventName = "status",
            keepEvents = TRUE
        )
    return(surv_data)
}

# setting 1
x <- lapply(1:500, \(...) {
    df <- gen_surv_data(
        surv_formula = "log(1.4)*x1 + log(0.6)*x2",
        surv_shape = 1,
        surv_scale = 1,
        censor_shape = 1,
        censor_scale = 0.125,
        sample_size = 1000
    )
    list(
        events = sum(df$status),
        cind = concordance(survival::coxph(Surv(obsTime, status) ~ x1 + x2, data = df)),
        surv = median(survfit(Surv(obsTime, status) ~ 1, data = df))
    )
})
mean(sapply(x, \(i) i$events))
mean(sapply(x, \(i) i$cind[[1]]))
mean(sapply(x, \(i) i$surv[[1]]), na.rm = T)


# setting 2
x <- lapply(1:500, \(...) {
    df <- gen_surv_data(
        surv_formula = "log(3)*x1 + log(1/3)*x2",
        surv_shape = 1,
        surv_scale = 0.3,
        censor_shape = 1,
        censor_scale = 0.3,
        sample_size = 1000
    )
    list(
        events = sum(df$status),
        cind = concordance(survival::coxph(Surv(obsTime, status) ~ x1 + x2, data = df)),
        surv = median(survfit(Surv(obsTime, status) ~ 1, data = df))
    )
})
mean(sapply(x, \(i) i$events))
mean(sapply(x, \(i) i$cind[[1]]))
mean(sapply(x, \(i) i$surv[[1]]), na.rm = T)
