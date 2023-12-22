library(tidyverse)
library(survival)

{
    source("R/utils.r")
    setting <- get_simulation_settings_surv()$sim9
    s <- simulate_dca_population_surv(
        sim_setting = setting,
        n_pop = 2e3, # ceiling(10000 / setting$event_fraction),
        thresholds = validate_thresholds(
            thresholds = c(0, .1)
        ),
        .seed = 6698,
        pred_time = 12,
        .verbose = TRUE
    )
    d <- s$df_pop
}

df <- d[1:1000, ]
df$outcomes <- Surv(df$obsTime, event=df$status)
df$predictor <- df$p_hat

dcurves_fit <- try(
    {
      dcurves::dca(
        outcomes ~ predictor,
        data = df,
        time = pred_time,
        thresholds = thresholds
      )
    },
    silent = TRUE
  )
dcurves_fit$dca
B <- 200
res <- replicate(
    B,
    {
        ix <- sample(nrow(df), replace = TRUE)
        .fit <- dcurves::dca(
            outcomes ~ predictor,
            data = df[ix, ],
            time = pred_time,
            thresholds = thresholds
        )
        .fit$dca %>% 
            dplyr::select(variable, label, threshold, net_benefit)
    },
    simplify = FALSE
)
res <- dplyr::bind_rows(res, .id = "id")

res %>% 
    group_by(variable, label, threshold) %>% 
    summarise(se = sd(net_benefit))
