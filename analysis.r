library(tidyverse)
library(tidybayes)
library(rstan)
library(survival)

"
source("R/utils.r"); par(mfrow = c(1,3)); sim_setting <- get_simulation_settings_surv()[[1]]; d <- gen_surv_data(shape = sim_setting$shape, scale = sim_setting$scale, true_beta = sim_setting$
    true_beta, max_follow_up = sim_setting$max_follow_up, sample_size = 1e4, seed=NULL); hist(d$survTime); rms::survest(rms::psm(survival::Surv(obsTime, status) ~ 1, data = d), time = 12); median
    (d$survTime); mean(d$survTime >= 12); x <- seq(0, 24, length=200); plot(x, 1 - pweibull(x, shape=sim_setting$shape, scale=sim_setting$scale), ylab="surv", xlab = "months", type='l', lwd = 3, 
    main = "baseline survival"); abline(v=12,h=.1); plot(x, 1 - ecdf(d$survTime)(x), main = "Survival", ylab = "surv", xlab = "months"); abline(v=12,h=.1)
"

# leukemia data
leukemia0 <- leukemia
pred_time <- 75
time_scaling <- pred_time
leukemia$time <- leukemia0$time / time_scaling

times <- seq(0, max(leukemia$time), length = 200)
wb <- rms::survest(rms::psm(Surv(time, status) ~ 1, data = leukemia), times = times)
wb <- data.frame(
    time = wb$time,
    conf.low = wb$lower,
    conf.high = wb$upper,
    estimate = wb$surv
)


standata <- list(
    n_events = sum(leukemia$status == 1L),
    n_censored = sum(leukemia$status == 0L),
    event_times = leukemia$time[leukemia$status == 1],
    censored_times = leukemia$time[leukemia$status == 0],
    sd_log_alpha = 10, # actually sd log alpha
    mean_mu = 5,
    sd_mu = 10,
    prior_only = 0
)
fit <- rstan::stan("extra/weibull.stan", data = standata)
fit

draws <- tidybayes::tidy_draws(fit) %>%
    sample_n(300) %>%
    select(.draw, mu, alpha) %>%
    mutate(
        sigma = exp(-mu / alpha),
        survival_functions = purrr::pmap(
            list(alpha, sigma),
            function(a, s) {
                function(t) {
                    exp(-(t / s)^a)
                }
            }
        ),
        posterior_surv = map(survival_functions, \(surv_fn) tibble(ts = times, s = surv_fn(times)))
    )

draws %>%
    select(.draw, posterior_surv) %>%
    unnest(posterior_surv) %>%
    group_by(ts) %>%
    median_qi(s, .width = c(0.5, 0.8, .95)) %>%
    ggplot(aes(x = ts * time_scaling, y = s, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    scale_fill_brewer() +
    geom_line(
        data = wb,
        aes(x = time * time_scaling, y = estimate),
        color = "red", inherit.aes = FALSE
    ) +
    geom_line(
        data = wb,
        aes(x = time * time_scaling, y = conf.low),
        color = "red", inherit.aes = FALSE, lty = 2
    ) +
    geom_line(
        data = wb,
        aes(x = time * time_scaling, y = conf.high),
        color = "red", inherit.aes = FALSE, lty = 2
    ) +
    ggtitle(
        paste0(
            "pred time: ", pred_time,
            " sd_log_alpha: ", standata$sd_log_alpha,
            " mean_mu: ", standata$mean_mu,
            " sd_mu: ", standata$sd_mu,
            " prior_only: ", standata$prior_only,
            " time scaling: ", time_scaling
        )
    ) +
    geom_vline(
        xintercept = pred_time
    ) +
    theme_bw()



# bladder2 data
bladder2$time <- bladder2$stop - bladder2$start
bladder20 <- bladder2


pred_time <- 35
time_scaling <- pred_time
bladder2$time <- bladder20$time / time_scaling

times <- seq(0, max(bladder2$time), length = 200)
wb <- rms::survest(rms::psm(Surv(time, event) ~ 1, data = bladder2), times = times)
wb <- data.frame(
    time = wb$time,
    conf.low = wb$lower,
    conf.high = wb$upper,
    estimate = wb$surv
)


standata <- list(
    n_events = sum(bladder2$event == 1L),
    n_censored = sum(bladder2$event == 0L),
    event_times = bladder2$time[bladder2$event == 1],
    censored_times = bladder2$time[bladder2$event == 0],
    sd_log_alpha = 10,
    mean_mu = 5,
    sd_mu = 10,
    prior_only = 0
)
fit <- rstan::stan("extra/weibull.stan", data = standata)
fit

draws <- tidybayes::tidy_draws(fit) %>%
    sample_n(300) %>%
    select(.draw, mu, alpha) %>%
    mutate(
        sigma = exp(-mu / alpha),
        posterior_surv = purrr::pmap(
            list(alpha, sigma),
            function(a, s) {
                tibble(ts = times, s = exp(-(times / s)^a))
            }
        )
    )

draws %>%
    select(.draw, posterior_surv) %>%
    unnest(posterior_surv) %>%
    group_by(ts) %>%
    median_qi(s, .width = c(0.5, 0.8, .95)) %>%
    ggplot(aes(x = ts * time_scaling, y = s, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    scale_fill_brewer() +
    geom_line(
        data = wb,
        aes(x = time * time_scaling, y = estimate),
        color = "red", inherit.aes = FALSE
    ) +
    geom_line(
        data = wb,
        aes(x = time * time_scaling, y = conf.low),
        color = "red", inherit.aes = FALSE, lty = 2
    ) +
    geom_line(
        data = wb,
        aes(x = time * time_scaling, y = conf.high),
        color = "red", inherit.aes = FALSE, lty = 2
    ) +
    ggtitle(
        paste0(
            "pred time: ", pred_time,
            " sd_log_alpha: ", standata$sd_log_alpha,
            " mean_mu: ", standata$mean_mu,
            " sd_mu: ", standata$sd_mu,
            " prior_only: ", standata$prior_only,
            " time scaling: ", time_scaling
        )
    ) +
    geom_vline(
        xintercept = pred_time
    ) +
    theme_bw()
