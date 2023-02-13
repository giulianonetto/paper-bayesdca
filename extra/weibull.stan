data {
    int<lower=0> n_events;
    int<lower=0> n_censored;
    vector[n_events] event_times;
    vector[n_censored] censored_times;
    real<lower=0> sd_log_alpha;
    real<lower=0> sd_mu;
    real mean_mu;
    int<lower=0, upper=1> prior_only;
}

parameters {
    real alpha_raw;
    real mu;
}

transformed parameters {
    real<lower=0> alpha;
    alpha = exp(sd_log_alpha * alpha_raw);
}

model {
    if (prior_only == 0) {
        event_times ~ weibull(alpha, exp(- mu / alpha));
        target += weibull_lccdf(censored_times | alpha, exp(- mu / alpha));
    }

    alpha_raw ~ normal(0.0, 1.0);
    mu ~ normal(mean_mu, sd_mu);
}
