simulate_data <- function(n) {
    true_lp <- rnorm(n, -1.7, sd = 1.05)
    obs_lp <- true_lp / 0.85 + 0.1
    y <- rbinom(n, size = 1, prob = plogis(true_lp))
    return(data.frame(y = y, obs_lp = obs_lp, true_lp = true_lp))
}

{
    d <- simulate_data(1e5)
    print(pROC::auc(pROC::roc(d$y, d$true_lp)))
    print(mean(d$y))
    plot(plogis(d$obs_lp)[1:1e3], plogis(d$true_lp)[1:1e3])
    abline(0, 1, col = "red")
    print(coef(glm(d$y ~ d$obs_lp, family = binomial)))
}
