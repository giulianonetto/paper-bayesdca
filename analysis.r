simulate_data <- function(n) {
    true_lp <- rnorm(n, -1.7, sd = 1.05)
    obs_lp <- true_lp / 0.85 + 0.1
    y <- rbinom(n, size = 1, prob = plogis(true_lp))
    return(data.frame(y = y, obs_lp = obs_lp, true_lp = true_lp))
}

{
    print(pROC::auc(pROC::roc(d$y, d$true_lp)))
    print(mean(d$y))
    plot(plogis(d$obs_lp)[1:1e3], plogis(d$true_lp)[1:1e3])
    abline(0, 1, col = "red")
    print(coef(glm(d$y ~ d$obs_lp, family = binomial)))
}

{
    n <- 1e6
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    true_lp <- -2 + 0.9 * x1 + 0.7 * x2 + rnorm(n)
    lp0 <- -1.6 + 0.74 * x1
    lp1 <- -1.8 + 0.8 * x1 + 0.8 * x2
    true_y <- rbinom(n, size = 1, prob = plogis(true_lp))
    y0 <- rbinom(n, size = 1, prob = plogis(lp0))
    y1 <- rbinom(n, size = 1, prob = plogis(lp1))
    data.frame(
        type = c("true_y", "y0", "y1"),
        prev = c(mean(true_y), mean(y0), mean(y1)),
        auc = c(
            pROC::auc(pROC::roc(true_y, true_lp, quiet = TRUE)),
            pROC::auc(pROC::roc(true_y, lp0, quiet = TRUE)),
            pROC::auc(pROC::roc(true_y, lp1, quiet = TRUE))
        )
    ) %>%
        dplyr::mutate_if(is.numeric, round, 2) %>%
        print()
    rms::val.prob(plogis(lp0)[1:1e3], true_y[1:1e3], xlab = "plogis(lp0)")
    rms::val.prob(plogis(lp1)[1:1e3], true_y[1:1e3], xlab = "plogis(lp1)")
}

nn <- 1e5
d <- data.frame(
    outcomes = true_y,
    lp0 = plogis(lp0),
    lp1 = plogis(lp1)
)
f <- bayesDCA::dca(d, thresholds = seq(0, 0.25, 0.01))
plot(f)
bayesDCA::compare_dca(f, type = "useful")
bayesDCA::compare_dca(f, type = "best")
bayesDCA::plot_delta(f, strategies = c("lp1", "lp0"), type = "pairwise") +
    ggplot2::coord_cartesian(ylim = c(-0.01, 0.03))
# take true delta NB for each threshold
# compute frequency of wrong decisions
# compute average NB loss
