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
    # rms::val.prob(plogis(lp0), true_y, xlab = "plogis(lp0)")
    # rms::val.prob(plogis(lp1), true_y, xlab = "plogis(lp1)")
    d <- data.frame(outcomes = true_y, phat0 = plogis(lp0), phat1 = plogis(lp1), true_p = plogis(true_lp))
    f <- bayesDCA::dca(d, thresholds = c(0.1, 0.15, 0.2, 0.25, 0.3))
    p1 <- bayesDCA::plot_delta(f, strategies = c("true_p", "phat1"), type = "pairwise") +
        ggplot2::coord_cartesian(ylim = c(-0.01, 0.03))
    p2 <- bayesDCA::plot_delta(f, strategies = c("phat1", "phat0"), type = "pairwise") +
        ggplot2::coord_cartesian(ylim = c(-0.01, 0.03))
    ggpubr::ggarrange(p1, p2, ncol = 2)
}

{
    n <- 1e5
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    x3 <- rnorm(n)
    x4 <- rnorm(n)
    true_lp <- -2.5 + 0.5 * x1 + -1 * x2 + -1.5 * x3 + 1.7 * x4
    lp0 <- -1.7 + 0.7 * x1 + 0.6 * x1^2 + -0.75 * x2
    lp1 <- -1.6 + 0.35 * x1 + -1 * x3
    y <- rbinom(n, size = 1, prob = plogis(true_lp))
    d <- data.frame(outcomes = y, phat0 = plogis(lp0), phat1 = plogis(lp1), true_p = plogis(true_lp))
    f <- bayesDCA::dca(d, thresholds = seq(0, 0.3, 0.005))
    .calib_phat0 <- get_calibration_binary(f$.data$outcomes, f$.data$phat0)
    .calib_phat1 <- get_calibration_binary(f$.data$outcomes, f$.data$phat1)
    data.frame(
        model = c("phat0", "phat1"),
        true_prevalence = rep(mean(f$.data$outcomes), 2),
        average_prediction = c(mean(f$.data$phat0), mean(f$.data$phat1)),
        calibration_oe = c(.calib_phat0$oe, .calib_phat1$oe),
        calibration_slope = c(.calib_phat0$slope, .calib_phat1$slope),
        auc = c(
            pROC::auc(pROC::roc(f$.data$outcomes, f$.data$phat0, quiet = TRUE)),
            pROC::auc(pROC::roc(f$.data$outcomes, f$.data$phat1, quiet = TRUE))
        )
    ) %>%
        dplyr::mutate_if(is.numeric, round, 2) %>%
        print()

    p2 <- bayesDCA::plot_delta(f, strategies = c("phat1", "phat0"), type = "pairwise") +
        ggplot2::coord_cartesian(ylim = c(-0.04, 0.1))
    p3 <- bayesDCA::plot_delta(f, strategies = c("phat0", "phat1"), type = "useful") +
        ggplot2::coord_cartesian(ylim = c(-0.04, 0.1)) +
        ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.25, 0.75))
    ggpubr::ggarrange(plot(f) + ggplot2::theme(legend.position = "top"), p2, p3, ncol = 3, nrow = 1)
    print(cor(lp0, lp1))
}

{
    n <- 1e5
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    x3 <- rnorm(n)
    x4 <- rnorm(n)
    # true_lp <- -2.5 + 0.5 * x1 + -1 * x2 + -1.5 * x3 + 1.7 * x4
    # lp0 <- -1.7 + 0.7 * x1 + 0.6 * x1^2 + -0.75 * x2
    # lp1 <- -1.6 + 0.35 * x1 + -1.1 * x3
    true_lp <- -2.5 + 0.5 * x1 + -0.5 * x2 + 1.5 * x3 + -1.7 * x4
    lp0 <- -1.85 + .5 * x1 + -.5 * x2
    lp1 <- -1.7 + 0.4 * x1 + 1.1 * x3
    y <- rbinom(n, size = 1, prob = plogis(true_lp))
    d <- data.frame(outcomes = y, phat0 = plogis(lp0), phat1 = plogis(lp1), true_p = plogis(true_lp))
    f <- bayesDCA::dca(d, thresholds = seq(0, 0.3, 0.005))
    .calib_phat0 <- get_calibration_binary(f$.data$outcomes, f$.data$phat0)
    .calib_phat1 <- get_calibration_binary(f$.data$outcomes, f$.data$phat1)
    data.frame(
        model = c("phat0", "phat1"),
        true_prevalence = rep(mean(f$.data$outcomes), 2),
        average_prediction = c(mean(f$.data$phat0), mean(f$.data$phat1)),
        calibration_oe = c(.calib_phat0$oe, .calib_phat1$oe),
        calibration_slope = c(.calib_phat0$slope, .calib_phat1$slope),
        auc = c(
            pROC::auc(pROC::roc(f$.data$outcomes, f$.data$phat0, quiet = TRUE)),
            pROC::auc(pROC::roc(f$.data$outcomes, f$.data$phat1, quiet = TRUE))
        )
    ) %>%
        dplyr::mutate_if(is.numeric, round, 2) %>%
        print()

    p2 <- bayesDCA::plot_delta(f, strategies = c("phat1", "phat0"), type = "pairwise") +
        ggplot2::coord_cartesian(ylim = c(-0.04, 0.1))
    p3 <- bayesDCA::plot_delta(f, strategies = c("phat0", "phat1"), type = "useful") +
        ggplot2::coord_cartesian(ylim = c(-0.04, 0.1)) +
        ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.25, 0.75))
    print(cor(lp0, lp1))
    ggpubr::ggarrange(plot(f) + ggplot2::theme(legend.position = "top"), p2, p3, ncol = 3, nrow = 1)
}
rms::val.prob(plogis(lp0)[1:1e5], y[1:1e5], xlab = "plogis(lp0)")
rms::val.prob(plogis(lp1)[1:1e5], y[1:1e5], xlab = "plogis(lp1)")

rms::val.prob(plogis(lp0), true_y, xlab = "plogis(lp0)")
rms::val.prob(plogis(lp1), true_y, xlab = "plogis(lp1)")

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


d <- read_tsv("output/sample-size-simulation//simulation_results_sample_size.tsv") %>%
    filter(
        # str_detect(comparison, "max"),
        map_lgl(threshold, ~ any(near(.x, seq(0.1, 0.3, by = 0.05))))
    ) %>%
    mutate(
        comparison = factor(
            comparison,
            levels = c(
                "phat0 vs max(treat all, treat none)",
                "phat1 vs max(treat all, treat none)",
                "true_p vs max(treat all, treat none)",
                "phat1 vs phat0",
                "true p vs phat1",
                "true p vs phat0"
            )
        )
    )
d %>%
    ggplot(aes(factor(threshold), prob, color = ordered(sample_size))) +
    ggbeeswarm::geom_quasirandom(dodge.width = 0.9) +
    facet_wrap(~comparison) +
    labs(
        x = "Threshold",
        y = "P(superior)",
        color = "Sample size"
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    geom_boxplot(
        alpha = 0
    ) +
    theme(legend.position = "top")

d %>%
    filter(
        str_detect(comparison, "max"),
        threshold %in% seq(0.1, 0.3, by = 0.05)
    ) %>%
    ggplot(aes(factor(threshold), delta, color = ordered(sample_size))) +
    ggbeeswarm::geom_quasirandom(dodge.width = 0.5) +
    facet_wrap(~comparison) +
    labs(
        x = "Threshold",
        y = latex2exp::TeX("$\\Delta_\\useful$"),
        color = "Sample size"
    )


{
    get_sample_size_dca <- function(n, include_true_p = FALSE) {
        # Population prevalence: 0.2
        ## baseline model
        ### avg pred 0.22
        ### calib o/e 0.92
        ### calib slope 0.8
        ### concordance 0.62
        ## new model
        ### avg pred 0.22
        ### calib o/e 0.93
        ### calib slope 0.86
        ### concordance 0.77
        x1 <<- rnorm(n)
        x2 <<- rnorm(n)
        x3 <<- rnorm(n)
        x4 <<- rnorm(n)
        true_lp <- -3 + 0.1 * x1 + 1.4 * x2 + -1.75 * x3 + 2 * x4
        lp0 <<- -1.6 + 0.12 * x1 + 1.55 * x2
        lp1 <<- -1.7 + 0.075 * x1 + -1.1 * x3
        y <<- rbinom(n, size = 1, prob = plogis(true_lp))
        d <- data.frame(outcomes = y, phat0 = plogis(lp0), phat1 = plogis(lp1))
        if (isTRUE(include_true_p)) {
            d$true_p <- plogis(true_lp)
        }
        f <- bayesDCA::dca(d, thresholds = seq(0, 0.3, 0.01))
        return(f)
    }

    run_expected_regret_simulation_single_run <- function(n, return_dca = FALSE) {
        # simulate data
        f <- get_sample_size_dca(n = n)
        # scenario 1: phat0 is default, needs to decide whether to implement phat1

        # compute run results for probabilities
        probs <- bayesDCA:::get_superiority_prob_plot_data_binary(f) %>%
            dplyr::select(threshold, decision_strategy, prob)
        # compute run results for deltas
        deltas <- bayesDCA:::get_delta_plot_data_binary(f) %>%
            dplyr::select(threshold, decision_strategy, delta := estimate)
        # nb
        nb <- f$summary$net_benefit %>%
            dplyr::select(threshold, decision_strategy := decision_strategy_name, nb := estimate)

        # combine results
        run_results <- dplyr::inner_join(
            probs,
            deltas,
            by = c("threshold", "decision_strategy")
        ) %>%
            dplyr::inner_join(
                nb,
                by = c("threshold", "decision_strategy")
            ) %>%
            dplyr::mutate(sample_size = n)

        if (isTRUE(return_dca)) {
            return(list(res = run_results, dca = f))
        }

        return(run_results)
    }
    f <- get_sample_size_dca(n = 1e6)
    .calib_phat0 <- get_calibration_binary(f$.data$outcomes, f$.data$phat0)
    .calib_phat1 <- get_calibration_binary(f$.data$outcomes, f$.data$phat1)
    data.frame(
        model = c("phat0", "phat1"),
        true_prevalence = rep(mean(f$.data$outcomes), 2),
        average_prediction = c(mean(f$.data$phat0), mean(f$.data$phat1)),
        calibration_oe = c(.calib_phat0$oe, .calib_phat1$oe),
        calibration_slope = c(.calib_phat0$slope, .calib_phat1$slope),
        auc = c(
            pROC::auc(pROC::roc(f$.data$outcomes, f$.data$phat0, quiet = TRUE)),
            pROC::auc(pROC::roc(f$.data$outcomes, f$.data$phat1, quiet = TRUE))
        )
    ) %>%
        dplyr::mutate_if(is.numeric, round, 2) %>%
        print()

    p2 <- bayesDCA::plot_delta(f, strategies = c("phat1", "phat0"), type = "pairwise") +
        ggplot2::coord_cartesian(ylim = c(NA, 0.05))
    p3 <- bayesDCA::plot_delta(f, strategies = c("phat0", "phat1"), type = "useful") +
        ggplot2::coord_cartesian(ylim = c(-0.04, 0.1)) +
        ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.25, 0.75))
    # n_sim <- 200
    # d <- purrr::map_df(
    #     seq_len(n_sim),
    #     ~ run_expected_regret_simulation_single_run(n = 500),
    #     .id = "run_id"
    # )
    # p4 <- d %>%
    #     pivot_wider(
    #         id_cols = c(threshold, run_id),
    #         names_from = decision_strategy,
    #         values_from = nb
    #     ) %>%
    #     filter(
    #         purrr::map_lgl(threshold, ~ any(dplyr::near(.x, seq(0.1, 0.3, by = 0.05))))
    #     ) %>%
    #     ggplot(aes(phat0, phat1, fill = ordered(threshold))) +
    #     geom_point(size = 2, alpha = 0.5, color = "black", pch = 21) +
    #     geom_abline() +
    #     geom_smooth(
    #         method = "lm", se = FALSE, aes(color = ordered(threshold)),
    #         show.legend = FALSE
    #     ) +
    #     ggpubr::stat_cor(aes(color = ordered(threshold), fill = NULL),
    #         show.legend = FALSE, geom = "label", label.x = .12, label.y = seq(.10, .03, length = 5)
    #     ) +
    #     labs(fill = NULL)
    print(cor(lp0, lp1))
    ggpubr::ggarrange(
        plot(f) + ggplot2::theme(legend.position = "top"),
        p4,
        p2, p3,
        widths = c(.4, .6),
        ncol = 2, nrow = 2
    )
}
