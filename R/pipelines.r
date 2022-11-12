run_bayes_vs_frequentist <- function(outdir = "output/bayes_vs_frequentist", thr_interval = 0.005, .seed = 123) {
    # Comparison with other packages ----
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 14))
    outdir <- str_path(outdir)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    thresholds <- validate_thresholds(thresholds = thresholds)

    ## GUSTO-I trial ----

    load(str_path("data/gusto.rda"))

    ### Small sample size ----

    set.seed(.seed)
    test_ix <- sample(nrow(gusto), 500)
    dev_data <- gusto[-test_ix, ]
    test_data <- gusto[test_ix, ]


    fit_gusto <- glm(
        day30 ~ age + sysbp + pulse + Killip,
        data = dev_data,
        family = binomial
    )

    test_data$phat1 <- predict(
        fit_gusto,
        newdata = test_data,
        type = "response"
    )

    plot_gusto <- test_data %>%
        plot_bdca_vs_rmda(
            outcomes = "day30",
            predictor = "phat1",
            thresholds = thresholds,
            bootstraps = 2e3
        ) +
        theme(legend.position = c(.7, .7)) +
        labs(title = NULL)

    ggsave(
        str_path("{outdir}/gusto.png"),
        plot_gusto,
        width = 8, height = 4.5, dpi = 600
    )


    ### Reasonable sample size ----

    set.seed(.seed)
    test_ix <- sample(nrow(gusto), 1500)
    dev_data <- gusto[-test_ix, ]
    test_data_large <- gusto[test_ix, ]


    fit_gusto <- glm(
        day30 ~ age + sysbp + pulse + Killip,
        data = dev_data,
        family = binomial
    )

    test_data_large$phat1 <- predict(
        fit_gusto,
        newdata = test_data_large,
        type = "response"
    )

    plot_gusto2 <- test_data_large %>%
        plot_bdca_vs_rmda(
            outcomes = "day30",
            predictor = "phat1",
            thresholds = thresholds,
            bootstraps = 2e3
        ) +
        theme(legend.position = c(.7, .7))

    ggsave(
        str_path("{outdir}/gusto-1500.png"),
        plot_gusto2,
        width = 8, height = 4.5, dpi = 600
    )

    # Compare dca

    misclass <- summary(
        OptimalCutpoints::optimal.cutpoints(
            X = "phat1",
            status = "day30",
            tag.healthy = 0,
            data = test_data_large,
            methods = c("MCT")
        )
    )
    t_misclass <- misclass$MCT$Global$optimal.cutoff$cutoff[1]

    df <- data.frame(
        outcomes = test_data_large$day30,
        model = test_data_large$phat1,
        classifier = as.numeric(test_data_large$phat1 > t_misclass),
        noisy_model = expit(
            qlogis(test_data_large$phat1) + rnorm(nrow(test_data_large), sd = 1.5)
        )
    )

    bdca <- bayesDCA::dca(df, cores = 4, thresholds = thresholds)

    p1 <- compare_dca(bdca, type = "best")
    ggsave(
        str_path("{outdir}/gusto-1500-compare-best.png"), p1,
        width = 10, height = 6, dpi = 600
    )
    p2 <- compare_dca(bdca, type = "useful")
    ggsave(
        str_path("{outdir}/gusto-1500-compare-useful.png"), p2,
        width = 10, height = 6, dpi = 600
    )
    p3 <- compare_dca(bdca, models_or_tests = c("model", "classifier"), type = "pairwise")
    ggsave(
        str_path("{outdir}/gusto-1500-compare-pairwise-model-vs-classifier.png"), p3,
        width = 10, height = 6, dpi = 600
    )

    p4 <- compare_dca(bdca, models_or_tests = c("model", "noisy_model"), type = "pairwise")
    ggsave(
        str_path("{outdir}/gusto-1500-compare-pairwise-model-vs-noisy-model.png"), p4,
        width = 10, height = 6, dpi = 600
    )

    output <- list(
        plot_gusto = plot_gusto,
        plot_gusto2 = plot_gusto2,
        p1 = p1, p2 = p2, p3 = p3, p4 = p4,
        thresholds = thresholds
    )

    return(output)
}

run_simulation_study <- function(thresholds, n_pop, .seed) {
    # Simulation section ----
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 14))
    outdir <- str_path("output/simulation_study")
    thresholds <- validate_thresholds(thresholds = thresholds)

    simulation_settings <- list(
        # AUC 0.65, prev 0.05 vs 0.057
        sim1 = list(
            true_beta = c(-3.1, -log(1.5), log(1.5)),
            beta_hat = c(-3.9, -log(1.5) * 3, log(1.5) * 3)
        ),
        # AUC 0.65, prev 0.3
        sim2 = list(
            true_beta = c(-0.9, -log(1.55), log(1.55)),
            beta_hat = c(-1.2, -log(1.55) * 3, log(1.55) * 3)
        ),
        # AUC 0.85, prev  0.05 TODO
        sim3 = list(
            true_beta = c(-1.3, -log(4.5), log(4.5)),
            beta_hat = c(-2.25, -log(4.5) * 3, log(4.5) * 3)
        ),
        # AUC 0.85, prev  0.3
        sim4 = list(
            true_beta = c(-1.3, -log(4.5), log(4.5)),
            beta_hat = c(-2.25, -log(4.5) * 3, log(4.5) * 3)
        )
    )
    true_beta <- simulation_settings$sim1$true_beta
    beta_hat <- simulation_settings$sim1$beta_hat
    d <- length(true_beta) - 1
    x <- cbind(1, matrix(rexp(n_pop * d, 1), ncol = d))
    true_p <- plogis(as.vector(x %*% true_beta))
    set.seed(122331)
    y <- rbinom(n_pop, 1, true_p)
    p_hat <- plogis(as.vector(x %*% beta_hat))
    simulation <- simulate_dca(
        n_pop = n_pop,
        thresholds = seq(0, 1, .1),
        true_beta = c(-2.1, -log(1.5), log(1.5)),
        beta_hat = c(-3, -log(1.5) * 3, log(1.5) * 3),
        events = 100,
        .seed = .seed
    )
    return(output)
}
