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

run_simulation_study <- function(thresholds, .seed) {
    # Simulation section ----
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 14))
    outdir <- str_path("output/simulation_study")
    thresholds <- validate_thresholds(thresholds = thresholds)
    n <- 2e6
    d <- 2

    ## simulate predictors, true p|x, and y|p
    set.seed(.seed)
    x <- cbind(1, matrix(rexp(n * d, 1), ncol = d))
    true_beta <- c(-2.1, -log(1.5), log(1.5))
    true_p <- plogis(as.vector(x %*% true_beta))
    set.seed(.seed)
    y <- rbinom(n, 1, true_p)
    ## calculate p_hat|x
    w <- 3
    beta_hat <- c(-3, true_beta[2] * w, true_beta[3] * w)
    p_hat <- plogis(as.vector(x %*% beta_hat))
    ## calculate true (approximate) NB|y,p_hat
    true_nb <- map_df(thresholds, ~ {
        tpr <- mean(
            p_hat > .x & y == 1
        )
        fpr <- mean(
            p_hat > .x & y == 0
        )
        tibble(
            .thr = .x,
            nb = tpr - fpr * (.x / (1 - .x))
        )
    })

    sample_size <- ceiling(100 / mean(true_p))
    set.seed(.seed)
    sample_ix <- sample(1:n, sample_size)
    df_sample <- data.frame(
        outcomes = y[sample_ix],
        model_predictions = p_hat[sample_ix]
    )

    .plot <- df_sample %>%
        plot_bdca_vs_rmda(
            outcomes = "outcomes",
            predictor = "model_predictions",
            thresholds = thresholds,
            bootstraps = 2e3
        ) +
        geom_line(
            data = true_nb,
            aes(x = .thr, y = nb, group = 1),
            color = "red", inherit.aes = FALSE
        ) +
        theme(legend.position = c(.7, .7)) +
        ggplot2::scale_y_continuous(expand = c(.275, 0))

    ggsave(
        str_path("{outdir}/simulation.png"),
        .plot,
        width = 7, height = 4, dpi = 600
    )
    output <- list(
        .plot = .plot,
        true_nb = true_nb,
        df_sample = df_sample,
        thresholds = thresholds
    )
    return(output)
}
