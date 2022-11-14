run_bayes_vs_frequentist <- function(outdir = "output/bayes_vs_frequentist", thresholds, .seed = 123) {
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

    plot_gusto <- plot_bdca_vs_rmda(
        dataset = test_data,
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

    plot_gusto2 <- plot_bdca_vs_rmda(
        dataset = test_data_large,
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

run_simulation_study <- function(n_sim, thresholds, n_pop,
                                 outdir, overwrite, .seed) {
    # Simulation section ----
    thresholds <- validate_thresholds(thresholds = thresholds)
    dir.create(
        outdir,
        showWarnings = FALSE,
        recursive = TRUE
    )

    simulation_settings <- get_simulation_settings()
    n_settings <- length(simulation_settings)
    set.seed(.seed)
    settings_seeds <- sample(1:1000, n_settings)
    simulation_results <- vector("list", n_settings * n_sim)
    results_ix <- 1
    for (i in 1:n_settings) {
        .setting <- simulation_settings[[i]]
        .setting_label <- names(simulation_settings)[i]
        .setting_seed <- settings_seeds[i]
        msg <- cli::col_br_magenta(
            paste0("Running simulation setting ", i, " with seed ", .setting_seed)
        )
        message(msg)
        set.seed(.setting_seed)
        sim_seeds <- sample(1:2e4, n_sim)
        for (j in 1:n_sim) {
            .run_seed <- sim_seeds[j]
            msg <- cli::col_br_green(paste0(
                "Run ", j, " with seed ", .run_seed, "\t(setting ", i, ")"
            ))
            message(msg)
            .label <- paste0(
                .setting_label, "_", .setting_seed,
                "-run", j, "_", .run_seed
            )
            .simulation_output <- simulate_dca(
                n_pop = n_pop,
                thresholds = thresholds,
                true_beta = .setting$true_beta,
                beta_hat = .setting$beta_hat,
                events = 100,
                .seed = .run_seed,
                .setting_label = .setting_label,
                .label = .label,
                result_path = str_path("{outdir}/tmp/{.label}.tsv"),
                overwrite = overwrite,
                .verbose = FALSE
            )
            simulation_results[[results_ix]] <- .simulation_output$result
            results_ix <- results_ix + 1
        }
    }

    simulation_results <- as.data.frame(dplyr::bind_rows(simulation_results))

    return(simulation_results)
}

#' Plot simulated DCA results (Results' subsection 02)
#' @param simulation_results Results from `run_simulation_study`.
#' @param outdir Path for output directory
#' @param global_simulation_seed Global seed used for simulation (from `_targets.R` file)
#' @import tidyverse
plot_simulation_results <- function(simulation_results, outdir, global_simulation_seed) {
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 14))
    .colors <- RColorBrewer::brewer.pal(3, "Dark2")
    names(.colors) <- c(
        "True NB", "Frequentist", "Bayesian"
    )
    .colors[".true_nb"] <- .colors[1]

    setting_labels_pretty <- purrr::map_chr(
        get_simulation_settings(),
        ~ paste0(
            "AUC ", .x$auc, ", prevalence ", round(.x$prev * 100), "%"
        )
    )

    df <- simulation_results %>%
        dplyr::filter(strategy == "Model-based decisions") %>%
        dplyr::mutate(
            setting_label = factor(
                setting_labels_pretty[setting_label],
                levels = setting_labels_pretty
            )
        ) %>%
        tibble::as_tibble()

    # point estimates are nearly identical
    p1 <- df %>%
        dplyr::filter(threshold <= 0.9) %>%
        dplyr::select(
            threshold, setting_label, simulation_label, .type, estimate, .true_nb
        ) %>%
        tidyr::pivot_wider(names_from = .type, values_from = estimate) %>%
        tidyr::pivot_longer(cols = c(.true_nb, Frequentist, Bayesian)) %>%
        dplyr::mutate(
            name = ifelse(
                name == ".true_nb",
                "True NB",
                name
            ),
            name = forcats::fct_relevel(
                name,
                "True NB", "Bayesian", "Frequentist"
            )
        ) %>%
        ggplot2::ggplot(ggplot2::aes(factor(threshold), value)) +
        ggplot2::geom_hline(
            yintercept = 0,
            lty = 2,
            alpha = 0.7,
            color = "gray40",
            linewidth = 0.5
        ) +
        ggplot2::geom_boxplot(
            ggplot2::aes(color = name),
            position = ggplot2::position_dodge(width = .8),
            width = .5, lwd = 0.55
        ) +
        ggplot2::facet_wrap(~setting_label, scales = "free_y") +
        ggplot2::scale_color_manual(values = .colors) +
        ggplot2::scale_x_discrete(breaks = scales::pretty_breaks(10)) +
        ggplot2::theme(
            legend.position = c(.1, .15)
        ) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "Net benefit",
            color = NULL
        )

    ggplot2::ggsave(
        str_path("{outdir}/point_estimates_distributions.png"),
        p1,
        width = 12, height = 6.5, dpi = 600
    )

    # 95% intervals coverage

    p2 <- df %>%
        dplyr::group_by(threshold, .type, setting_label) %>%
        dplyr::summarise(
            cov = mean(truth_within_interval),
            .groups = "drop"
        ) %>%
        dplyr::mutate(
            .type = factor(
                .type,
                levels = c("Bayesian", "Frequentist")
            )
        ) %>%
        dplyr::arrange(.type) %>%
        ggplot2::ggplot(ggplot2::aes(factor(threshold), cov, color = .type)) +
        ggplot2::geom_hline(
            yintercept = 0.95,
            lty = 2,
            alpha = 0.7,
            color = "gray40",
            linewidth = 0.5
        ) +
        ggplot2::geom_line(
            ggplot2::aes(linetype = .type, group = .type)
        ) +
        ggplot2::geom_point(
            # position = position_dodge(width = .025),
            ggplot2::aes(shape = .type, size = .type),
            stroke = 1.5
        ) +
        ggplot2::facet_wrap(~setting_label) +
        ggplot2::scale_y_continuous(
            labels = scales::label_percent(),
            limits = c(0.75, 1)
        ) +
        ggplot2::scale_color_manual(values = .colors) +
        ggplot2::scale_shape_manual(values = c(19, 19)) +
        ggplot2::scale_x_discrete(
            breaks = scales::pretty_breaks(10)
        ) +
        ggplot2::scale_size_manual(values = c(2.5, 1.25)) +
        ggplot2::theme(
            legend.position = c(.12, .135)
        ) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "Empirical covarage\n(95% uncertainty intervals)",
            color = NULL,
            size = NULL,
            linetype = NULL,
            shape = NULL
        )

    ggplot2::ggsave(
        str_path("{outdir}/empirical_coverage.png"),
        p2,
        width = 9, height = 5.5, dpi = 600
    )

    return(list(point_estimates = p1, coverage = p2))
}

#' Run Baysian DCA case study (Results' subsection 03)
#' @param .seed Global rng seed.
run_case_study <- function(.seed) {
    import::from(magrittr, `%>%`)

    d <- readr::read_csv(
        str_path("data/12916_2019_1425_MOESM1_ESM.csv"),
        show_col_types = FALSE
    ) %>%
        dplyr::select(
            outcomes := outcome,
            model_predictions := pred
        )

    ## noise predictor
    set.seed(.seed)
    z <- rnorm(nrow(d), mean = 0, sd = 2)

    d <- d %>%
        dplyr::mutate(
            binary_test = as.numeric(model_predictions > 0.5),
            model_predictions2 = plogis(
                qlogis(model_predictions) + log(5) * z
            )
        )


    # Fit Bayesian Decision Curve Analysis ----

    fit <- bayesDCA::dca(d, cores = 4)

    return(fit)
}


plot_case_study_results <- function(fit, outdir) {
    library(patchwork)

    # Plot DCA ----

    .labels <- list(
        "model_predictions" = "ADNEX",
        "binary_test" = "Binary test/Classifier",
        "model_predictions2" = "ADNEX updated"
    )
    dca_plot <- bayesDCA:::plot.BayesDCAList(
        fit,
        labels = .labels
    ) +
        ggplot2::theme(
            legend.position = c(0.2, 0.35),
            legend.text = ggplot2::element_text(size = 9)
        )

    ggplot2::ggsave(
        str_path("{outdir}/dca.png"),
        dca_plot,
        width = 8, height = 4.5, dpi = 600
    )


    # Interrogate posteriors ----
    plots_best <- bayesDCA::compare_dca(
        fit,
        plot_list = TRUE
    )
    plots_useful <- bayesDCA::compare_dca(
        fit,
        plot_list = TRUE,
        type = "useful", labels = .labels
    )

    p1 <- plots_useful$prob_better
    p2 <- plots_useful$delta +
        ggplot2::guides(
            color = ggplot2::guide_legend(title = NULL)
        ) +
        ggplot2::theme(
            legend.position = c(0.75, 0.25),
            legend.text = ggplot2::element_text(size = 8),
            legend.key.size = ggplot2::unit(0.5, "cm")
        )
    p3 <- plots_best$prob_better
    p4 <- plots_best$delta

    posterior_iterrogation_plots <- ((p1 / p2) | (p3 / p4)) +
        patchwork::plot_annotation(tag_levels = "A")

    ggplot2::ggsave(
        str_path("{outdir}/dca-posterior-interrogation.png"),
        posterior_iterrogation_plots,
        width = 10, height = 6.75, dpi = 600
    )

    # EVPI ----

    evpi_plot <- bayesDCA:::plot_evpi(fit)
    ggplot2::ggsave(
        str_path("{outdir}/evpi.png"),
        evpi_plot,
        width = 8, height = 4.5, dpi = 600
    )

    output <- list(
        dca = dca_plot,
        posterior_interrogation = posterior_iterrogation_plots,
        evpi = evpi_plot
    )
    return(output)
}
