run_gusto_trial_example <- function(outdir, thresholds, .seed = 123) {
    # Comparison with other packages ----
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 18))
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
        bootstraps = 2e3,
        show_informative_prior = FALSE
    ) +
        theme(legend.position = c(.7, .7)) +
        labs(title = NULL)

    ggsave(
        str_path("{outdir}/gusto.png"),
        plot_gusto,
        width = 8, height = 4.5, dpi = 600
    )

    output <- list(
        plot_gusto = plot_gusto,
        thresholds = thresholds
    )

    return(output)
}

run_simulation_study <- function(n_sim, thresholds, n_pop,
                                 outdir, overwrite, .seed, .workers = 2, .verbose = FALSE) {
    simulation_results_file <- str_path("{outdir}/simulation_results.tsv")
    if (file.exists(simulation_results_file) && isFALSE(overwrite)) {
        msg <- cli::col_br_red("Simulation results axist and will not be overwritten")
        message(msg)
        simulation_results <- readr::read_tsv(
            simulation_results_file,
            show_col_types = FALSE
        )
        return(simulation_results)
    }
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
    simulation_results <- vector("list", n_settings)
    populations_summaries <- vector("list", n_settings)
    results_ix <- 1
    for (i in 1:n_settings) {
        .setting <- simulation_settings[[i]]
        .setting_label <- names(simulation_settings)[i]
        .setting_seed <- settings_seeds[i]
        msg <- cli::col_br_magenta(
            paste0("Running simulation setting ", i, " with seed ", .setting_seed)
        )
        message(msg)
        # simulate population data
        setting_population <- simulate_dca_population(
            true_beta = .setting$true_beta,
            beta_hat = .setting$beta_hat,
            n_pop = n_pop,
            thresholds = thresholds,
            .seed = .setting_seed,
            .verbose = .verbose
        )
        populations_summaries[[i]] <- list(
            setting_id = i,
            prevalence = mean(setting_population$y),
            average_prediction = mean(setting_population$p_hat),
            calibration_oe = setting_population$calibration_oe,
            calibration_slope = setting_population$calibration_slope,
            auc = setting_population$auc,
            n = length(setting_population$y)
        )
        # simulate samples for DCA
        df_sample_list <- get_setting_sample_list(
            events = 100, # sample size corresponds to expected n events
            n_sim = n_sim, # number of simulated samples
            population_data = setting_population,
            .setting_seed = .setting_seed,
            .setting_label = .setting_label
        )
        .true_prevalence <- mean(setting_population$true_p)
        .true_nb <- setting_population$true_nb
        plan(multisession, workers = .workers)
        run_df <- furrr::future_map_dfr(
            df_sample_list,
            function(df_sample) {
                .run_seed <- df_sample$.run_seed[1]
                .run_label <- df_sample$simulation_run_label[1]
                if (isTRUE(.verbose)) {
                    msg <- cli::col_br_green(paste0(
                        "Run ", df_sample$.run_id[1], " with run seed ", .run_seed, "\t(setting ", i, ")"
                    ))
                    message(msg)
                }
                .simulation_output <- run_dca_simulation(
                    df_sample = df_sample,
                    thresholds = thresholds,
                    true_nb = .true_nb,
                    true_prevalence = .true_prevalence,
                    .setting_label = .setting_label,
                    .run_label = .run_label,
                    result_path = str_path("{outdir}/tmp/{.run_label}.tsv"),
                    overwrite = overwrite,
                    .verbose = .verbose
                )
                return(.simulation_output$result)
            },
            .options = furrr::furrr_options(seed = .setting_seed)
        )
        plan(sequential)
        simulation_results[[i]] <- run_df
    }

    simulation_results <- as.data.frame(dplyr::bind_rows(simulation_results))
    readr::write_tsv(
        simulation_results,
        simulation_results_file
    )
    readr::write_tsv(
        dplyr::bind_rows(populations_summaries),
        str_path("{outdir}/populations_summaries.tsv")
    )

    return(simulation_results)
}

#' Plot simulated DCA results (Results' subsection 02)
#' @param simulation_results Results from `run_simulation_study`.
#' @param outdir Path for output directory
#' @param global_simulation_seed Global seed used for simulation (from `_targets.R` file)
#' @import tidyverse
plot_simulation_results <- function(simulation_results, outdir, global_simulation_seed,
                                    surv = FALSE, estimation_types = NULL, .colors = NULL) {
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 18))
    dir.create(
        outdir,
        showWarnings = FALSE,
        recursive = TRUE
    )
    if (is.null(estimation_types)) {
        estimation_types <- c(
            "Bayesian", "Bayesian 2", "Frequentist"
        )
    }

    if (is.null(.colors)) {
        .colors <- c(
            "True NB" = "#1B9E77",
            "Bayesian" = "#7570B3",
            "Bayesian 2" = "red",
            "Frequentist" = "#D95F02"
        )
    }

    stopifnot(all(estimation_types %in% names(.colors)))

    .colors <- .colors[names(.colors) %in% c(estimation_types, "True NB")]
    n_types <- length(estimation_types)
    simulation_results <- simulation_results[
        simulation_results$.type %in% estimation_types,
    ]

    # make sure NA intervals get NA truth_within_interval
    simulation_results <- simulation_results %>%
        dplyr::mutate(
            truth_within_interval = ifelse(
                is.na(estimate) | is.na(.lower) | is.na(.upper),
                NA,
                truth_within_interval
            )
        )

    if (isFALSE(surv)) {
        setting_labels_pretty <- purrr::map_chr(
            get_simulation_settings(),
            ~ paste0(
                "AUC ", .x$auc, ", prevalence ", round(.x$prev * 100), "%"
            )
        )
        plot_height <- 6.5
    } else {
        setting_labels_pretty <- purrr::map_chr(
            get_simulation_settings_surv(),
            ~ paste0(
                "C-statistic ", .x$concord,
                ", 1-year survival ",
                round(.x$one_year_survival_rate * 100),
                "%"
            )
        )
        plot_height <- 8.5
    }

    df <- simulation_results %>%
        dplyr::filter(strategy == "Model-based decisions") %>%
        dplyr::mutate(
            setting_label = factor(
                setting_labels_pretty[setting_label],
                levels = setting_labels_pretty
            )
        ) %>%
        tibble::as_tibble()

    thresholds <- sort(unique(df$threshold))

    # point estimates are nearly identical
    p1 <- df %>%
        dplyr::filter(threshold <= .75) %>%
        dplyr::select(
            threshold, setting_label, simulation_run_label, .type, estimate, .true_nb
        ) %>%
        tidyr::pivot_wider(names_from = .type, values_from = estimate) %>%
        tidyr::pivot_longer(cols = dplyr::any_of(c(".true_nb", estimation_types))) %>%
        dplyr::mutate(
            name = ifelse(
                name == ".true_nb",
                "True NB",
                name
            ),
            name = factor(as.character(name), levels = names(.colors)),
            thr_factor = factor(
                threshold,
                levels = thresholds,
                labels = get_thr_label(thresholds)
            )
        ) %>%
        dplyr::filter(!is.na(value)) %>%
        ggplot2::ggplot(ggplot2::aes(thr_factor, value)) +
        ggplot2::geom_hline(
            yintercept = 0,
            lty = 2,
            alpha = 0.7,
            color = "gray40",
            linewidth = 0.5
        ) +
        ggplot2::geom_boxplot(
            data = . %>%
                dplyr::filter(
                    stringr::str_detect(
                        stringr::str_to_lower(name),
                        "true nb",
                        negate = TRUE
                    )
                ),
            ggplot2::aes(color = name),
            position = ggplot2::position_dodge(width = .75),
            width = .5, lwd = 0.55,
            show.legend = FALSE
        ) +
        ggplot2::geom_segment(
            data = . %>%
                dplyr::filter(
                    stringr::str_detect(
                        stringr::str_to_lower(name),
                        "true nb",
                        negate = FALSE
                    )
                ),
            ggplot2::aes(
                y = value, yend = value,
                x = as.numeric(thr_factor) - 0.5,
                xend = as.numeric(thr_factor) + 0.5,
                color = name
            ),
            linewidth = 0.9
        ) +
        ggplot2::facet_wrap(~setting_label, scales = "free_y") +
        ggplot2::scale_color_manual(values = .colors) +
        ggplot2::theme(
            legend.position = "top"
        ) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "Net benefit",
            color = NULL
        )

    ggplot2::ggsave(
        str_path("{outdir}/point_estimates_distributions.png"),
        p1,
        width = 15, height = plot_height, dpi = 600
    )

    # 95% intervals coverage
    n_types <- df %>%
        dplyr::filter(threshold <= .75, !is.na(truth_within_interval)) %>%
        dplyr::pull(.type) %>%
        dplyr::n_distinct()

    if (n_types == 1) {
        point_sizes <- 1
    } else if (n_types == 2) {
        point_sizes <- c(1, 2 / 3)
    } else {
        point_sizes <- n_types:1
    }

    p2 <- df %>%
        dplyr::filter(threshold <= .75, !is.na(truth_within_interval)) %>%
        dplyr::group_by(threshold, .type, setting_label) %>%
        dplyr::summarise(
            cov = list(binom::binom.wilson(sum(truth_within_interval), n())),
            .groups = "drop"
        ) %>%
        tidyr::unnest_wider(cov) %>%
        dplyr::rename(
            cov := mean,
            cov_lower := lower,
            cov_upper := upper
        ) %>%
        dplyr::mutate(
            .type = factor(
                .type,
                levels = sort(estimation_types)
            ),
            thr_factor = factor(
                threshold,
                levels = thresholds,
                labels = get_thr_label(thresholds)
            )
        ) %>%
        dplyr::arrange(.type) %>%
        ggplot2::ggplot(ggplot2::aes(thr_factor, cov, color = .type)) +
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
        ggplot2::geom_pointrange(
            ggplot2::aes(
                ymin = cov_lower, ymax = cov_upper,
                shape = .type, size = .type
            ),
            stroke = 1.5
        ) +
        ggplot2::facet_wrap(~setting_label) +
        ggplot2::scale_y_continuous(
            labels = scales::label_percent()
        ) +
        ggplot2::scale_color_manual(values = .colors) +
        ggplot2::scale_shape_manual(values = rep(19, n_types)) +
        ggplot2::scale_size_manual(values = point_sizes) +
        ggplot2::theme(
            legend.position = c(.12, .135)
        ) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "Empirical coverage\n(95% uncertainty intervals)",
            color = NULL,
            size = NULL,
            linetype = NULL,
            shape = NULL
        )

    ggplot2::ggsave(
        str_path("{outdir}/empirical_coverage.png"),
        p2,
        width = 14, height = plot_height, dpi = 600
    )

    # failure rate
    p2.2 <- df %>%
        dplyr::filter(threshold <= .75) %>%
        dplyr::group_by(threshold, .type, setting_label) %>%
        dplyr::summarise(
            failure = list(binom::binom.wilson(sum(is.na(truth_within_interval)), n())),
            .groups = "drop"
        ) %>%
        tidyr::unnest_wider(failure) %>%
        dplyr::rename(
            failure := mean,
            failure_lower := lower,
            failure_upper := upper
        ) %>%
        dplyr::mutate(
            .type = factor(
                .type,
                levels = sort(estimation_types)
            ),
            thr_factor = factor(
                threshold,
                levels = thresholds,
                labels = get_thr_label(thresholds)
            )
        ) %>%
        dplyr::arrange(.type) %>%
        ggplot2::ggplot(ggplot2::aes(thr_factor, failure, color = .type)) +
        ggplot2::geom_line(
            ggplot2::aes(linetype = .type, group = .type)
        ) +
        ggplot2::geom_pointrange(
            ggplot2::aes(
                ymin = failure_lower, ymax = failure_upper,
                shape = .type, size = .type
            ),
            stroke = 1.5
        ) +
        ggplot2::facet_wrap(~setting_label) +
        ggplot2::scale_y_continuous(
            labels = scales::label_percent()
        ) +
        ggplot2::scale_color_manual(values = .colors) +
        ggplot2::scale_shape_manual(values = rep(19, n_types)) +
        ggplot2::scale_size_manual(values = point_sizes) +
        ggplot2::theme(
            legend.position = c(0.08, .85)
        ) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "Failed estimation frequency",
            color = NULL,
            size = NULL,
            linetype = NULL,
            shape = NULL
        )

    ggplot2::ggsave(
        str_path("{outdir}/empirical_failure_rate.png"),
        p2.2,
        width = 14, height = plot_height, dpi = 600
    )

    # raw error
    p3 <- df %>%
        dplyr::filter(threshold <= .75) %>%
        dplyr::mutate(abs_error = estimate - .true_nb) %>%
        dplyr::select(
            threshold, setting_label, simulation_run_label, .type, abs_error
        ) %>%
        tidyr::pivot_wider(names_from = .type, values_from = abs_error) %>%
        tidyr::pivot_longer(cols = dplyr::any_of(estimation_types)) %>%
        dplyr::mutate(
            name = factor(
                as.character(name),
                levels = names(.colors)
            ),
            thr_factor = factor(
                threshold,
                levels = thresholds,
                labels = get_thr_label(thresholds)
            )
        ) %>%
        dplyr::filter(!is.na(value)) %>%
        ggplot2::ggplot(ggplot2::aes(thr_factor, value)) +
        ggplot2::geom_hline(
            yintercept = 0,
            lty = 2,
            alpha = 0.7,
            color = "gray40",
            linewidth = 0.5
        ) +
        ggplot2::geom_boxplot(
            data = . %>%
                dplyr::filter(
                    stringr::str_detect(
                        stringr::str_to_lower(name),
                        "true nb",
                        negate = TRUE
                    )
                ),
            ggplot2::aes(color = name),
            position = ggplot2::position_dodge(width = .75),
            width = .5, lwd = 0.55
        ) +
        ggplot2::facet_wrap(~setting_label, scales = "free_y") +
        ggplot2::scale_color_manual(values = .colors) +
        ggplot2::theme(
            legend.position = "top"
        ) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "Estimated NB - True NB",
            color = NULL
        )

    if (isTRUE(surv)) {
        p3 <- p3 +
            ggplot2::scale_y_continuous(breaks = c(-.5, 0, .5)) +
            ggplot2::coord_cartesian(ylim = c(-.7, .7))
    }

    ggplot2::ggsave(
        str_path("{outdir}/point_estimates_error.png"),
        p3,
        width = 15, height = plot_height, dpi = 600
    )

    # MAPE -- scale absolute errors by maximum achievable NB
    if (isTRUE(surv)) {
        df$ape <- abs(df$abs_error) / df$.true_nb
        .ylim <- c(0, NA)
        .title <- "Mean Absolute Percentage Error (MAPE)"
    } else {
        df$ape <- abs(df$abs_error) / df$true_prevalence
        .ylim <- c(0, .2)
        .title <- "Average Percentage Error relative to max NB"
    }
    # plot mape
    p4 <- df %>%
        dplyr::filter(threshold <= .75) %>%
        dplyr::select(
            threshold, setting_label, simulation_run_label, .type, ape
        ) %>%
        tidyr::pivot_wider(names_from = .type, values_from = ape) %>%
        tidyr::pivot_longer(cols = dplyr::any_of(estimation_types)) %>%
        dplyr::mutate(
            name = factor(
                as.character(name),
                levels = names(.colors)
            ),
            thr_factor = factor(
                threshold,
                levels = thresholds,
                labels = get_thr_label(thresholds)
            )
        ) %>%
        dplyr::filter(!is.na(value)) %>%
        dplyr::group_by(thr_factor, setting_label, name) %>%
        dplyr::summarise(
            mape = ggpubr::mean_ci(value),
            .groups = "drop"
        ) %>%
        tidyr::unnest_wider(mape) %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = thr_factor,
                y = y,
                ymin = ymin,
                ymax = ymax,
                color = name,
                fill = name
            )
        ) +
        ggplot2::geom_col(
            position = ggplot2::position_dodge(width = 0.6),
            width = 1e-3,
            lty = "dotted",
            alpha = 0.2,
            show.legend = FALSE
        ) +
        ggplot2::geom_pointrange(
            position = ggplot2::position_dodge(width = 0.6),
            size = 1, linewidth = 1.2,
            pch = 21, color = "gray20"
        ) +
        ggplot2::facet_wrap(~setting_label, scales = "free_y") +
        ggplot2::scale_color_manual(values = .colors) +
        ggplot2::scale_fill_manual(values = .colors) +
        ggplot2::scale_y_continuous(labels = scales::percent, limits = .ylim) +
        ggplot2::theme(
            legend.position = c(0.075, 0.9)
        ) +
        ggplot2::labs(
            x = "Decision threshold",
            title = .title,
            y = NULL,
            color = NULL,
            fill = NULL
        )

    ggplot2::ggsave(
        str_path("{outdir}/average_error.png"),
        p4,
        width = 15, height = plot_height, dpi = 600
    )

    # ci width
    p5 <- df %>%
        dplyr::filter(threshold <= .75, !is.na(truth_within_interval)) %>%
        dplyr::mutate(ci_width = .upper - .lower) %>%
        dplyr::select(
            threshold, setting_label, simulation_run_label, .type, ci_width
        ) %>%
        tidyr::pivot_wider(names_from = .type, values_from = ci_width) %>%
        tidyr::pivot_longer(cols = dplyr::any_of(estimation_types)) %>%
        dplyr::mutate(
            name = factor(
                as.character(name),
                levels = names(.colors)
            ),
            thr_factor = factor(
                threshold,
                levels = thresholds,
                labels = get_thr_label(thresholds)
            )
        ) %>%
        ggplot2::ggplot(ggplot2::aes(thr_factor, value)) +
        ggplot2::geom_hline(
            yintercept = 0,
            lty = 2,
            alpha = 0.7,
            color = "gray40",
            linewidth = 0.5
        ) +
        ggplot2::geom_boxplot(
            data = . %>%
                dplyr::filter(
                    stringr::str_detect(
                        stringr::str_to_lower(name),
                        "true nb",
                        negate = TRUE
                    )
                ),
            ggplot2::aes(color = name),
            position = ggplot2::position_dodge(width = .75),
            width = .5, lwd = 0.55
        ) +
        ggplot2::facet_wrap(~setting_label, scales = "free_y") +
        ggplot2::scale_color_manual(values = .colors) +
        ggplot2::theme(
            legend.position = "top"
        ) +
        ggplot2::labs(
            x = "Decision threshold",
            y = "95% Uncertainty interval width",
            color = NULL
        )

    ggplot2::ggsave(
        str_path("{outdir}/ci_width.png"),
        p5,
        width = 15, height = plot_height, dpi = 600
    )

    # runtime
    p6 <- df %>%
        dplyr::filter(threshold <= .75, !is.na(estimate)) %>%
        dplyr::select(.type, setting_label, runtime, simulation_run_label) %>%
        unique() %>%
        ggplot2::ggplot(ggplot2::aes(.type, runtime)) +
        ggplot2::geom_boxplot(alpha = 0) +
        ggbeeswarm::geom_quasirandom() +
        ggplot2::facet_wrap(~setting_label) +
        ggplot2::labs(x = NULL, y = "Runtime (seconds)")

    ggplot2::ggsave(
        str_path("{outdir}/runtime.png"),
        p6,
        width = 15, height = plot_height, dpi = 600
    )

    return(list(point_estimates = p1, coverage = p2, abs_error = p3, mape = p4, ci_width = p5, runtime = p6))
}

#' Run Baysian DCA case study (Results' subsection 03)
#' @param .seed Global rng seed.
run_case_study <- function(thresholds, .seed) {
    import::from(magrittr, `%>%`)

    d <- readr::read_csv(
        str_path("data/12916_2019_1425_MOESM1_ESM.csv"),
        show_col_types = FALSE
    ) %>%
        dplyr::select(
            outcomes := outcome,
            model_predictions := pred
        )

    soc_test_se <- 0.81
    soc_test_sp <- 0.88

    .seed <- 12112022 # for reproducibility of example percentages in the paper
    set.seed(.seed)
    d <- d %>%
        rowwise() %>%
        dplyr::mutate(
            binary_test = ifelse(
                outcomes == 1,
                rbinom(n = 1, size = 1, soc_test_se),
                rbinom(n = 1, size = 1, 1 - soc_test_sp)
            )
        )


    # Fit Bayesian Decision Curve Analysis ----

    fit <- bayesDCA::dca(d, thresholds = thresholds)

    return(fit)
}


plot_case_study_results <- function(fit, outdir) {
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 18))
    dir.create(
        outdir,
        showWarnings = FALSE,
        recursive = TRUE
    )
    library(patchwork)
    library(ggdist)

    # Plot DCA ----

    .labels <- list(
        "model_predictions" = "ADNEX",
        "binary_test" = "Hypothetical Standard of Care test"
    )
    dca_plot <- bayesDCA:::plot.BayesDCA(
        fit,
        labels = .labels,
        linewidth = 1.25
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
    plots_useful <- bayesDCA::compare_dca(
        fit,
        plot_list = TRUE,
        type = "useful", labels = .labels,
        linewidth = 1.25
    )
    plots_best <- bayesDCA::compare_dca(
        fit,
        plot_list = TRUE,
        linewidth = 1.25
    )

    # tweak plots
    p1 <- plots_useful$prob_better +
        ggplot2::guides(
            color = ggplot2::guide_legend(title = NULL)
        ) +
        ggplot2::theme(
            legend.position = c(0.65, 0.25),
            legend.text = ggplot2::element_text(size = 8),
            legend.key.size = ggplot2::unit(0.5, "cm")
        ) +
        ggplot2::labs(tag = "A")

    p2 <- plots_useful$delta
    suppressMessages({
        suppressWarnings({
            p2_zoom <- p2 +
                ggplot2::coord_cartesian(
                    ylim = c(-.001, .025),
                    xlim = c(.04, .08)
                ) +
                ggplot2::scale_y_continuous(
                    breaks = c(0, 0.01, 0.02),
                    labels = scales::label_number(0.01)
                ) +
                ggplot2::scale_x_continuous(
                    breaks = seq(0, .1, .01),
                    labels = scales::label_percent(1)
                ) +
                ggplot2::guides(color = "none", fill = "none") +
                ggplot2::theme_bw(base_size = 7) +
                ggplot2::labs(subtitle = NULL) +
                ggplot2::theme(
                    plot.background = ggplot2::element_rect(
                        colour = "black", fill = NA, linewidth = 1
                    ),
                    axis.title.x = ggplot2::element_text(size = 5)
                )
        })
    })
    p2 <- p2 +
        ggplot2::labs(tag = "B") +
        patchwork::inset_element(
            p2_zoom,
            left = .575,
            bottom = .02,
            top = .45,
            right = .98,
            ignore_tag = TRUE
        )

    p3 <- plots_best$prob_better +
        ggplot2::labs(tag = "C")
    p4 <- plots_best$delta +
        ggplot2::labs(tag = "D") +
        ggplot2::coord_cartesian(ylim = c(-0.085, 0.085))

    posterior_iterrogation_plots <- (
        (p1 | p2) / (p3 | p4)
    ) & ggplot2::theme(
        plot.tag = ggplot2::element_text(face = "bold")
    )

    ggplot2::ggsave(
        str_path("{outdir}/dca-posterior-interrogation.png"),
        posterior_iterrogation_plots,
        width = 10, height = 6.75, dpi = 600
    )

    # Pairwise ----
    strategies <- c(
        "model_predictions", "binary_test"
    )
    suppressMessages({
        suppressWarnings({
            pairwise_delta_plot <- bayesDCA::plot_delta(
                fit,
                type = "pairwise",
                strategies = strategies,
                labels = .labels[1:2],
                linewidth = 1.25
            ) +
                ggplot2::coord_cartesian(ylim = c(-.1, .1)) +
                ggplot2::labs(subtitle = "Net benefit differences across all thresholds", tag = "A")
        })
    })
    nb1 <- fit$fit$distributions$net_benefit[[strategies[1]]]
    nb2 <- fit$fit$distributions$net_benefit[[strategies[2]]]
    .delta <- nb1 - nb2
    ix_41 <- dplyr::near(fit$thresholds, 0.41)
    .delta_41 <- .delta[, ix_41]

    prob_superior_any <- round(mean(.delta_41 > 0.0) * 100)
    prob_superior_mcid <- round(mean(.delta_41 > 0.01) * 100)
    suppressMessages({
        suppressWarnings({
            pairwise_41_plot <- data.frame(
                x = .delta_41
            ) %>%
                ggplot2::ggplot(ggplot2::aes(x = x)) +
                ggdist::stat_halfeye(
                    slab_color = "gray45",
                    ggplot2::aes(
                        fill = ggplot2::after_stat(x > 0)
                    ),
                    interval_size_range = c(1.5, 3),
                    normalize = "none", adjust = 2
                ) +
                ggplot2::scale_y_continuous(limits = c(0, 50)) +
                ggplot2::labs(
                    x = latex2exp::TeX("$\\Delta_{NB}$"),
                    y = "Posterior density", tag = "B",
                    subtitle = "Difference at 41% threshold"
                ) +
                ggplot2::guides(fill = "none") +
                ggplot2::scale_fill_manual(
                    values = c("#bdd5ed", "#4788c4")
                ) +
                ggplot2::geom_vline(
                    xintercept = 0,
                    lty = 2,
                    color = "gray30"
                ) +
                ggplot2::geom_segment(
                    x = 1e-3, y = 45,
                    xend = 0.035, yend = 45,
                    lineend = "round",
                    linejoin = "round",
                    linewidth = 1.5,
                    arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "inches")),
                    colour = "gray30"
                ) +
                ggplot2::geom_segment(
                    x = -1e-3, y = 45,
                    xend = -0.035, yend = 45,
                    lineend = "round",
                    linejoin = "round",
                    linewidth = 1.5,
                    arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "inches")),
                    colour = "gray30" # Also accepts "red", "blue' etc
                ) +
                ggplot2::annotate(
                    "text",
                    x = .035, y = 32.5, size = 4, color = "#4788c4",
                    label = latex2exp::TeX(paste0("$P(\\Delta_{NB} > 0)$ = ", prob_superior_any, "%"))
                ) +
                ggplot2::annotate(
                    "text",
                    hjust = 0, fontface = "bold",
                    x = 0.0025, y = 48, size = 4, color = "gray30",
                    label = "Favors ADNEX"
                ) +
                ggplot2::annotate(
                    "text",
                    hjust = 1, fontface = "bold",
                    x = -0.0025, y = 48, size = 4, color = "gray30",
                    label = "Favors SoC"
                ) +
                ggplot2::coord_cartesian(xlim = c(-0.05, 0.05))
        })
    })

    pairwise_plot <- (pairwise_delta_plot | pairwise_41_plot) +
        patchwork::plot_annotation(
            title = "ADNEX versus Hypothetical Standard of Care test",
            theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
        )
    ggplot2::ggsave(
        str_path("{outdir}/pairwise-41.png"),
        pairwise_plot,
        width = 12, height = 4.5, dpi = 600
    )

    # EVPI ----

    evpi_plot <- bayesDCA:::plot_evpi(fit, linewidth = 1.25) +
        ggplot2::labs(subtitle = "Expected Value of Perfect Information (EVPI)")
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

run_simulation_study_surv <- function(n_sim, thresholds, n_pop,
                                      outdir, overwrite, .seed,
                                      pred_time = 12, .workers = 2, .verbose = FALSE) {
    simulation_results_file <- str_path("{outdir}/simulation_results_surv.tsv")
    if (file.exists(simulation_results_file) && isFALSE(overwrite)) {
        msg <- cli::col_br_red("Simulation results (survival) exist and will not be overwritten")
        message(msg)
        simulation_results <- readr::read_tsv(
            simulation_results_file,
            show_col_types = FALSE
        )
        return(simulation_results)
    }
    # Simulation section ----
    thresholds <- validate_thresholds(thresholds = thresholds)
    dir.create(
        outdir,
        showWarnings = FALSE,
        recursive = TRUE
    )

    simulation_settings <- get_simulation_settings_surv()
    n_settings <- length(simulation_settings)
    set.seed(.seed)
    settings_seeds <- sample(1:1000, n_settings)
    simulation_results <- vector("list", n_settings)
    populations_summaries <- vector("list", n_settings)
    for (i in 1:n_settings) {
        .setting <- simulation_settings[[i]]
        .setting_label <- names(simulation_settings)[i]
        .setting_seed <- settings_seeds[i]
        msg <- cli::col_br_magenta(
            paste0("Running survival simulation setting ", i, " with seed ", .setting_seed)
        )
        message(msg)
        # simulate population data
        setting_population <- simulate_dca_population_surv(
            sim_setting = .setting,
            n_pop = n_pop,
            thresholds = thresholds,
            .seed = .setting_seed,
            pred_time = pred_time,
            .verbose = .verbose
        )
        populations_summaries[[i]] <- list(
            setting_id = i,
            event_rate = setting_population$event_rate,
            average_prediction = setting_population$average_prediction,
            calibration_oe = setting_population$calibration_oe,
            calibration_slope = setting_population$calibration_slope,
            concordance = setting_population$concordance,
            n = n_pop
        )
        # simulate samples for DCA
        df_sample_list <- get_setting_sample_list_surv(
            events = 100, # sample size corresponds to expected 100 events
            n_sim = n_sim, # number of simulated samples
            population_data = setting_population,
            .setting_seed = .setting_seed,
            .setting_label = .setting_label
        )
        .true_incidence <- mean(setting_population$df_pop$survTime <= pred_time) # actually death rate
        plan(multisession, workers = .workers)
        run_df <- furrr::future_map_dfr(
            df_sample_list,
            function(df_sample) {
                .run_seed <- df_sample$.run_seed[1]
                .run_label <- df_sample$simulation_run_label[1]
                if (isTRUE(.verbose)) {
                    msg <- cli::col_br_green(paste0(
                        "Run ", df_sample$.run_id[1], " with run seed ", .run_seed, "\t(setting ", i, ")"
                    ))
                    message(msg)
                }
                .simulation_output <- run_dca_simulation_surv(
                    df_sample = df_sample,
                    thresholds = thresholds,
                    pred_time = pred_time,
                    true_nb = setting_population$true_nb,
                    true_incidence = .true_incidence,
                    .setting_label = .setting_label,
                    .run_label = .run_label,
                    result_path = str_path("{outdir}/tmp/{.run_label}.tsv"),
                    overwrite = overwrite,
                    .verbose = .verbose
                )
                return(.simulation_output$result)
            },
            .options = furrr::furrr_options(seed = .setting_seed)
        )
        plan(sequential)
        simulation_results[[i]] <- run_df
    }

    simulation_results <- as.data.frame(dplyr::bind_rows(simulation_results))
    readr::write_tsv(
        simulation_results,
        simulation_results_file
    )
    readr::write_tsv(
        dplyr::bind_rows(populations_summaries),
        str_path("{outdir}/populations_summaries.tsv")
    )

    return(simulation_results)
}

merge_survival_simulation_plots <- function(survival_simulation_plots, outdir) { # nolint
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 18))
    outdir <- str_path(outdir)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    library(patchwork)
    p1 <- survival_simulation_plots$abs_error +
        ggplot2::labs(x = NULL) +
        ggplot2::facet_wrap(~setting_label) +
        ggplot2::theme(legend.position = c(0.12, 0.85))
    p1$data <- p1$data %>%
        dplyr::filter(stringr::str_detect(setting_label, "C-statistic 0.95", negate = TRUE)) %>%
        dplyr::filter(stringr::str_detect(setting_label, "1-year survival 10%"))
    p2 <- survival_simulation_plots$coverage +
        ggplot2::theme(
            legend.position = c(.12, .3)
        )
    p2$data <- p2$data %>%
        dplyr::filter(stringr::str_detect(setting_label, "C-statistic 0.95", negate = TRUE)) %>%
        dplyr::filter(stringr::str_detect(setting_label, "1-year survival 10%"))

    final_fig <- p1 / p2
    ggplot2::ggsave(
        str_path("{outdir}/survival-simulation-final-figure.png"),
        final_fig,
        width = 11, height = 9, dpi = 600
    )
    return(final_fig)
}

plot_informative_priors_ppc <- function(thresholds, outdir) {
    # not that informative
    fit_informative <- bayesDCA::dca(
        data.frame(outcomes = 1, model = 1),
        thresholds = thresholds,
        threshold_varying_prior = TRUE,
        prior_only = TRUE
    )
    ppc <- bayesDCA::plot_ppc(fit_informative) +
        patchwork::plot_annotation(tag_levels = "A") &
        ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold"))
    ggsave(
        str_path("{outdir}/ppc.png"),
        ppc,
        width = 14, height = 7.5, dpi = 600
    )
    # informative
    fit_informative2 <- bayesDCA::dca(
        data.frame(outcomes = 1, model = 1),
        thresholds = thresholds,
        threshold_varying_prior = TRUE,
        ignorance_region_cutpoints = NULL,
        prior_only = TRUE,
        prev_prior_mean = 0.3,
        prev_prior_sample_size = 50,
        max_sens_prior_sample_size = 5
    )
    ppc2 <- bayesDCA::plot_ppc(fit_informative2) +
        patchwork::plot_annotation(tag_levels = "A") &
        ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold"))
    ggsave(
        str_path("{outdir}/ppc2.png"),
        ppc2,
        width = 14, height = 7.5, dpi = 600
    )
}


#' Run EVPI simulations
#'
#' Code adapted from Sadatsafavi et al (2023) [DOI: 10.1177/0272989X231178317]
#' by Giuliano Netto Flores Cruz on June 27, 2023.
#' Original source code URL:
#' https://github.com/resplab/papercode/blob/main/voipred_ex/CaseStudy.Rmd
run_evpi_simulation <- function(
    n_sim,
    outdir) {
    import::from(magrittr, `%>%`)

    data(gusto, package = "predtools")
    gusto$kill <- (as.numeric(gusto$Killip) > 1) * 1
    gusto$Y <- gusto$day30
    data_us <- gusto[gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15), ]
    data_other <- gusto[!gusto$regl %in% c(1, 7, 9, 10, 11, 12, 14, 15), ]
    dev_data <- data_other

    model <- glm(
        Y ~ age + miloc + pmi + kill + pmin(sysbp, 100) + pulse,
        data = dev_data,
        family = binomial(link = "logit")
    )

    res <- sim_by_size(
        data_us = data_us,
        model = model,
        n_sim = n_sim,
        sample_sizes = c(250, 500, 1000, 2000, 4000, 8000, 16000, Inf)
    )
    pres <- res %>%
        dplyr::group_by(method, sample_size) %>%
        dplyr::summarise(
            val1 = mean(val1),
            val2 = mean(val2),
            val3 = mean(val3),
            val4 = mean(val4),
            .groups = "drop"
        ) %>%
        as.data.frame()
    .levels <- c(
        "Asymptotic",
        "Ordinary bootstrap",
        "Bayesian bootstrap",
        "BayesDCA (uniform)",
        "BayesDCA (informative)"
    )
    .colors <- c(
        "orange",
        "red",
        "blue",
        "gray40",
        "gray10"
    ) %>%
        setNames(.levels)

    ggplot2::theme_set(ggplot2::theme_bw(base_size = 16))
    .plot <- pres %>%
        tidyr::pivot_longer(cols = contains("val")) %>%
        dplyr::filter(
            (name == "val1") |
                (name == "val2" & sample_size <= 8000) |
                (name == "val3" & sample_size <= 4000) |
                (name == "val4" & sample_size <= 4000)
        ) %>%
        dplyr::mutate(
            thr = c("val1" = 0.01, "val2" = 0.02, "val3" = 0.05, "val4" = 0.1)[name],
            thr = paste0("threshold = ", thr),
            method = dplyr::case_match(
                method,
                "asy" ~ .levels[1],
                "BB" ~ .levels[2],
                "OB" ~ .levels[3],
                "full_bayes" ~ .levels[4],
                "full_bayes_informative" ~ .levels[5]
            ),
            method = factor(
                method,
                levels = .levels
            )
        ) %>%
        ggplot2::ggplot(
            ggplot2::aes(sample_size, value, color = method)
        ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::facet_wrap(
            ~thr,
            scales = "free"
        ) +
        ggplot2::scale_color_manual(
            values = .colors,
            name = NULL
        ) +
        ggplot2::labs(
            x = "Sample size",
            y = "EVPI"
        ) +
        ggplot2::theme(
            legend.position = "top"
        ) +
        ggplot2::scale_y_continuous(
            labels = (
                function(x) {
                    ifelse(
                        x == 0,
                        0,
                        scales::scientific_format()(x)
                    )
                }),
            limits = c(0, NA)
        )

    ggplot2::ggsave(
        str_path("{outdir}/evpi-simulation.png"),
        .plot,
        width = 14,
        height = 7.5,
        dpi = 600
    )
}


run_bayes_risk_simulation <- function(n_sim, outdir, overwrite = FALSE) {
    simulation_results_file <- str_path("{outdir}/simulation_results_bayes_risk.tsv")
    if (file.exists(simulation_results_file) && isFALSE(overwrite)) {
        msg <- cli::col_br_red("Simulation results (bayes risk) exist and will not be overwritten")
        message(msg)
        bayes_risk_simulation_results <- readr::read_tsv(
            simulation_results_file,
            show_col_types = FALSE
        )
    } else {
        msg <- cli::col_br_red(str_glue("Starting bayes risk simulation with n_sim={n_sim}"))
        message(msg)
        # ensure outdir exists
        dir.create(
            outdir,
            showWarnings = FALSE,
            recursive = TRUE
        )
        # each sim setting is a sample size
        expected_events <- c(50, 100, 200)
        sample_size_list <- expected_events / 0.2
        message(cli::col_br_red(
            paste0(
                "sample_size_list: ",
                paste0(sample_size_list, collapse = ", ")
            )
        ))
        bayes_risk_simulation_results <- vector("list", length(sample_size_list))
        # run simulation for each sample size
        for (i in seq_along(sample_size_list)) {
            msg <- cli::col_br_magenta(
                paste0("Running bayes risk simulation with sample_size=", sample_size_list[i])
            )
            message(msg)
            bayes_risk_simulation_results[[i]] <- purrr::map_df(
                seq_len(n_sim),
                ~ run_bayes_risk_simulation_single_run(n = sample_size_list[[i]]),
                .id = "run_id"
            )
        }
        # combine and save results
        bayes_risk_simulation_results <- dplyr::bind_rows(
            bayes_risk_simulation_results,
            .id = "setting_id"
        ) %>%
            dplyr::select(
                threshold, decision_strategy, sample_size,
                setting_id, run_id,
                prob, delta
            )
        readr::write_tsv(
            bayes_risk_simulation_results,
            simulation_results_file
        )
    }

    # compute grouond truth
    msg <- cli::col_br_magenta("Computing ground truth")
    message(msg)
    dca_ground_truth <- get_sample_size_dca(n = 2e6)
    dca_ground_truth_plot <- plot(
        dca_ground_truth,
        strategies = c("phat0", "phat1"),
        labels = c("phat0" = "model A", "phat1" = "model B")
    ) +
        ggplot2::scale_x_continuous(
            breaks = seq(0, 0.3, by = 0.05),
            labels = scales::label_percent(1)
        ) +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::theme(
            legend.position = "inside",
            legend.position.inside = c(0.85, 0.8),
            legend.text = ggplot2::element_text(size = 11)
        ) +
        ggplot2::labs(
            subtitle = "True decision curves"
        )

    ggplot2::ggsave(
        str_path("{outdir}/ground_truth.png"),
        dca_ground_truth_plot,
        width = 9,
        height = 5.5,
        dpi = 600
    )

    dca_ground_truth_delta_plot <- bayesDCA:::plot_delta(
        dca_ground_truth,
        type = "pairwise",
        strategies = c("phat1", "phat0")
    ) +
        ggplot2::scale_x_continuous(
            breaks = seq(0, 0.3, by = 0.05),
            labels = scales::label_percent(1)
        ) +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::theme(
            legend.position = "inside",
            legend.position.inside = c(0.85, 0.8),
            legend.text = ggplot2::element_text(size = 11)
        ) +
        ggplot2::labs(
            subtitle = "True delta NB (Model B vs model A)"
        ) +
        ggplot2::coord_cartesian(ylim = c(-0.05, 0.05)) +
        ggplot2::scale_y_continuous(breaks = seq(-0.05, 0.05, by = 0.025))

    ggplot2::ggsave(
        str_path("{outdir}/ground_truth_delta_nb.png"),
        dca_ground_truth_delta_plot,
        width = 9,
        height = 5.5,
        dpi = 600
    )

    # summarise results
    results_summary <- bayes_risk_simulation_results %>%
        dplyr::group_by(
            threshold, decision_strategy, sample_size, setting_id
        ) %>%
        dplyr::summarise(
            delta_above0 = mean(delta > 0),
            p_above50 = mean(prob > 0.5),
            p_above55 = mean(prob > 0.55),
            p_above60 = mean(prob > 0.6),
            p_above65 = mean(prob > 0.65),
            p_above70 = mean(prob > 0.7),
            p_above75 = mean(prob > 0.75),
            p_above80 = mean(prob > 0.8),
            p_above85 = mean(prob > 0.85),
            p_above90 = mean(prob > 0.9),
            p_above95 = mean(prob > 0.95),
            .groups = "drop"
        ) %>%
        dplyr::arrange(
            decision_strategy, sample_size, setting_id, threshold
        )

    # save results summary
    readr::write_tsv(
        results_summary,
        str_path("{outdir}/results_summary.tsv")
    )

    # save some ground truth metrics
    groud_truth_data <- dca_ground_truth$.data
    .calib_phat0 <- get_calibration_binary(groud_truth_data$outcomes, groud_truth_data$phat0)
    .calib_phat1 <- get_calibration_binary(groud_truth_data$outcomes, groud_truth_data$phat1)
    perf_summaries <- data.frame(
        model = c("phat0", "phat1"),
        true_prevalence = rep(mean(groud_truth_data$outcomes), 2),
        average_prediction = c(mean(groud_truth_data$phat0), mean(groud_truth_data$phat1)),
        calibration_oe = c(.calib_phat0$oe, .calib_phat1$oe),
        calibration_slope = c(.calib_phat0$slope, .calib_phat1$slope),
        auc = c(
            pROC::auc(pROC::roc(groud_truth_data$outcomes, groud_truth_data$phat0, quiet = TRUE)),
            pROC::auc(pROC::roc(groud_truth_data$outcomes, groud_truth_data$phat1, quiet = TRUE))
        )
    )

    readr::write_tsv(
        perf_summaries,
        str_path("{outdir}/ground_truth_simulation_metrics.tsv")
    )


    # post-process results
    rlang::inform("Plotting results")
    .odds <- function(p) p / (1 - p)
    # here the correct decision is always implement (threaholds >= 0.1 at least)
    risk_results <- results_summary %>%
        dplyr::filter(
            purrr::map_lgl(threshold, ~ any(dplyr::near(.x, seq(0.1, 0.3, by = 0.05))))
        ) %>%
        # implementation probability is the probability that a given implementation
        # decision rule spits out "implement" instead of "not implement"
        tidyr::pivot_longer(
            cols = dplyr::contains("above"),
            values_to = "implement_prob",
            names_to = "implementation_decision_rule"
        ) %>%
        dplyr::filter(
            implementation_decision_rule %in% c(
                "delta_above0",
                paste0("p_above", c(50, 80, 90))
            )
        ) %>%
        dplyr::mutate(
            # one for each risk aversion profile
            unscaled_risk_0.5 = if_else(
                decision_strategy == "phat0",
                implement_prob * .odds(0.5),
                1 - implement_prob
            ),
            unscaled_risk_0.8 = if_else(
                decision_strategy == "phat0",
                implement_prob * .odds(0.8),
                1 - implement_prob
            ),
            unscaled_risk_0.9 = if_else(
                decision_strategy == "phat0",
                implement_prob * .odds(0.9),
                1 - implement_prob
            )
        ) %>%
        tidyr::pivot_longer(
            cols = dplyr::contains("unscaled_risk"),
            names_to = "risk_aversion_parameter",
            names_transform = \(xx) as.numeric(stringr::str_extract(xx, "\\d+\\.\\d+")),
            values_to = "unscaled_risk"
        ) %>%
        dplyr::mutate(
            setting_label = if_else(
                decision_strategy == "phat1",
                "Implement",
                paste0("Not implement\n(", risk_aversion_parameter, ")")
            ),
            sample_size = factor(
                paste0("N=", sample_size),
                levels = paste0("N=", sort(unique(sample_size)))
            ),
            implementation_decision_rule = dplyr::case_match(
                implementation_decision_rule,
                "delta_above0" ~ "Point estimate > 0",
                "p_above50" ~ "P(best) > 50%",
                "p_above80" ~ "P(best) > 80%",
                "p_above90" ~ "P(best) > 90%"
            ),
            implementation_decision_rule = factor(
                implementation_decision_rule,
                levels = c(
                    "Point estimate > 0",
                    "P(best) > 50%",
                    "P(best) > 80%",
                    "P(best) > 90%"
                )
            )
        )

    freq_risk <- risk_results %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = threshold, y = unscaled_risk,
                fill = implementation_decision_rule
            )
        ) +
        ggplot2::geom_col(position = ggplot2::position_dodge()) +
        ggplot2::facet_grid(setting_label ~ sample_size, scales = "free") +
        ggplot2::scale_fill_brewer(palette = "Dark2") +
        ggplot2::labs(
            x = "Decision threshold",
            y = "Frequentist Risk (unscaled)",
            fill = "Implementation\ndecision rule"
        ) +
        ggplot2::theme_bw(base_size = 18)

    ggplot2::ggsave(
        str_path("{outdir}/frequentist_risk.png"),
        freq_risk,
        width = 17, height = 8.5, dpi = 600
    )

    bayes_risk_results <- risk_results %>%
        dplyr::select(threshold, sample_size, implementation_decision_rule, setting_label, unscaled_risk) %>%
        unique() %>%
        tidyr::pivot_wider(
            id_cols = c(threshold, sample_size, implementation_decision_rule),
            names_from = setting_label,
            values_from = unscaled_risk
        ) %>%
        dplyr::mutate(
            bayes_risk_uniform_0.5 = Implement + `Not implement\n(0.5)`,
            bayes_risk_uniform_0.8 = Implement + `Not implement\n(0.8)`,
            bayes_risk_uniform_0.9 = Implement + `Not implement\n(0.9)`,
            bayes_risk_skeptical_0.5 = Implement * 0.25 + `Not implement\n(0.5)` * 0.75,
            bayes_risk_skeptical_0.8 = Implement * 0.25 + `Not implement\n(0.8)` * 0.75,
            bayes_risk_skeptical_0.9 = Implement * 0.25 + `Not implement\n(0.9)` * 0.75,
            bayes_risk_optimistic_0.5 = Implement * 0.75 + `Not implement\n(0.5)` * 0.25,
            bayes_risk_optimistic_0.8 = Implement * 0.75 + `Not implement\n(0.8)` * 0.25,
            bayes_risk_optimistic_0.9 = Implement * 0.75 + `Not implement\n(0.9)` * 0.25
        ) %>%
        dplyr::select(
            threshold, sample_size, implementation_decision_rule,
            dplyr::contains("bayes_risk")
        ) %>%
        tidyr::pivot_longer(
            cols = dplyr::contains("bayes_risk"),
            names_to = "risk_aversion_parameter",
            values_to = "bayes_risk",
            names_transform = \(xx) paste0(
                stringr::str_extract(xx, "\\d+\\.\\d+"),
                "_",
                stringr::str_extract(xx, "uniform|skeptical|optimistic")
            )
        ) %>%
        dplyr::mutate(
            prior = stringr::str_extract(risk_aversion_parameter, "uniform|skeptical|optimistic"),
            risk_aversion_parameter = as.numeric(stringr::str_extract(risk_aversion_parameter, "\\d+\\.\\d+"))
        )
    bayes_risk_plot <- bayes_risk_results %>%
        dplyr::filter(prior == "uniform") %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = threshold, y = bayes_risk,
                fill = implementation_decision_rule
            )
        ) +
        ggplot2::geom_col(position = ggplot2::position_dodge()) +
        ggplot2::facet_grid(
            risk_aversion_parameter ~ sample_size,
            scales = "free",
            labeller = ggplot2::labeller(
                risk_aversion_parameter = function(x) {
                    dplyr::if_else(
                        x == 0.5,
                        "Risk-neutral\n\u03b3 = 0.5",
                        paste0("Risk-averse\n\u03b3 = ", x)
                    )
                }
            )
        ) +
        ggplot2::scale_fill_brewer(palette = "Dark2") +
        ggplot2::labs(
            x = "Decision threshold",
            y = "Bayes Risk",
            fill = "Implementation\ndecision rule"
        ) +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans())

    ggplot2::ggsave(
        str_path("{outdir}/bayes_risk.png"),
        bayes_risk_plot,
        width = 17, height = 8, dpi = 600
    )

    bayes_risk_plot_skeptical <- bayes_risk_results %>%
        dplyr::filter(prior == "skeptical") %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = threshold, y = bayes_risk,
                fill = implementation_decision_rule
            )
        ) +
        ggplot2::geom_col(position = ggplot2::position_dodge()) +
        ggplot2::facet_grid(risk_aversion_parameter ~ sample_size, scales = "free") +
        ggplot2::scale_fill_brewer(palette = "Dark2") +
        ggplot2::labs(
            x = "Decision threshold",
            y = "Bayes Risk (unscaled)",
            fill = "Implementation\ndecision rule"
        ) +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans())

    ggplot2::ggsave(
        str_path("{outdir}/bayes_risk_skeptical.png"),
        bayes_risk_plot_skeptical,
        width = 17, height = 8, dpi = 600
    )

    bayes_risk_plot_optimistic <- bayes_risk_results %>%
        dplyr::filter(prior == "optimistic") %>%
        ggplot2::ggplot(
            ggplot2::aes(
                x = threshold, y = bayes_risk,
                fill = implementation_decision_rule
            )
        ) +
        ggplot2::geom_col(position = ggplot2::position_dodge()) +
        ggplot2::facet_grid(risk_aversion_parameter ~ sample_size, scales = "free") +
        ggplot2::scale_fill_brewer(palette = "Dark2") +
        ggplot2::labs(
            x = "Decision threshold",
            y = "Bayes Risk (unscaled)",
            fill = "Implementation\ndecision rule"
        ) +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans())

    ggplot2::ggsave(
        str_path("{outdir}/bayes_risk_optimistic.png"),
        bayes_risk_plot_optimistic,
        width = 17, height = 8, dpi = 600
    )
}
