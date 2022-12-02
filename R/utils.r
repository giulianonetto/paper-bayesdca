import::from(rmda, rmda_data = dcaData)
import::from(dcurves, dcurves_dca = dca)
if (!require(bayesDCA)) {
  devtools::install_github("giulianonetto/bayesdca")
  require(bayesDCA)
}

#' Get maximum threshold allowed
validate_thresholds <- function(thresholds) {
  pmax(
    pmin(thresholds, 0.99),
    0
  )
}

#' Get path from `here::here` using `stringr::str_glue`
#'
#' @param .path A path potentially with glue-style
#' notation, e.g., "my/{path}/"
#' @return Character path as defined by `here::here` after
#' applying `stringr::str_glue` to the input path.
#' @examples
#' outdir <- path("R/output")
#' dir.create(outdir, recursive = TRUE)
#' path("{outdir}/my-new-file.pdf")
str_path <- function(.path, .envir = parent.frame()) {
  here::here(stringr::str_glue(.path, .envir = .envir))
}

#' Inverse logit (expit) function
#' @param x Real value (or vector) to be transformed into probability scale
#' @return Probability (or vector of)
expit <- function(x) {
  exp(x) / (1 + exp(x))
}

#' Compare DCA from `bayesDCa` and `rmda` packages
#'
#' @param dataset `data.frame` with outcomes and predictor of outcome
#' @param outcomes outcome variable (character string with column name)
#' @param predictor predicted probabilities (character string with column name)
#' @param bootstraps Number (int) of bootstrap samples for rmda
#' @param treat_all_rmda Logical indicating wether to plot Treat all from rmda (defaults to FALSE).
#' @param show_informative_prior Whether to show results from BayesDCA with informative prior.
#' @param refresh Refresh value for `rstan::sampling` (defaults to 0).
#' @param cores Number of cores for `bayesDCA::dca`. Defaults to 1.
#' @param thresholds Numeric vector (between 0 and 1) of thresholds for DCA.
#' @importFrom magrittr %>%
compare_bdca_vs_rmda <- function(dataset, outcomes,
                                 predictor, thresholds,
                                 treat_all_rmda = FALSE, bootstraps = 500,
                                 show_informative_prior = FALSE,
                                 refresh = 0, .quiet = FALSE, cores = 1) {
  df <- data.frame(
    outcomes = dataset[[outcomes]],
    predictor = dataset[[predictor]]
  )
  thresholds <- validate_thresholds(thresholds = thresholds)
  # Estimate decision curves
  if (isFALSE(.quiet)) {
    msg <- cli::col_blue("Estimating DCA with bayesDCA (default)")
    message(msg)
  }
  bdca_fit <- bayesDCA::dca(
    df,
    thresholds = thresholds,
    refresh = refresh,
    constant_prior = TRUE,
    cores = cores
  )

  if (isTRUE(show_informative_prior)) {
    if (isFALSE(.quiet)) {
      msg <- cli::col_blue("Estimating DCA with bayesDCA (thr-varying prior)")
      message(msg)
    }
    bdca_fit_informative <- bayesDCA::dca(
      df,
      thresholds = thresholds,
      refresh = refresh,
      constant_prior = FALSE,
      cores = cores
    )
  }

  if (isFALSE(.quiet)) {
    msg <- cli::col_blue("Estimating DCA with rmda")
    message(msg)
  }
  rmda_fit <- rmda::decision_curve(
    outcomes ~ predictor,
    data = df,
    fitted.risk = TRUE,
    thresholds = thresholds,
    bootstraps = bootstraps
  )

  # get results into standardized data.frames
  if (isFALSE(.quiet)) {
    msg <- cli::col_blue("Plotting results")
    message(msg)
  }
  res_bdca <- dplyr::bind_rows(
    bdca_fit$summary$net_benefit %>%
      dplyr::select(
        threshold, estimate,
        .lower := `2.5%`, .upper := `97.5%`
      ) %>%
      dplyr::mutate(
        .type = "Bayesian",
        strategy = "Model-based decisions"
      ),
    bdca_fit$summary$treat_all %>%
      dplyr::select(
        threshold, estimate,
        .lower := `2.5%`, .upper := `97.5%`
      ) %>%
      dplyr::mutate(
        .type = "Bayesian",
        strategy = "Treat all"
      )
  )

  res_rmda <- rmda_fit$derived.data %>%
    dplyr::filter(model %in% c("All", "outcomes ~ predictor")) %>%
    dplyr::mutate(
      strategy = ifelse(
        model == "All",
        "Treat all",
        "Model-based decisions"
      ),
      .type = "Frequentist"
    ) %>%
    dplyr::select(
      threshold := thresholds,
      estimate := NB,
      .lower := NB_lower,
      .upper := NB_upper,
      .type, strategy
    )

  res_all <- bind_rows(res_bdca, res_rmda)

  if (isTRUE(show_informative_prior)) {
    res_bdca_informative <- dplyr::bind_rows(
      bdca_fit_informative$summary$net_benefit %>%
        dplyr::select(
          threshold, estimate,
          .lower := `2.5%`, .upper := `97.5%`
        ) %>%
        dplyr::mutate(
          .type = "Bayesian (informative prior)",
          strategy = "Model-based decisions"
        ),
      bdca_fit_informative$summary$treat_all %>%
        dplyr::select(
          threshold, estimate,
          .lower := `2.5%`, .upper := `97.5%`
        ) %>%
        dplyr::mutate(
          .type = "Bayesian (informative prior)",
          strategy = "Treat all"
        )
    )
    res_all <- dplyr::bind_rows(res_all, res_bdca_informative)
  }
  return(res_all)
}

#' Plot bayesDCA vs rmda comparison
#'
#' Plots DCA from `bayesDCa` and `rmda` packages
#' for visual comparison
#'
#' @param comparison output from `compare_bdca_vs_rmda`.
#' May also leave as `NULL` and pass arguments for the
#' `compare_bdca_vs_rmda` function which is run internally.
#' @param dataset `data.frame` with outcomes and predictor of outcome
#' @param outcomes outcome variable (character string with column name)
#' @param predictor predicted probabilities (character string with column name)
#' @param bootstraps Number`data.frame` with outcomes and predictor of outcome (int) of bootstrap samples for rmda
#' @param show_informative_prior Whether to show results from informative prior Bayesian DCA
#' @param treat_all_rmda Logical indicating wether to plot Treat all from rmda (defaults to FALSE).
#' @param refresh Refresh value for `rstan::sampling` (defaults to 0).
#' @param cores Number of cores for `bayesDCA::dca`. Defaults to 1.
#' @thresholds Numeric vector (between 0 and 1) of thresholds for DCA.
#' @importFrom magrittr %>%
plot_bdca_vs_rmda <- function(comparison = NULL,
                              dataset = NULL, outcomes = NULL,
                              predictor = NULL, thresholds = NULL,
                              treat_all_rmda = FALSE, bootstraps = 500,
                              show_informative_prior = FALSE,
                              refresh = 0, cores = 1, .quiet = FALSE) {
  import::from(magrittr, `%>%`)

  if (is.null(comparison)) {
    comparison <- compare_bdca_vs_rmda(
      dataset = dataset, outcomes = outcomes,
      predictor = predictor, thresholds = thresholds,
      treat_all_rmda = FALSE, bootstraps = 500,
      show_informative_prior = show_informative_prior,
      refresh = 0, cores = cores, .quiet = FALSE
    )
  }
  # get plot helper objects
  max_estimate <- max(comparison$estimate)
  .ymin <- ifelse(
    max_estimate > 0.02,
    -0.02,
    -max_estimate
  )

  .cols <- c(
    "Model-based decisions.Bayesian" = "#1B9E77",
    "Model-based decisions.Frequentist" = "red",
    "Treat all.Bayesian" = "gray40"
  )
  .labels <- c(
    "Model-based decisions (Bayesian)",
    "Model-based decisions (Frequentist)",
    "Treat all"
  )

  if (isTRUE(treat_all_rmda)) {
    .cols["Treat all.Frequentist"] <- "red"
    .labels[3] <- "Treat all (Bayesian)"
    .labels[4] <- "Treat all (Frequentist)"
  } else {
    comparison <- comparison %>%
      dplyr::filter(
        !(strategy == "Treat all" & .type == "Frequentist")
      )
  }

  .plot <- comparison %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = threshold, y = estimate,
        ymin = .lower, ymax = .upper,
        color = interaction(strategy, .type),
        fill = interaction(strategy, .type)
      )
    ) +
    ggplot2::geom_ribbon(alpha = 0.2, aes(color = NULL)) +
    ggplot2::geom_line(aes(lty = .type), show.legend = FALSE) +
    ggplot2::geom_hline(
      yintercept = 0,
      lty = "longdash",
      alpha = 0.5
    ) +
    ggplot2::coord_cartesian(ylim = c(.ymin, NA)) +
    ggplot2::scale_color_manual(
      values = .cols,
      labels = .labels
    ) +
    ggplot2::scale_fill_manual(
      values = .cols,
      labels = .labels
    ) +
    ggplot2::scale_linetype_manual(
      values = c("solid", "longdash")
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(1)
    ) +
    ggplot2::labs(
      color = NULL, fill = NULL, lty = NULL,
      y = "Net benefit", x = "Decision threshold"
    )

  return(.plot)
}

#' Compute net benefit for given threshold
#' @param y
#' @param pred
#' @param thr
compute_nb <- function(y, pred, thr) {
  tp_ix <- pred >= thr & y == 1
  fp_ix <- pred >= thr & y == 0
  tpr <- mean(tp_ix)
  fpr <- mean(fp_ix)
  tibble::tibble(
    .thr = thr,
    nb = tpr - fpr * (thr / (1 - thr))
  )
}

#' Simulate population data for DCA simulation (Results's subsection 02)
#'
simulate_dca_population <- function(true_beta, beta_hat, n_pop, thresholds, .seed, .verbose = FALSE) {
  thresholds <- validate_thresholds(thresholds = thresholds)
  msg <- cli::col_blue(
    paste0(
      "Simulating population DCA data with population seed = ", .seed
    )
  )
  message(msg)
  d <- length(true_beta) - 1
  ## simulate predictors x, even prob true_p|x, and outcome y|true_p
  set.seed(.seed)
  x <- cbind(1, matrix(rexp(n_pop * d, 1), ncol = d))
  true_p <- plogis(as.vector(x %*% true_beta))
  set.seed(.seed)
  y <- rbinom(n_pop, 1, true_p)
  ## calculate predictions p_hat|x
  p_hat <- plogis(as.vector(x %*% beta_hat))
  if (isTRUE(.verbose)) {
    msg <- cli::col_blue(
      paste0(
        "\tTrue prevalence: ", round(mean(true_p), 3),
        "\n\tPrevalence implied by model: ", round(mean(p_hat), 3)
      )
    )
    message(msg)
  }

  ## calculate true NB|y,p_hat
  true_nb <- map_df(thresholds, ~ {
    compute_nb(y = y, pred = p_hat, thr = .x)
  })

  output <- list(
    y = y,
    x = x,
    p_hat = p_hat,
    true_p = true_p,
    thresholds = thresholds,
    n_pop = n_pop,
    true_nb = true_nb,
    true_beta = true_beta,
    beta_hat = beta_hat,
    population_seed = .seed
  )

  return(output)
}

get_setting_sample_list <- function(events, n_sim, population_data, .setting_seed, .setting_label) {
  sample_size <- ceiling(events / mean(population_data$true_p))
  set.seed(.setting_seed)
  sim_seeds <- sample(1:2e6, n_sim)

  df_sample_list <- lapply(1:n_sim, function(j) {
    .run_seed <- sim_seeds[j]
    .run_label <- paste0(
      .setting_label, "_", .setting_seed,
      "-run", j, "_", .run_seed
    )
    set.seed(.run_seed)
    sample_ix <- sample(population_data$n_pop, sample_size)
    .df_sample <- data.frame(
      outcomes = population_data$y[sample_ix],
      model_predictions = population_data$p_hat[sample_ix],
      setting_label = .setting_label,
      simulation_run_label = .run_label,
      .setting_seed = .setting_seed,
      .run_seed = .run_seed,
      .run_id = j
    )
    return(.df_sample)
  })

  return(df_sample_list)
}

#' Simulate DCA (Results' subsection 02)
#' @param population_data Output from `simulate_dca_population`
#' @param thresholds DCA decision thresholds
#' @param events Expected number of events in DCA sample
#' @param .seed RNG seed.
#' @param raw_data Whether to return the raw data. Defaults to FALSE.
#' @param .plot Whether to plot simulation. Defaults to FALSE.
#' @param result_path If given, the simulation result will be written as a .tsv file to the specified path.
#' @param .run_label If given, it will be appended to result as a "simulation_label" column.
#' @param .setting_label If given, it will be appended to result as a "setting_label" column.
#' @param overwrite If TRUE, any existing file in `result_path` will be overwritten by a new simualtion.
#' @param cores Number of cores for bayesDCA. Defaults to 1.
#' @param .verbose If TRUE more info is printed.
run_dca_simulation <- function(df_sample,
                               thresholds,
                               true_nb,
                               true_prevalence,
                               raw_data = FALSE, .plot = FALSE,
                               result_path = NULL, .run_label = NULL,
                               .setting_label = NULL,
                               overwrite = FALSE, cores = 1, .verbose = FALSE) {
  if (!is.null(result_path)) {
    if (file.exists(result_path) && isFALSE(overwrite)) {
      msg <- cli::col_red(paste0(
        "Skipping simulation, cannot overwrite: ", result_path
      ))
      message(msg)
      output <- list(
        result = read_tsv(result_path, show_col_types = FALSE),
        thresholds = thresholds
      )
      return(output)
    } else {
      dir.create(
        dirname(result_path),
        showWarnings = FALSE,
        recursive = TRUE
      )
    }
  }

  dca_comparison <- compare_bdca_vs_rmda(
    dataset = df_sample,
    outcomes = "outcomes",
    predictor = "model_predictions",
    thresholds = thresholds,
    show_informative_prior = TRUE,
    bootstraps = 1e3,
    .quiet = TRUE,
    cores = cores
  )

  result <- left_join(
    dca_comparison,
    true_nb %>%
      dplyr::select(
        threshold := .thr,
        .true_nb := nb
      ),
    by = "threshold"
  ) %>%
    dplyr::mutate(
      abs_error = abs(estimate - .true_nb),
      truth_within_interval = .true_nb >= .lower & .true_nb <= .upper,
      true_prevalence = true_prevalence,
      .simulation_seed = .seed
    )

  if (!is.null(.setting_label)) {
    result$setting_label <- unique(df_sample$setting_label)
  }
  if (!is.null(.run_label)) {
    result$simulation_run_label <- unique(df_sample$simulation_run_label)
  }


  if (!is.null(result_path)) {
    write_tsv(
      result,
      result_path,
    )
  }

  output <- list(
    result = result,
    thresholds = thresholds
  )

  if (isTRUE(.plot)) {
    sim_plot <- plot_bdca_vs_rmda(
      comparison = dca_comparison
    ) +
      geom_line(
        data = true_nb,
        aes(x = .thr, y = nb, group = 1),
        color = "red", inherit.aes = FALSE
      ) +
      theme(legend.position = c(.7, .7)) +
      ggplot2::scale_y_continuous(expand = c(.275, 0))
    output[["plot"]] <- sim_plot
  }

  if (isTRUE(raw_data)) {
    output[["true_nb"]] <- true_nb
    output[["df_sample"]] <- df_sample
  }

  return(output)
}

#' Get simulation settings (Results' subsection 2)
get_simulation_settings <- function() {
  simulation_settings <- list(
    # AUC 0.65, prev 0.01 vs 0.0096
    sim1 = list(
      true_beta = c(-4.75, -log(1.5), log(1.5)),
      beta_hat = c(-5, -log(1.5) * 1.25, log(1.5) * 1.25),
      auc = 0.65,
      prev = 0.01
    ),
    # AUC 0.65, prev 0.05 vs 0.06
    sim2 = list(
      true_beta = c(-3.1, -log(1.5), log(1.5)),
      beta_hat = c(-3.9, -log(1.5) * 3, log(1.5) * 3),
      auc = 0.65,
      prev = 0.05
    ),
    # AUC 0.65, prev 0.3 vs 0.3
    sim3 = list(
      true_beta = c(-0.9, -log(1.55), log(1.55)),
      beta_hat = c(-1.2, -log(1.55) * 3, log(1.55) * 3),
      auc = 0.65,
      prev = 0.3
    ),
    # AUC 0.85, prev  0.01 vs 0.009
    sim4 = list(
      true_beta = c(-5.6, -log(2.57), log(2.57)),
      beta_hat = c(-6.9, -log(2.57) * 1.5, log(2.57) * 1.5),
      auc = 0.85,
      prev = 0.01
    ),
    # AUC 0.85, prev  0.05 vs 0.06
    sim5 = list(
      true_beta = c(-3.755, -log(2.95), log(2.95)),
      beta_hat = c(-7.3, -log(2.95) * 3, log(2.95) * 3),
      auc = 0.85,
      prev = 0.05
    ),
    # AUC 0.85, prev  0.3 vs 0.32
    sim6 = list(
      true_beta = c(-1.3, -log(4.5), log(4.5)),
      beta_hat = c(-2.25, -log(4.5) * 3, log(4.5) * 3),
      auc = 0.85,
      prev = 0.3
    )
  )

  return(simulation_settings)
}

test_sim_setting <- function(.sim, n = 1e5, .return = FALSE) {
  x <- cbind(1, rexp(n, 1), rexp(n, 1))
  p <- as.vector(
    plogis(x %*% .sim$true_beta)
  )
  phat <- as.vector(
    plogis(x %*% .sim$beta_hat)
  )
  y <- rbinom(n, 1, p)
  suppressMessages({
    suppressWarnings({
      a <- round(pROC::auc(y, p)[1], 2)
      b <- round(pROC::auc(y, phat)[1], 2)
    })
  })
  print(stringr::str_glue("Claimed AUC: {.sim$auc}"))
  print(stringr::str_glue("Claimed prev: {.sim$prev}"))
  print(stringr::str_glue("AUC p: {a}"))
  print(stringr::str_glue("AUC phat: {b}"))
  print(stringr::str_glue("Prevalence y: {round(mean(y), 3)}"))
  print(stringr::str_glue("Prevalence phat: {round(mean(phat), 3)}"))
  print(stringr::str_glue("High preds prop: {round(mean(phat > .3), 5)}"))
  print(stringr::str_glue("MAX p: {round(max(p), 2)}"))
  print(stringr::str_glue("MAX phat: {round(max(phat), 2)}"))
  if (isTRUE(.return)) {
    return(list(p = p, phat = phat, y = y))
  }
}
