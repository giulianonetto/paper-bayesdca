import::from(rmda, rmda_data = dcaData)
import::from(dcurves, dcurves_dca = dca)
if (!require(bayesDCA)) {
  stop("Make sure `bayesDCA` is installed.")
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
#' @param thresholds Numeric vector (between 0 and 1) of thresholds for DCA.
#' @importFrom magrittr %>%
compare_bdca_vs_rmda <- function(dataset, outcomes,
                                 predictor, thresholds,
                                 treat_all_rmda = FALSE, bootstraps = 500,
                                 show_informative_prior = FALSE,
                                 refresh = 0, .quiet = FALSE) {
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
  t0 <- proc.time()["elapsed"]
  bdca_fit <- bayesDCA::dca(
    df,
    thresholds = thresholds
  )
  time_bayes <- proc.time()["elapsed"] - t0

  if (isTRUE(show_informative_prior)) {
    if (isFALSE(.quiet)) {
      msg <- cli::col_blue("Estimating DCA with bayesDCA (thr-varying prior)")
      message(msg)
    }
    t0 <- proc.time()["elapsed"]
    bdca_fit_informative <- bayesDCA::dca(
      df,
      thresholds = thresholds
    )
    time_bayes2 <- proc.time()["elapsed"] - t0
  }

  if (isFALSE(.quiet)) {
    msg <- cli::col_blue("Estimating DCA with rmda")
    message(msg)
  }
  t0 <- proc.time()["elapsed"]
  rmda_fit <- rmda::decision_curve(
    outcomes ~ predictor,
    data = df,
    fitted.risk = TRUE,
    thresholds = thresholds,
    bootstraps = bootstraps
  )
  time_freq <- proc.time()["elapsed"] - t0

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
  ) %>%
    dplyr::mutate(runtime = time_bayes)

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
    ) %>%
    dplyr::mutate(runtime = time_freq)

  res_all <- dplyr::bind_rows(res_bdca, res_rmda)

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
    ) %>%
      dplyr::mutate(runtime = time_bayes2)

    res_all <- dplyr::bind_rows(res_all, res_bdca_informative)
  }
  return(res_all)
}

#' Extract posterior summaries for `BayesDCASurv` object
summarise_bdca_surv_pars <- function(fit) {
  surv <- rstan::summary(fit$fit, pars = "St_positives")$summary %>%
    dplyr::as_tibble(rownames = "par_name") %>%
    dplyr::select(-c(se_mean, sd, n_eff, Rhat)) %>% # nolint
    dplyr::rename(surv_estimate := mean) %>% # nolint
    dplyr::mutate(
      threshold_ix = stringr::str_extract(par_name, "\\d+\\]") %>% # nolint
        stringr::str_remove(string = ., pattern = "\\]") %>%
        as.integer(),
      strategy_ix = stringr::str_extract(par_name, "\\[\\d+") %>% # nolint
        stringr::str_remove(string = ., pattern = "\\[") %>%
        as.integer(),
      threshold = fit$thresholds[threshold_ix],
      decision_strategy_name = fit$strategies[strategy_ix],
      par_name = stringr::str_extract(par_name, "\\w+")
    ) %>%
    dplyr::select(
      threshold, decision_strategy_name, surv_estimate,
      surv_lower := `2.5%`, surv_upper := `97.5%`
    )
  pos <- fit$summary$positivity %>%
    dplyr::select(
      threshold, decision_strategy_name,
      pos_estimate := estimate,
      pos_lower := `2.5%`, pos_upper := `97.5%`
    )
  dplyr::inner_join(surv, pos, by = c("threshold", "decision_strategy_name"))
}

#' Get results from bayesDCA surv
get_results_surv_bdca <- function(fit, type_label) {
  nb_summary <- fit$summary$net_benefit %>%
    dplyr::select(
      threshold, estimate,
      .lower := `2.5%`, .upper := `97.5%`
    ) %>%
    dplyr::mutate(
      .type = type_label,
      strategy = "Model-based decisions"
    )
  nb_pars_summary <- summarise_bdca_surv_pars(fit)
  tall_summary <- fit$summary$treat_all %>%
    dplyr::select(
      threshold, estimate,
      .lower := `2.5%`, .upper := `97.5%`
    ) %>%
    dplyr::mutate(
      .type = type_label,
      strategy = "Treat all"
    )
  tall_pars_summary <- fit$summary$overall_surv %>%
    dplyr::select(
      surv_estimate := estimate,
      surv_lower := `2.5%`, surv_upper := `97.5%`
    ) %>%
    dplyr::mutate(
      pos_estimate = 1, pos_lower = 1, pos_upper = 1
    )
  res_bdca <- dplyr::bind_rows(
    dplyr::left_join(nb_summary, nb_pars_summary, by = "threshold"),
    dplyr::bind_cols(tall_summary, tall_pars_summary)
  )

  return(res_bdca)
}


#' Compare DCA (Survival) from `bayesDCa` and `dcurves` packages
#'
#' @param dataset `data.frame` with outcomes and predictor of outcome
#' @param outcomes outcome variable (character string with column name)
#' @param predictor predicted probabilities (character string with column name)
#' @param pred_time Prediction time horizon
#' @param treat_all_rmda Logical indicating wether to plot Treat all from rmda (defaults to FALSE).
#' @param refresh Refresh value for `rstan::sampling` (defaults to 0).
#' @param cores Number of cores for `bayesDCA::dca`. Defaults to 1.
#' @param thresholds Numeric vector (between 0 and 1) of thresholds for DCA.
#' @importFrom magrittr %>%
compare_bdca_vs_dcurves <- function(dataset, outcomes,
                                    predictor, thresholds,
                                    pred_time,
                                    treat_all_rmda = FALSE, bootstraps = 500,
                                    refresh = 0, .quiet = FALSE, cores = 1) {
  df <- data.frame(
    outcomes = dataset[[outcomes]],
    predictor = dataset[[predictor]]
  )
  .event_times <- unname(df[["outcomes"]][, 1])[unname(df[["outcomes"]][, 2]) == 1L]
  thresholds <- validate_thresholds(thresholds = thresholds)
  # Estimate decision curves
  if (isFALSE(.quiet)) {
    msg <- cli::col_blue("Estimating DCA with bayesDCA (default)")
    message(msg)
  }

  t0 <- proc.time()["elapsed"]
  bdca_fit <- try(
    {
      bayesDCA::dca_surv(
        df,
        thresholds = thresholds,
        prediction_time = pred_time,
        keep_fit = TRUE,
        shape_prior_pars = c(5, 0, 1.5),
        refresh = refresh,
        cores = cores
      )
    },
    silent = TRUE
  )
  time_bayes <- proc.time()["elapsed"] - t0

  if (isFALSE(.quiet)) {
    msg <- cli::col_blue("Estimating DCA with dcurves")
    message(msg)
  }

  t0 <- proc.time()["elapsed"]
  dcurves_fit <- try(
    {
      fit_bootstrap_dcurves(
        formula = outcomes ~ predictor,
        data = df,
        time = pred_time,
        thresholds = thresholds,
        B = bootstraps
      )
    },
    silent = TRUE
  )
  time_freq <- proc.time()["elapsed"] - t0

  # get results into standardized data.frames
  if (isFALSE(.quiet)) {
    msg <- cli::col_blue("Plotting results")
    message(msg)
  }

  if (isFALSE(inherits(bdca_fit, "try-error"))) {
    res_bdca <- get_results_surv_bdca(fit = bdca_fit, type_label = "Bayesian") %>%
      dplyr::mutate(runtime = time_bayes)
  } else {
    res_bdca <- data.frame()
  }

  if (isFALSE(inherits(dcurves_fit, "try-error"))) {
    res_dcurves <- dcurves_fit %>%
      dplyr::filter(variable %in% c("all", "predictor")) %>%
      dplyr::mutate(
        strategy = ifelse(
          variable == "all",
          "Treat all",
          "Model-based decisions"
        ),
        .type = "Frequentist",
        .lower := net_benefit - 1.96 * se,
        .upper := net_benefit + 1.96 * se,
      ) %>%
      dplyr::select(
        threshold,
        estimate := net_benefit,
        .lower, .upper,
        .type, strategy
      ) %>%
      dplyr::mutate(runtime = time_freq)
  } else {
    res_dcurves <- data.frame()
  }

  res_all <- dplyr::bind_rows(res_bdca, res_dcurves)
  return(res_all)
}

#' Fit bootstrapped version of dcurves
#' @param .formula
#' @param data
#' @param time
#' @param thresholds
fit_bootstrap_dcurves <- function(
    formula,
    data,
    time,
    thresholds,
    B = 500) {
  res <- lapply(
    seq_len(B),
    function(i) {
      ix <- sample(nrow(data), replace = TRUE)
      x <- try(
        {
          dcurves::dca(
            formula,
            data = data[ix, ],
            time = time,
            thresholds = thresholds
          )$dca %>%
            dplyr::select(variable, label, threshold, net_benefit)
        },
        silent = TRUE
      )
      if (isFALSE(inherits(x, "try-error"))) {
        return(x)
      } else {
        return(NULL)
      }
    }
  )
  res <- dplyr::bind_rows(res, .id = "id")

  se <- res %>%
    dplyr::group_by(variable, label, threshold) %>%
    dplyr::summarise(
      se = sd(net_benefit, na.rm = TRUE),
      na_prop = mean(is.na(net_benefit)),
      .groups = "drop"
    )
  if (any(se$na_prop > 0.5)) {
    print(se)
    msg <- "more than 50% NAs in dcurves bootstrap for some case."
    logger::log_warn(cli::col_br_red(msg))
    if (any(se$na_prop >= 0.95)) {
      se$se[se$na_prop >= 0.95] <- NA_real_ # bootstrap failed
    }
  }

  se <- se %>% dplyr::select(-na_prop)
  estimate <- dcurves::dca(
    formula,
    data = data,
    time = time,
    thresholds = thresholds
  )$dca

  output <- dplyr::inner_join(estimate, se, by = c("variable", "label", "threshold"))
  return(output)
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
#' @thresholds Numeric vector (between 0 and 1) of thresholds for DCA.
#' @importFrom magrittr %>%
plot_bdca_vs_rmda <- function(comparison = NULL,
                              dataset = NULL, outcomes = NULL,
                              predictor = NULL, thresholds = NULL,
                              treat_all_rmda = FALSE, bootstraps = 500,
                              show_informative_prior = FALSE,
                              refresh = 0, .quiet = FALSE) {
  import::from(magrittr, `%>%`)

  if (is.null(comparison)) {
    comparison <- compare_bdca_vs_rmda(
      dataset = dataset, outcomes = outcomes,
      predictor = predictor, thresholds = thresholds,
      treat_all_rmda = FALSE, bootstraps = bootstraps,
      show_informative_prior = show_informative_prior,
      refresh = 0, .quiet = FALSE
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
    "Model-based decisions.Bayesian" = "#7570B3",
    "Model-based decisions.Frequentist" = "#D95F02",
    "Treat all.Bayesian" = "gray20"
  )
  .labels <- c(
    "Model-based decisions.Bayesian" = "Model-based decisions (Bayesian)",
    "Model-based decisions.Frequentist" = "Model-based decisions (Frequentist)",
    "Treat all.Bayesian" = "Treat all"
  )

  if (isTRUE(treat_all_rmda)) {
    .cols["Treat all.Frequentist"] <- "red"
    .labels["Treat all.Bayesian"] <- "Treat all (Bayesian)"
    .labels["Treat all.Frequentist"] <- "Treat all (Frequentist)"
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
    ggplot2::geom_ribbon(alpha = 0.3, aes(color = NULL)) +
    ggplot2::geom_line(aes(lty = .type), show.legend = FALSE) +
    ggplot2::geom_hline(
      yintercept = 0,
      lty = "longdash",
      alpha = 0.5
    ) +
    ggplot2::coord_cartesian(ylim = c(.ymin, NA)) +
    ggplot2::scale_color_manual(
      values = .cols,
      labels = .labels,
      breaks = names(.cols)
    ) +
    ggplot2::scale_fill_manual(
      values = .cols,
      labels = .labels,
      breaks = names(.cols)
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

#' Get calibration for binary outcomes
#' @param y
#' @param pred
get_calibration_binary <- function(y, pred) {
  oe <- mean(y) / mean(pred)
  lp <- qlogis(pred)
  ix <- !is.infinite(lp) & !is.na(lp)
  lp <- lp[ix]
  y <- y[ix]
  slope <- coef(glm(y ~ 1 + lp, family = binomial(link = "logit")))[[2]]
  return(list(oe = oe, slope = slope))
}

#' Get performance (calibration + concordance) for survival outcomes
#' @param surv_time The true simulated survival times
#' @param pred Predicted 12-month event probability from example model being validated
get_performance_surv <- function(surv_time, obs_time, status, lp, predicted_risk, pred_time = 12) {
  event_rate <- mean(surv_time < pred_time)
  oe <- event_rate / mean(predicted_risk)
  fit <- survival::coxph(survival::Surv(obs_time, status) ~ lp)
  slope <- coef(fit)[[1]]
  concordance <- survival:::concordance(fit)[[1]]
  time_roc_oracle <- pROC::auc(
    surv_time > pred_time,
    predicted_risk,
    quiet = TRUE
  )[[1]]

  output <- list(
    oe = oe, slope = slope, concordance = concordance,
    time_roc = NA, time_roc_oracle = time_roc_oracle
  )
  return(output)
}

#' Compute net benefit for given threshold
#' @param surv_time The true simulated survival times
#' @param surv_hat Predicted 1-year event probability from example model being validated
#' @param thr Decision threshold
#' @param pred_time Prediction time horizon
compute_nb_surv <- function(surv_time, p_hat, thr, pred_time = 1) {
  pos <- mean(p_hat >= thr) # prob of positive prediction
  cond_surv <- mean(surv_time[p_hat >= thr] >= pred_time) # conditional survival
  tibble::tibble(
    .thr = thr,
    positivity = pos,
    conditional_survival = cond_surv,
    nb = (1 - cond_surv) * pos - cond_surv * pos * (thr / (1 - thr))
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

  ## calculate performance metrics
  calibration <- get_calibration_binary(y = y, pred = p_hat)
  auc <- pROC::auc(y, p_hat, quiet = TRUE)[[1]]

  output <- list(
    y = y,
    x = x,
    p_hat = p_hat,
    true_p = true_p,
    thresholds = thresholds,
    n_pop = n_pop,
    true_nb = true_nb,
    calibration_oe = calibration$oe,
    calibration_slope = calibration$slope,
    auc = auc,
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

get_setting_sample_list_surv <- function(events, n_sim, population_data, .setting_seed, .setting_label) {
  sample_size <- ceiling(events / mean(population_data$df_pop$status))
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
    d <- population_data$df_pop[sample_ix, ]
    .df_sample <- data.frame(
      survTime = d$survTime,
      obsTime = d$obsTime,
      status = d$status,
      model_predictions = d$p_hat,
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

#' Simulate DCA (Results' subsection 2.5.1)
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
#' @param .verbose If TRUE more info is printed.
run_dca_simulation <- function(df_sample,
                               thresholds,
                               true_nb,
                               true_prevalence,
                               raw_data = FALSE, .plot = FALSE,
                               result_path = NULL, .run_label = NULL,
                               .setting_label = NULL,
                               overwrite = FALSE, .verbose = FALSE) {
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
    bootstraps = 500,
    .quiet = TRUE
  )

  result <- dplyr::left_join(
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
    readr::write_tsv(
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
      ggplot2::geom_line(
        data = true_nb,
        aes(x = .thr, y = nb, group = 1),
        color = "red", inherit.aes = FALSE
      ) +
      ggplot2::theme(legend.position = c(.7, .7)) +
      ggplot2::scale_y_continuous(expand = c(.275, 0))
    output[["plot"]] <- sim_plot
  }

  if (isTRUE(raw_data)) {
    output[["true_nb"]] <- true_nb
    output[["df_sample"]] <- df_sample
  }

  return(output)
}

#' Simulate DCA Surv (Results' subsection 2.5.2)
#' @param df_sample Output from `simulate_dca_population`
#' @param thresholds DCA decision thresholds
#' @param true_nb Expected number of events in DCA sample
#' @param true_incidence
#' @param .seed RNG seed.
#' @param raw_data Whether to return the raw data. Defaults to FALSE.
#' @param .plot Whether to plot simulation. Defaults to FALSE.
#' @param result_path If given, the simulation result will be written as a .tsv file to the specified path.
#' @param .run_label If given, it will be appended to result as a "simulation_label" column.
#' @param .setting_label If given, it will be appended to result as a "setting_label" column.
#' @param overwrite If TRUE, any existing file in `result_path` will be overwritten by a new simualtion.
#' @param cores Number of cores for bayesDCA. Defaults to 1.
#' @param .verbose If TRUE more info is printed.
run_dca_simulation_surv <- function(df_sample,
                                    thresholds,
                                    true_nb,
                                    true_incidence,
                                    pred_time,
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

  dca_comparison <- df_sample %>%
    dplyr::filter(obsTime > 0) %>%
    dplyr::mutate(
      outcomes = survival::Surv(obsTime, status),
    ) %>%
    dplyr::select(outcomes, model_predictions) %>%
    compare_bdca_vs_dcurves(
      outcomes = "outcomes",
      predictor = "model_predictions",
      thresholds = thresholds,
      pred_time = pred_time,
      .quiet = TRUE,
      cores = cores
    )

  result <- dplyr::left_join(
    dca_comparison,
    true_nb %>%
      dplyr::select(
        threshold := .thr,
        .true_nb := nb,
        .true_pos := positivity,
        .true_surv := conditional_survival
      ),
    by = "threshold"
  ) %>%
    dplyr::mutate(
      abs_error = abs(estimate - .true_nb),
      truth_within_interval = .true_nb >= .lower & .true_nb <= .upper,
      truth_within_interval = ifelse(
        is.na(truth_within_interval),
        FALSE,
        truth_within_interval
      ),
      truth_within_interval_surv = .true_surv >= surv_lower & .true_surv <= surv_upper,
      truth_within_interval_pos = .true_pos >= pos_lower & .true_pos <= pos_upper,
      error_surv = surv_estimate - .true_surv,
      error_pos = pos_estimate - .true_pos,
      true_incidence = true_incidence,
      .simulation_seed = .seed
    ) %>%
    dplyr::select(
      threshold, estimate, .lower, .upper, .type, strategy,
      .true_nb, abs_error, truth_within_interval,
      dplyr::contains("surv_"), dplyr::contains("_surv"),
      dplyr::contains("pos_"), dplyr::contains("_pos"), runtime,
      dplyr::everything()
    )

  if (!is.null(.setting_label)) {
    result$setting_label <- unique(df_sample$setting_label)
  }
  if (!is.null(.run_label)) {
    result$simulation_run_label <- unique(df_sample$simulation_run_label)
  }

  if (!is.null(result_path)) {
    readr::write_tsv(
      result,
      result_path,
    )
  }

  output <- list(
    result = result,
    thresholds = thresholds
  )

  if (isTRUE(.plot)) {
    .cols <- c(
      "Model-based decisions.Bayesian" = "#1B9E77",
      "Model-based decisions.Frequentist" = "red",
      "Treat all.Bayesian" = "gray40"
    )
    .labels <- c(
      "Model-based decisions.Bayesian" = "Model-based decisions (Bayesian)",
      "Model-based decisions.Frequentist" = "Model-based decisions (Frequentist)",
      "Treat all.Bayesian" = "Treat all"
    )

    sim_plot <- dca_comparison %>%
      dplyr::filter(
        !(strategy == "Treat all" & .type == "Frequentist")
      ) %>%
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
      ggplot2::coord_cartesian(ylim = c(-.01, NA)) +
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
      ) +
      ggplot2::geom_line(
        data = true_nb,
        aes(x = .thr, y = nb, group = 1),
        color = "blue", inherit.aes = FALSE
      ) +
      ggplot2::theme(legend.position = c(.7, .7)) +
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

#' Get simulation survival settings (Results' subsection 3)
get_simulation_settings_surv <- function() {
  # gamma is shape
  # lambda is scale
  simulation_settings <- list(
    sim1 = list(
      true_beta = c(log(1.3), log(0.7)),
      beta_hat = c(log(1.3), log(0.7)) * 1.01,
      lambda = 0.12,
      gamma = 1.22,
      max_follow_up = 12 * 2,
      event_fraction = 0.77,
      one_year_survival_rate = 0.1,
      concord = 0.6,
      median_surv = 4 # months
    ),
    sim2 = list(
      true_beta = c(log(1.3), log(0.7)),
      beta_hat = c(log(1.3), log(0.7)) * 1.01,
      lambda = 0.12,
      gamma = 1.07,
      max_follow_up = 12 * 2,
      event_fraction = 0.7,
      one_year_survival_rate = 0.2,
      concord = 0.6,
      median_surv = 5 # months
    ),
    sim3 = list(
      true_beta = c(log(1.3), log(0.7)),
      beta_hat = c(log(1.3), log(0.7)) * 1.01,
      lambda = 0.12,
      gamma = 0.7,
      max_follow_up = 12 * 2,
      event_fraction = 0.46,
      one_year_survival_rate = 0.5,
      concord = 0.6,
      median_surv = 12 # months
    ),
    sim4 = list(
      true_beta = c(log(1.95), log(0.05)),
      beta_hat = c(log(1.95), log(0.05)) * 1.25,
      lambda = 4e-4,
      gamma = 4.6,
      max_follow_up = 12 * 2,
      event_fraction = 0.75,
      one_year_survival_rate = 0.1,
      concord = 0.9,
      median_surv = 5 # months
    ),
    sim5 = list(
      true_beta = c(log(1.95), log(0.05)),
      beta_hat = c(log(1.95), log(0.05)) * 1.25,
      lambda = 4e-4,
      gamma = 4,
      max_follow_up = 12 * 2,
      event_fraction = 0.67,
      one_year_survival_rate = 0.2,
      concord = 0.9,
      median_surv = 6 # months
    ),
    sim6 = list(
      true_beta = c(log(1.95), log(0.05)),
      beta_hat = c(log(1.95), log(0.05)) * 1.25,
      lambda = 4e-4,
      gamma = 2.9,
      max_follow_up = 12 * 2,
      event_fraction = 0.44,
      one_year_survival_rate = 0.5,
      concord = 0.9,
      median_surv = 12 # months
    ),
    sim7 = list(
      true_beta = c(log(1.95), log(0.001)),
      beta_hat = c(log(1.95), log(0.001)) * 1.25,
      lambda = 4e-4,
      gamma = 6.5,
      max_follow_up = 12 * 2,
      event_fraction = 0.75,
      one_year_survival_rate = 0.1,
      concord = 0.95,
      median_surv = 3 # months
    ),
    sim8 = list(
      true_beta = c(log(1.95), log(0.001)),
      beta_hat = c(log(1.95), log(0.001)) * 1.25,
      lambda = 4e-4,
      gamma = 5.4,
      max_follow_up = 12 * 2,
      event_fraction = 0.67,
      one_year_survival_rate = 0.2,
      concord = 0.95,
      median_surv = 4 # months
    ),
    sim9 = list(
      true_beta = c(log(1.95), log(0.001)),
      beta_hat = c(log(1.95), log(0.001)) * 1.25,
      lambda = 4e-4,
      gamma = 3.1,
      max_follow_up = 12 * 2,
      event_fraction = 0.44,
      one_year_survival_rate = 0.5,
      concord = 0.95,
      median_surv = 12 # months
    )
  )

  return(simulation_settings)
}

#' Generate survival data using simsurv package
#'
gen_surv_data <- function(lambda,
                          gamma,
                          true_beta,
                          max_follow_up,
                          sample_size,
                          seed) {
  predictors <- data.frame(
    id = seq_len(sample_size),
    x1 = rnorm(sample_size),
    x2 = rnorm(sample_size)
  )
  df <- simsurv::simsurv(
    dist = "weibull",
    lambdas = c(lambda),
    gammas = c(gamma),
    x = predictors,
    betas = setNames(true_beta, c("x1", "x2")),
    maxt = .Machine$double.max - 1e-9,
    interval = c(0, .Machine$double.max),
    seed = seed
  )
  df <- merge(df, predictors)
  df$censorTime <- max_follow_up * runif(sample_size)
  df$survTime <- df$eventtime
  df$obsTime <- pmin(df$censorTime, df$survTime)
  df$status <- as.numeric(df$censorTime > df$survTime)
  df <- df[
    c(
      "id", "x1", "x2",
      "survTime", "censorTime",
      "obsTime", "status"
    )
  ]

  return(df)
}

#' Generate survival data using simstudy package
#'
gen_surv_data_old <- function(surv_formula,
                              surv_shape,
                              surv_scale,
                              censor_shape,
                              censor_scale,
                              sample_size) {
  import::from(magrittr, `%>%`)
  def_x <- simstudy::defData(
    varname = "x1",
    dist = "exponential",
    formula = 1
  ) %>%
    simstudy::defData(
      varname = "x2",
      dist = "exponential",
      formula = 1
    )
  def_surv <- simstudy::defSurv(
    varname = "survTime",
    formula = surv_formula,
    scale = surv_shape,
    shape = surv_scale
  ) %>%
    simstudy::defSurv(
      varname = "censorTime",
      shape = censor_shape,
      scale = censor_scale
    )

  surv_data <- simstudy::genData(
    n = sample_size,
    def_x
  ) %>%
    simstudy::genSurv(
      def_surv
      # timeName = "obsTime",
      # censorName = "censorTime",
      # eventName = "status",
      # keepEvents = TRUE
    ) %>%
    tibble::as_tibble()
  # tibble::as_tibble() %>%
  # dplyr::select(dplyr::contains("Time"), status, x1, x2)
  return(surv_data)
}


#' Simulate population data for DCA simulation (Results's subsection 02)
#'
simulate_dca_population_surv <- function(sim_setting, n_pop, thresholds, .seed,
                                         pred_time = 12, .verbose = FALSE, rm_zeros = FALSE) {
  thresholds <- validate_thresholds(thresholds = thresholds)
  msg <- cli::col_blue(
    paste0(
      "Simulating population DCA data with population seed = ", .seed
    )
  )
  message(msg)
  set.seed(.seed)
  df_pop <- gen_surv_data(
    lambda = sim_setting$lambda,
    gamma = sim_setting$gamma,
    true_beta = sim_setting$true_beta,
    max_follow_up = sim_setting$max_follow_up,
    sample_size = n_pop,
    seed = .seed
  )

  n_zeros_surv <- sum(df_pop$survTime == 0)
  n_zeros_cens <- sum(df_pop$censorTime == 0)
  n_zeros_obs <- sum(df_pop$obsTime == 0)
  if (n_zeros_obs > 0) {
    msg <- cli::col_br_red(
      paste0(
        "Zero-type issues detected: ",
        " surv ", n_zeros_surv,
        " cens ", n_zeros_cens,
        " obs ", n_zeros_obs
      )
    )
    message(msg)
  }
  if (n_zeros_obs > 10) {
    msg <- stringr::str_glue(
      "Number of zeroed observations cannot exceed 10 (got {n_zeros_obs})."
    )
    stop(msg)
  }
  if (isTRUE(rm_zeros)) {
    df_pop <- df_pop %>%
      dplyr::filter(obsTime > 0)
    n_pop <- nrow(df_pop)
  }

  if (isTRUE(.verbose)) {
    msg <- cli::col_blue(
      paste0(
        "\tTrue median surv (months): ", round(median(df_pop$survTime)),
        "\n\t1-year survival rate (%): ", round(mean(df_pop$survTime > 12) * 100),
        "\n\tSurvival rate at prediction time (%): ", round(mean(df_pop$survTime > pred_time) * 100)
      )
    )
    message(msg)
  }

  # placeholder coxph fit object
  ix <- seq_len(min(nrow(df_pop), 1e5))
  fit_hat <- survival::coxph(survival::Surv(obsTime, status) ~ x1 + x2, data = df_pop[ix, ])
  fit_hat$coef <- fit_hat$coefficients <- sim_setting$beta_hat

  df_pop$p_hat <- 1 - predict(
    fit_hat,
    newdata = df_pop %>% dplyr::mutate(obsTime = pred_time),
    type = "survival"
  )
  df_pop$lp_hat <- predict(
    fit_hat,
    newdata = df_pop %>% dplyr::mutate(obsTime = pred_time),
    type = "lp"
  )

  fit_true <- survival::coxph(survival::Surv(obsTime, status) ~ x1 + x2, data = df_pop[ix, ])
  fit_true$coef <- fit_true$coefficients <- sim_setting$true_beta

  df_pop$true_p <- 1 - predict(
    fit_true,
    newdata = df_pop %>% dplyr::mutate(obsTime = pred_time),
    type = "survival"
  )
  df_pop$true_lp <- predict(
    fit_true,
    newdata = df_pop %>% dplyr::mutate(obsTime = pred_time),
    type = "lp"
  )

  ## calculate true NB|y,p_hat
  true_nb <- purrr::map_df(thresholds, ~ {
    compute_nb_surv(
      surv_time = df_pop$survTime,
      p_hat = df_pop$p_hat,
      pred_time = pred_time,
      thr = .x
    )
  })
  performance <- get_performance_surv(
    surv_time = df_pop$survTime,
    obs_time = df_pop$obsTime,
    status = df_pop$status,
    lp = df_pop$lp_hat,
    predicted_risk = df_pop$p_hat,
    pred_time = pred_time
  )

  performance_true_beta <- get_performance_surv(
    surv_time = df_pop$survTime,
    obs_time = df_pop$obsTime,
    status = df_pop$status,
    lp = df_pop$true_lp,
    predicted_risk = df_pop$true_p,
    pred_time = pred_time
  )

  output <- list(
    df_pop = df_pop,
    true_nb = true_nb,
    thresholds = thresholds,
    n_pop = n_pop,
    true_beta = sim_setting$true_beta,
    beta_hat = sim_setting$beta_hat,
    sim_setting = sim_setting,
    n_zeros_surv = n_zeros_surv,
    n_zeros_cens = n_zeros_cens,
    calibration_oe = performance$oe,
    calibration_slope = performance$slope,
    concordance = performance$concordance,
    event_rate = mean(df_pop$survTime <= pred_time),
    average_prediction = mean(df_pop$p_hat),
    population_seed = .seed,
    performance = performance,
    performance_true_beta = performance_true_beta
  )

  return(output)
}


get_thr_label <- function(x) {
  labs <- character(length(x))
  for (i in seq_along(x)) {
    if (x[i] <= 1e-9) {
      labs[i] <- "0%"
    } else if (x[i] < 0.01) {
      labs[i] <- paste0(round(x[i] * 100, 2), "%")
    } else {
      labs[i] <- paste0(round(x[i] * 100), "%")
    }
  }
  return(labs)
}

#' EVPI for external validation
#' (adapts `predtools::evpi_val()` to allow a
#' "Full Bayes" method based on the BayesDCA R package.)
evpi_val2 <- function(Y,
                      pi,
                      method = c("bootstrap", "bayesian_bootstrap", "asymptotic", "full_bayes"),
                      n_sim = 1000,
                      threshold_varying_prior = FALSE,
                      max_sens_prior_mean = NULL,
                      min_sens_prior_mean = NULL,
                      max_sens_prior_sample_size = NULL,
                      prev_prior_mean = 0.5,
                      prev_prior_sample_size = 2,
                      prior_p = NULL,
                      zs = (0:99) / 100,
                      weights = NULL) {
  n <- length(Y)
  if (method == "asymptotic") {
    if (is.null(weights)) {
      weights <- rep(1, n)
    }
    ENB_perfect <- ENB_current <- rep(0, length(zs))
    for (j in seq_along(zs)) {
      NB_model <- sum(weights * (pi > zs[j]) * (Y - (1 - Y) * zs[j] / (1 - zs[j]))) / n
      NB_all <- sum(weights * (Y - (1 - Y) * zs[j] / (1 - zs[j]))) / n
      parms <- predtools::calc_NB_moments(Y, pi, zs[j], weights)
      if (is.na(parms[5])) {
        ENB_perfect[j] <- ENB_current[j] <- max(0, NB_model, NB_all)
      } else {
        if (parms[5] > 0.999999) {
          parms[5] <- 0.999999
        }
        if (parms[5] < -0.999999) {
          parms[5] <- -0.999999
        }
        tryCatch(
          {
            ENB_perfect[j] <- predtools::mu_max_trunc_bvn(
              parms[1],
              parms[2], parms[3], parms[4], parms[5]
            )
          },
          error = function(cond) {
            return(NULL)
          }
        )
        ENB_current[j] <- max(0, NB_model, NB_all)
      }
    }
    return(data.frame(
      z = zs, ENB_perfect = ENB_perfect,
      ENB_current = ENB_current, EVPIv = ENB_perfect - ENB_current
    ))
  }
  NB_model <- NB_all <- matrix(0, n_sim, ncol = length(zs))
  if (method == "bootstrap" || method == "bayesian_bootstrap") {
    Bayesian_bootstrap <- method == "bayesian_bootstrap"
    for (i in 1:n_sim) {
      w_x <- bootstrap(n, Bayesian_bootstrap, weights = weights)
      for (j in seq_along(zs)) {
        NB_model[i, j] <- sum(w_x * (pi > zs[j]) * (Y - (1 - Y) * zs[j] / (1 - zs[j]))) / n
        NB_all[i, j] <- sum(w_x * (Y - (1 - Y) * zs[j] / (1 - zs[j]))) / n
      }
    }
  } else if (method == "full_bayes") {
    priors <- bayesDCA:::.get_prior_parameters(
      thresholds = zs,
      threshold_varying_prior = threshold_varying_prior,
      prior_p = prior_p,
      ignorance_region_cutpoints = NULL,
      max_sens_prior_mean = max_sens_prior_mean,
      min_sens_prior_mean = min_sens_prior_mean,
      max_sens_prior_sample_size = max_sens_prior_sample_size,
      prev_prior_mean = prev_prior_mean,
      prev_prior_sample_size = prev_prior_sample_size,
      n_strategies = 1
    )
    prev <- rbeta(n_sim, shape1 = sum(Y) + priors$p1, shape2 = n - sum(Y) + priors$p2)
    for (j in seq_along(zs)) {
      tp <- sum(Y == 1L & pi > zs[j])
      fn <- sum(Y == 1L & pi <= zs[j])
      tn <- sum(Y == 0L & pi <= zs[j])
      fp <- sum(Y == 0L & pi > zs[j])
      se_j <- rbeta(n_sim, shape1 = tp + priors$Se1[j, 1], shape2 = fn + priors$Se2[j, 1])
      sp_j <- rbeta(n_sim, shape1 = tn + priors$Sp1[j, 1], shape2 = fp + priors$Sp2[j, 1])
      NB_model[, j] <- prev * se_j - (1 - prev) * (1 - sp_j) * zs[j] / (1 - zs[j])
      NB_all[, j] <- prev * 1 - (1 - prev) * (1 - 0) * zs[j] / (1 - zs[j])
    }
  } else {
    stop("Method ", method, " is not recognized.")
  }
  ENB_model <- ENB_all <- ENB_perfect <- ENB_current <- EVPIv <- p_useful <- rep(NA, length(zs))
  for (i in seq_along(zs)) {
    ENB_model[i] <- mean(NB_model[, i])
    ENB_all[i] <- mean(NB_all[, i])
    ENB_perfect[i] <- mean(pmax(
      NB_model[, i], NB_all[, i],
      0
    ))
    ENB_current[i] <- max(ENB_model[i], ENB_all[i], 0)
    EVPIv[i] <- ENB_perfect[i] - ENB_current[i]
    p_useful[i] <- mean((pmax(
      NB_model[, i], NB_all[, i],
      0
    ) - NB_model[, i]) == 0)
  }
  data.frame(
    z = zs, ENB_model = ENB_model, ENB_all = ENB_all,
    ENB_current = ENB_current, ENB_perfect = ENB_perfect,
    EVPIv = EVPIv, p_useful = p_useful
  )
}


#' Code adapted from Sadatsafavi et al (2023) [DOI: 10.1177/0272989X231178317]
#' by Giuliano Netto Flores Cruz on June 27, 2023.
#' Originally used for Figure 5 in the cited paper.
#' Original source code URL:
#' https://github.com/resplab/papercode/blob/main/voipred_ex/CaseStudy.Rmd
sim_by_size <- function(data_us,
                        model,
                        n_sim = 1000,
                        sample_sizes = c(250, 500, 1000, 2000, 4000, 8000, 16000, Inf),
                        zs = c(0.01, 0.02, 0.05, 0.1)) {
  set.seed(1)
  out <- data.frame(method = character(), sample_size = integer())
  for (i in seq_along(zs)) {
    out[paste0("val", i)] <- double()
  }

  index <- 1

  for (s in seq_along(sample_sizes)) {
    sample_size <- sample_sizes[s]
    cat("Sample size: ", sample_size, "\n")
    res_bb <- res_ob <- res_ll <- rep(0, length(zs))

    for (i in 1:n_sim) {
      if (is.infinite(sample_size)) {
        this_data <- data_us
        this_data$pi <- predict(model, newdata = this_data, type = "response")
        sample_size <- dim(this_data)[1]
      } else {
        repeat {
          this_data <- data_us[sample(1:(dim(data_us)[1]), sample_size, FALSE), ] # nolint
          this_data$pi <- predict(model, newdata = this_data, type = "response")
          if (min(this_data$pi) < min(zs)) {
            break
          } else {
            cat("bad external validation sample.")
          }
        }
      }

      bad <- FALSE
      tmp <- predtools::evpi_val(
        Y = this_data$Y, pi = this_data$pi,
        method = "bayesian_bootstrap", zs = zs
      )
      if (is.null(tmp)) bad <- TRUE
      out[index, "method"] <- "BB"
      out[index, "sample_size"] <- sample_size
      out[index, c("val1", "val2", "val3", "val4")] <- tmp$EVPIv
      index <- index + 1

      tmp <- predtools::evpi_val(
        Y = this_data$Y, pi = this_data$pi,
        method = "bootstrap", zs = zs
      )
      if (is.null(tmp)) bad <- TRUE
      out[index, "method"] <- "OB"
      out[index, "sample_size"] <- sample_size
      out[index, c("val1", "val2", "val3", "val4")] <- tmp$EVPIv
      index <- index + 1

      tmp <- predtools::evpi_val(
        Y = this_data$Y, pi = this_data$pi,
        method = "asymptotic", zs = zs
      )
      if (is.null(tmp)) bad <- TRUE
      out[index, "method"] <- "asy"
      out[index, "sample_size"] <- sample_size
      out[index, c("val1", "val2", "val3", "val4")] <- tmp$EVPIv
      index <- index + 1

      tmp <- evpi_val2(
        Y = this_data$Y, pi = this_data$pi,
        method = "full_bayes", zs = zs
      )
      if (is.null(tmp)) bad <- TRUE
      out[index, "method"] <- "full_bayes"
      out[index, "sample_size"] <- sample_size
      out[index, c("val1", "val2", "val3", "val4")] <- tmp$EVPIv
      index <- index + 1

      tmp <- evpi_val2(
        Y = this_data$Y, pi = this_data$pi,
        method = "full_bayes",
        zs = zs,
        threshold_varying_prior = TRUE,
        max_sens_prior_mean = 0.95,
        min_sens_prior_mean = 0.5,
        max_sens_prior_sample_size = 10,
        prev_prior_mean = 0.5,
        prev_prior_sample_size = 2
      )
      if (is.null(tmp)) bad <- TRUE
      out[index, "method"] <- "full_bayes_informative"
      out[index, "sample_size"] <- sample_size
      out[index, c("val1", "val2", "val3", "val4")] <- tmp$EVPIv
      index <- index + 1

      if (bad) {
        index <- index - 5
        i <- i - 1
        message("bad")
      }
    }
  }

  return(out)
}

run_sample_size_simulation_single_setting <- function(n_sim, n) {
  setting_results <- vector("list", length = n_sim)
  for (j in seq_len(n_sim)) {
    setting_results[[j]] <- run_sample_size_simulation_single_run(n = n)
    setting_results[[j]]$run_id <- paste0("run_", j)
  }
  return(dplyr::bind_rows(setting_results))
}

run_sample_size_simulation_single_run <- function(n) {
  # simulate data
  ## baseline model
  ### calib int -0.01
  ### calib slope 0.99
  ### concordance 0.69
  ## new model
  ### calib int -0.15
  ### calib slope 0.86
  ### concordance 0.74
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  true_lp <- -2 + 0.9 * x1 + 0.7 * x2 + rnorm(n)
  lp0 <- -1.6 + 0.74 * x1
  lp1 <- -1.8 + 0.8 * x1 + 0.8 * x2
  y <- rbinom(n, size = 1, prob = plogis(true_lp))
  d <- data.frame(outcomes = y, phat0 = plogis(lp0), phat1 = plogis(lp1), true_p = plogis(true_lp))
  f <- bayesDCA::dca(d, thresholds = seq(0, 0.3, 0.01))
  # compute run results for probabilities
  pbest <- bayesDCA::plot_superiority_prob(f)$data %>%
    dplyr::mutate(comparison = paste0(decision_strategy, " vs best competitor")) %>%
    dplyr::select(threshold, comparison, prob)
  puseful <- bayesDCA::plot_superiority_prob(f, type = "useful")$data %>%
    dplyr::mutate(comparison = paste0(decision_strategy, " vs max(treat all, treat none)")) %>%
    dplyr::select(threshold, comparison, prob)
  phat1_vs_phat0 <- bayesDCA::plot_superiority_prob(f, type = "pairwise", strategies = c("phat1", "phat0"))$data %>%
    dplyr::mutate(comparison = paste0("phat1 vs phat0")) %>%
    dplyr::select(threshold, comparison, prob)
  phat1_vs_true_p <- bayesDCA::plot_superiority_prob(f, type = "pairwise", strategies = c("phat1", "true_p"))$data %>%
    dplyr::mutate(comparison = paste0("phat1 vs true p")) %>%
    dplyr::select(threshold, comparison, prob)
  phat0_vs_true_p <- bayesDCA::plot_superiority_prob(f, type = "pairwise", strategies = c("phat0", "true_p"))$data %>%
    dplyr::mutate(comparison = paste0("phat0 vs true p")) %>%
    dplyr::select(threshold, comparison, prob)
  # compute run results for delta estimates
  delta_best <- bayesDCA::plot_delta(f)$data %>%
    dplyr::mutate(comparison = paste0(decision_strategy, " vs best competitor")) %>%
    dplyr::select(threshold, comparison, delta := estimate)
  delta_useful <- bayesDCA::plot_delta(f, type = "useful")$data %>%
    dplyr::mutate(comparison = paste0(decision_strategy, " vs max(treat all, treat none)")) %>%
    dplyr::select(threshold, comparison, delta := estimate)
  delta_phat1_vs_phat0 <- bayesDCA::plot_delta(f, type = "pairwise", strategies = c("phat1", "phat0"))$data %>%
    dplyr::mutate(comparison = paste0("phat1 vs phat0")) %>%
    dplyr::select(threshold, comparison, delta := estimate)
  delta_phat1_vs_true_p <- bayesDCA::plot_delta(f, type = "pairwise", strategies = c("phat1", "true_p"))$data %>%
    dplyr::mutate(comparison = paste0("phat1 vs true p")) %>%
    dplyr::select(threshold, comparison, delta := estimate)
  delta_phat0_vs_true_p <- bayesDCA::plot_delta(f, type = "pairwise", strategies = c("phat0", "true_p"))$data %>%
    dplyr::mutate(comparison = paste0("phat0 vs true p")) %>%
    dplyr::select(threshold, comparison, delta := estimate)

  # combine results
  results_prob <- dplyr::bind_rows(
    pbest, puseful, phat1_vs_phat0, phat1_vs_true_p, phat0_vs_true_p
  )
  results_delta <- dplyr::bind_rows(
    delta_best, delta_useful, delta_phat1_vs_phat0, delta_phat1_vs_true_p, delta_phat0_vs_true_p
  )
  run_results <- dplyr::inner_join(
    results_prob, results_delta,
    by = c("threshold", "comparison")
  )
  return(run_results)
}
