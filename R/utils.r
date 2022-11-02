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

#' Plot bayesDCA vs rmda
#'
#' Plots DCA from `bayesDCa` and `rmda` packages
#' for visual comparison
#'
#' @param dataset
#' @param outcomes outcome variable (character string with column name)
#' @param predictor predicted probabilities (character string with column name)
#' @param bootstraps Number (int) of bootstrap samples for rmda
#' @param treat_all_rmda Logical indicating wether to plot Treat all from rmda (defaults to FALSE).
#' @param refresh Refresh value for `rstan::sampling` (defaults to 0).
#' @param cores Number of cores for `bayesDCA::dca`
#' @thresholds Numeric vector (between 0 and 1) of thresholds for DCA.
#' @importFrom magrittr %>%
plot_bdca_vs_rmda <- function(dataset, outcomes, predictor, thresholds,
                              treat_all_rmda = FALSE, bootstraps = 500,
                              refresh = 0, cores = 4) {
  import::from(magrittr, `%>%`)
  df <- data.frame(
    outcomes = dataset[[outcomes]],
    predictor = dataset[[predictor]]
  )
  thresholds <- pmin(thresholds, 0.999) %>%
    pmax(0.001)
  # Estimate decision curves
  msg <- cli::col_blue("Estimating DCA with bayesDCA")
  message(msg)
  bdca_fit <- bayesDCA::dca(
    df,
    thresholds = thresholds,
    refresh = refresh,
    cores = cores
  )
  msg <- cli::col_blue("Estimating DCA with rmda")
  message(msg)
  rmda_fit <- rmda::decision_curve(
    outcomes ~ predictor,
    data = df,
    fitted.risk = TRUE,
    thresholds = thresholds,
    bootstraps = bootstraps
  )

  # get results into standardized data.frames
  msg <- cli::col_blue("Plotting results")
  message(msg)
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

  # get plot helper objects
  max_estimate <- max(res_all$estimate)
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
    res_all <- res_all %>%
      dplyr::filter(
        !(strategy == "Treat all" & .type == "Frequentist")
      )
  }

  .plot <- res_all %>%
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

#' Get subsamples
#' @param n_event Number of events
#' @param n_nonevent Number of non-events
#' @df Data to subsample
#' @y Outcome variable name
get_subsamples <- function(n_event, n_nonevent, df, y) {
  df[
    c(
      sample(which(df[[y]] == 1), n_event),
      sample(which(df[[y]] == 0), n_nonevent)
    ),
  ]
}
