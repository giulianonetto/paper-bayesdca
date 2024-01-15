{
  x <- simulate_dca_population_surv(
          sim_setting = get_simulation_settings_surv()$sim1,
          n_pop = 1e3,
          thresholds = c(0.1, 0.5, 0.75),
          .seed = 547,
          .verbose = TRUE
  )
  print(x$sim_setting$concord)
  cat("\n\n")
  print(data.frame(x$performance))
  cat("\n\n")
  print(data.frame(x$performance_true_beta))
}

