import::from(pROC, auc)
n <- 5e5
x <- cbind(1, rnorm(n), rnorm(n))
true_beta <- c(-5.65, -log(2.95), log(2.95))
beta_hat <- c(-7, -log(2.95) * 1.5, log(2.95) * 1.5)
p <- as.vector(plogis(x %*% true_beta))
phat <- as.vector(plogis(x %*% beta_hat))
y <- rbinom(n, 1, p)


test_sim_setting(get_simulation_settings()$sim1)
