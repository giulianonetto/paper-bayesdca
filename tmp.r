B <- 2e3
d <- 50
N <- 100
p <- rbeta(B, 0.5 + d, 0.5 + N - d)
ts <- seq(0, 0.99, 0.01)
odds <- \(x) x / (1 - x)
get_se_mu <- \(t, shift = 0.45, slope = 0.025) {
  plogis(shift + slope * qlogis(1 - t)^3)
}

{
  plot(ts, get_se_mu(ts), ylab = "Prior Se mean", ylim = c(0, 1), type = "l", lwd = 5)
  axis(1, tck = 1, col.ticks = "light gray")
  axis(1, tck = -0.015, col.ticks = "black")
  axis(2, tck = 1, col.ticks = "light gray", lwd.ticks = "1")
  axis(2, tck = -0.015)
  Hmisc::minor.tick(nx = 5, ny = 2, tick.ratio = 0.5)
  par(new = TRUE)
  plot(ts, get_se_mu(ts), ylab = "Prior Se mean", ylim = c(0, 1), type = "l", lwd = 5)
  abline(h = c(.05, .1, .6, .95), lty = 2, col = "red")
}

get_sp_mu <- \(t, shift = 0.45, slope = 0.02) {
  plogis(shift + slope * qlogis(t)^3)
}

{
  plot(ts, get_sp_mu(ts), ylab = "Prior Sp mean", ylim = c(0, 1), type = "l", lwd = 5)
  axis(1, tck = 1, col.ticks = "light gray")
  axis(1, tck = -0.015, col.ticks = "black")
  axis(2, tck = 1, col.ticks = "light gray", lwd.ticks = "1")
  axis(2, tck = -0.015)
  Hmisc::minor.tick(nx = 5, ny = 2, tick.ratio = 0.5)
  par(new = TRUE)
  plot(ts, get_sp_mu(ts), ylab = "Prior Sp mean", ylim = c(0, 1), type = "l", lwd = 5)
  abline(h = c(.05, .1, .6, .95), lty = 2, col = "red")
}
nb <- sapply(ts, \(t) {
  cat(paste0("t=", t, "\n"))
  smpl_size <- 2
  se_mu <- get_se_mu(t)
  sp_mu <- get_sp_mu(t)
  # se_mu <- (1 - t)*0.7
  # sp_mu <- t*0.7
  cat(paste0("t=", t, " Se mu:", se_mu, "\tSp mu: ", sp_mu, "\n"))
  alpha_se <- se_mu * smpl_size
  beta_se <- (1 - se_mu) * smpl_size
  # cat(paste0("\tSe pars: ", alpha_se, "  ", beta_se, "\n"))
  se <- rbeta(B, alpha_se, beta_se)
  alpha_sp <- sp_mu * smpl_size
  beta_sp <- (1 - sp_mu) * smpl_size
  # cat(paste0("\tSp pars: ", alpha_sp, "  ", beta_sp, "\n"))
  sp <- rbeta(B, alpha_sp, beta_sp)
  se * p - (1 - sp) * (1 - p) * (t / (1 - t))
})

par(mfrow = c(1, 2))
plot(ts, colMeans(nb), type = "l", ylim = c(-1, 1), lwd = 5)
abline(h = 0, lwd = 2, col = "red")
points(ts, apply(nb, 2, quantile, probs = 0.975), type = "l")
points(ts, apply(nb, 2, quantile, probs = 0.025), type = "l")
points(
  ts, get_se_mu(ts) * d / N - (1 - d / N) * (1 - get_sp_mu(ts)) * odds(ts),
  type = "l", col = "red", lty = 2,
  lwd = 2,
  main = "Better prior"
)

# uniform prior


nb <- sapply(ts, \(t) {
  cat(paste0("t=", t, "\n"))
  smpl_size <- 2
  se_mu <- 1 / 2
  sp_mu <- 1 / 2
  # se_mu <- (1 - t)*0.7
  # sp_mu <- t*0.7
  cat(paste0("t=", t, " Se mu:", se_mu, "\tSp mu: ", sp_mu, "\n"))
  alpha_se <- se_mu * smpl_size
  beta_se <- (1 - se_mu) * smpl_size
  # cat(paste0("\tSe pars: ", alpha_se, "  ", beta_se, "\n"))
  se <- rbeta(B, alpha_se, beta_se)
  alpha_sp <- sp_mu * smpl_size
  beta_sp <- (1 - sp_mu) * smpl_size
  # cat(paste0("\tSp pars: ", alpha_sp, "  ", beta_sp, "\n"))
  sp <- rbeta(B, alpha_sp, beta_sp)
  se * p - (1 - sp) * (1 - p) * (t / (1 - t))
})

plot(
  ts, colMeans(nb),
  type = "l", ylim = c(-1, 1), lwd = 5,
  main = "Uniform prior"
)
abline(h = 0, lwd = 2, col = "red")
points(ts, apply(nb, 2, quantile, probs = 0.975), type = "l")
points(ts, apply(nb, 2, quantile, probs = 0.025), type = "l")
points(
  ts, 1 / 2 * d / N - (1 - d / N) * (1 - 1 / 2) * odds(ts),
  type = "l", col = "red", lty = 2,
  lwd = 2
)
