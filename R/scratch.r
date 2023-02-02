
x1 <- rnorm(1e6, 0.046, sd = 0.012) * 3
x2 <- rnorm(1e6, 0.12, sd = 0.0061)
d1 <- x1 - .02
d2 <- x2 - .01
m1 <- round(mean(d1), 2)
m2 <- round(mean(d2), 2)
s1 <- round(quantile(d1, c(.025, .975)), 2)
s2 <- round(quantile(d2, c(.025, .975)), 2)
lab1 <- stringr::str_glue(
    "Sepsis model\n{m1} (95% Cr.I. {s1[1]}\u2012{s1[2]})"
)
lab2 <- stringr::str_glue(
    "Thromboembolism model\n{m2} (95% Cr.I. {s2[1]}\u2012{s2[2]})"
)
png("example.png", res = 300, width = 3600 * .8, height = 2400 * .8)
hist(d1, xlab = "Net Benefit differences against corresponding SoC", main = NULL, cex.lab = 1.5)
hist(d2, col = "red", add = TRUE)
legend(
    x = c(-.09, -.09),
    y = c(200000, 200000),
    col = c("gray80", "red"),
    pch = c(19, 19),
    legend = c(lab1, lab2),
    cex = 1.2, y.intersp = 2
)
dev.off()
