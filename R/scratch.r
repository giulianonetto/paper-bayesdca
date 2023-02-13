
N <- 5e4
x1 <- rnorm(N, 0.046, sd = 0.012) * 3
x2 <- rnorm(N, 0.12, sd = 0.0061)
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
hist(d1, xlab = "Net Benefit differences against corresponding SoC", main = NULL, cex.lab = 1.5,probability=T)
hist(d2, col = "red", add = TRUE, probability=T)
legend(
    x = c(-.09, -.09),
    y = c(200000, 200000),
    col = c("gray80", "red"),
    pch = c(19, 19),
    legend = c(lab1, lab2),
    cex = 1.2, y.intersp = 2
)
dev.off()
library(tidyverse)
library(ggdist)
data.frame(
    d = c(d1, d2),
    l = rep(
        c(lab1, lab2),
        each = N
    ),
    f = rep(
        c(lab1, lab2),
        each = N
    ) %>% str_extract("Thromboembolism model|Sepsis model")
) %>%
    ggplot(aes(x = d, fill = l)) +
    geom_density() +
    facet_grid(rows = vars(f), scales = "free_y") +
    scale_x_continuous(breaks = scales::pretty_breaks(5)) +
    theme_bw(base_size = 16) +
    theme(
        legend.spacing.y = unit(1, "cm"),
        # legend.position = c(0.8, 0.),
        legend.background = element_rect(fill = "white", colour = NA)
    ) +
    guides(fill = guide_legend(byrow = TRUE)) +
    labs(
        x = "Net Benefit differences against corresponding SoC",
        y = "Posterior Density",
        fill = NULL
    ) +
    scale_fill_brewer(palette = "Dark2")
ggsave("example.png", dpi = 300, width = 10, height = 6)
