library(tidyverse)
theme_set(theme_bw(base_size = 14))
.colors <- RColorBrewer::brewer.pal(3, "Dark2")
names(.colors) <- c(
    "True NB", "Frequentist", "Bayesian"
)
.colors[".true_nb"] <- .colors[1]

setting_labels_pretty <- map_chr(
    get_simulation_settings(),
    ~ paste0(
        "AUC ", .x$auc, ", prevalence ", round(.x$prev * 100), "%"
    )
)

df <- results_02_subsection %>%
    filter(strategy == "Model-based decisions") %>%
    mutate(
        setting_label = factor(
            setting_labels_pretty[setting_label],
            levels = setting_labels_pretty
        )
    ) %>%
    as_tibble()

# point estimates are nearly identical
df %>%
    select(threshold, setting_label, simulation_label, .type, estimate, .true_nb) %>%
    pivot_wider(names_from = .type, values_from = estimate) %>%
    pivot_longer(cols = c(.true_nb, Frequentist, Bayesian)) %>%
    mutate(
        name = ifelse(
            name == ".true_nb",
            "True NB",
            name
        ),
        name = fct_relevel(
            name,
            "True NB", "Bayesian", "Frequentist"
        )
    ) %>%
    ggplot(aes(factor(threshold), value)) +
    geom_boxplot(aes(color = name), position = position_dodge(width = .8), width = .5, lwd = 2) +
    facet_wrap(~setting_label, scales = "free") +
    scale_color_manual(values = .colors) +
    labs(
        x = "Decision threshold",
        y = "Net benefit",
        color = NULL
    )

# 95% intervals coverage

df %>%
    group_by(threshold, .type, setting_label) %>%
    summarise(
        cov = mean(truth_within_interval)
    ) %>%
    ggplot(aes(threshold, cov, color = .type)) +
    geom_point(
        # position = position_dodge(width = .025),
        aes(shape = .type, size = .type), stroke = 1.5
    ) +
    geom_line(
        aes(linetype = .type)
    ) +
    facet_wrap(~setting_label) +
    scale_y_continuous(
        labels = scales::label_percent(),
        limits = c(0.8, 1)
    ) +
    scale_color_manual(values = .colors) +
    scale_shape_manual(values = c(21, 19)) +
    scale_size_manual(values = c(4, 2)) +
    labs(
        x = "Decision threshold",
        y = "Empirical covarage\n(95% uncertainty intervals)",
        color = NULL,
        size = NULL,
        linetype = NULL,
        shape = NULL
    )
