library(tidyverse)
library(bayesDCA)

d <- read_csv(
    "data/12916_2019_1425_MOESM1_ESM.csv",
    show_col_types = FALSE
) %>%
    dplyr::select(
        outcomes := outcome,
        model_predictions := pred
    )

# noise predictor
set.seed(12112022)
z <- rnorm(nrow(d), mean = 0, sd = 2)

d <- d %>%
    dplyr::mutate(
        binary_test = as.numeric(model_predictions > 0.5),
        model_predictions2 = plogis(
            qlogis(model_predictions) + log(5) * z
        )
    )

fit <- dca(d, cores = 4, refresh = 1)

.labels <- list(
    "model_predictions" = "ADNEX",
    "binary_test" = "Binary test/Classifier",
    "model_predictions2" = "ADNEX updated"
)
plot(fit, labels = .labels) +
    ggplot2::theme(
        legend.position = c(0.2, 0.35),
        legend.text = ggplot2::element_text(size = 9)
    )

ggsave(
    "output/case-study/dca.png",
    width = 8, height = 4.5, dpi = 600
)

plots_best <- compare_dca(fit, plot_list = TRUE)
plots_useful <- compare_dca(fit, plot_list = TRUE, type = "useful")
