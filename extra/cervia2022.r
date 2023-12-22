
# Code for manuscript entitled: "Immunoglobulin signature predicts risk of post-acute COVID-19 syndrome", Cervia et al., Nature Communications, 2021

# Legend:

# 1. Import packages and data
# 2. Data wrangling and grouping
# 3. Functions for data analysis
# 4. Analysis of total immunoglobulin levels
#   a) IgM
#   b) IgG1
#   c) IgG3
#   d) IgA
#   e) IgG2
#   f) IgG4
# 5. Development of prediction model - Variables
#   a) Linearity of variables
#   b) Age
#   c) Number of symptoms during primary infection
#   d) History of asthma bronchiale
#   e) COVID-19 severity
# 6. Development of prediction model - PACS score
#   a) Primary infection
#   b) 6-months follow-up
#   c) Decision curve analysis
#   d) Interaction plots
# 7. Code for implementation of PACS score




library("tidyverse") # version 1.3.1
library("dcurves") # version 0.2.0


# Data wrangling and grouping----

# Import Data

xdata <- readxl::read_excel(here::here("data/41467_2021_27797_MOESM5_ESM.xlsx"), na = "NA") # Enter excel file name

# Factors

xdata$Sampling_month <- factor(xdata$Sampling_month,
    levels = c("ONE", "SIX", "TWELVE")
)

xdata$COVID <- factor(xdata$COVID,
    levels = c("Healthy", "COVID")
)

xdata$SevMax2 <- factor(xdata$SevMax2,
    levels = c("Healthy", "Mild", "Severe")
)

xdata$SevMax6 <- factor(xdata$SevMax6,
    levels = c("Healthy", "Asymptomatic", "Mild", "Pneu", "Sev Pneu", "Mild ARDS", "Mod ARDS", "Sev ARDS")
)

xdata$Hospitalized <- factor(xdata$Hospitalized,
    levels = c("NO", "YES")
)

xdata$PACS <- factor(xdata$PACS,
    levels = c("NO", "YES")
)


# Add Healhty  to six month measurements
healthy <- xdata %>%
    filter(SevMax2 == "Healthy") %>%
    mutate(Sampling_month = "SIX", Complete_OST = "NO")
xdata <- rbind(xdata, healthy)



# Create groups
xdata_selectR <- xdata[xdata$Resampling %in% c("YES"), ] # followed-up patients
xdata_selectR2 <- xdata[xdata$Resampling_andH %in% c("YES"), ] # followed-up patients + healthy controls

xdata_selectONE <- xdata[xdata$Sampling_month %in% c("ONE"), ] # Primary infection
xdata_selectRONE <- xdata_selectR[xdata_selectR$Sampling_month %in% c("ONE"), ]
xdata_selectR2ONE <- xdata_selectR2[xdata_selectR2$Sampling_month %in% c("ONE"), ]

xdata_selectRSIX <- xdata_selectR[xdata_selectR$Sampling_month %in% c("SIX"), ] # 6-month follow-up
xdata_selectR2SIX <- xdata_selectR2[xdata_selectR2$Sampling_month %in% c("SIX"), ]

xdata_selectR2OS <- xdata_selectR2[xdata_selectR2$Sampling_month %in% c("ONE", "SIX"), ] # Primary infection & 6-month follow-up

xdata_selectRnH <- xdata_selectR[xdata_selectR$SevMax2 %in% c("Mild", "Severe"), ] # Without Healthy controls
xdata_selectRONEnH <- xdata_selectRONE[xdata_selectRONE$SevMax2 %in% c("Mild", "Severe"), ]

# PACS score----

PACS_model <- glm(PACS ~ scale(Age) + Nr_Symptoms + Asthma + IgM * IgG3,
    xdata_selectRONEnH,
    family = binomial(link = "logit"), x = TRUE
)

xdata_selectRONEnH$PC1 <- xdata_selectRONEnH_2 %>%
    select(contains("Ig")) %>%
    prcomp() %>%
    pluck("x") %>%
    data.frame() %>%
    select(PC1) %>%
    unlist()
PACS_model <- glm(PACS ~ scale(Age) + Nr_Symptoms + Asthma + IgM * IgG3,
    xdata_selectRONEnH,
    family = binomial(link = "logit"), x = TRUE
)
PACS_model2 <- glm(PACS ~ rms::rcs(Age, 3) + Nr_Symptoms + Asthma,
    xdata_selectRONEnH,
    family = binomial(link = "logit"), x = TRUE
)

PACS_predictions <- predict(PACS_model, xdata_selectRONEnH, type = "response")
per_roc_LC1_Asthma_IgG3_IgM_IgA1_Sym <- pROC::roc(PACS ~ PACS_predictions, xdata_selectRONEnH, ci = TRUE)


# __Primary infection----
PACS_model_sh <- shrink::shrink(PACS_model, type = "global", method = "dfbeta")

# __Decision curve analysis----

xdata_selectRONEnH_2 <-
    xdata_selectRONEnH %>%
    mutate(
        prob_PACS = predict(PACS_model, xdata_selectRONEnH, type = "response"),
        prob_PACS2 = predict(PACS_model2, xdata_selectRONEnH, type = "response"),
        prob_PACS_shrunk = predict(PACS_model_sh, xdata_selectRONEnH, type = "response")
    )

dcurves::dca(
    PACS ~
        prob_PACS,
    xdata_selectRONEnH_2
) %>%
    plot(smooth = TRUE)

# bayesDCA ----

bdca <- xdata_selectRONEnH_2 %>%
    mutate(outcomes = as.numeric(PACS == "YES")) %>%
    select(outcomes, prob_PACS) %>%
    bayesDCA::dca(
        thresholds = c(1e-10, 0.01, seq(0.02, 1, .02)),
        cores = 4,
        refresh = 1
    )

bdca2 <- xdata_selectRONEnH_2 %>%
    mutate(outcomes = as.numeric(PACS == "YES")) %>%
    select(outcomes, prob_PACS) %>%
    bayesDCA::dca(
        thresholds = c(1e-10, 0.01, seq(0.02, 1, .02)),
        cores = 4,
        refresh = 1,
        external_prevalence_data = c(62, 160)
    )
compare_dca(bdca2)

bdca3 <- xdata_selectRONEnH_2 %>%
    mutate(outcomes = as.numeric(PACS == "YES")) %>%
    select(outcomes, prob_PACS) %>%
    bayesDCA::dca(
        thresholds = c(1e-10, 0.01, seq(0.02, 1, .02)),
        cores = 4,
        refresh = 1,
        external_prevalence_data = c(216, 395)
    )
compare_dca(bdca3)
