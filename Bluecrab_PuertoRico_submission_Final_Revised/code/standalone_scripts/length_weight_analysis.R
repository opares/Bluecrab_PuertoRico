# Weight-length relationships of C. sapidus and C. bocourti (Section 3.6; Figs. 8-9)
#
# Cleaned-up version of the weight-length models in the original working script
# (codeforpublication_1stpaper.R). Input: data/analysis_ready/weight_length_crabs.csv,
# the size records exactly as analysed, with the original maturity and ovigerous codes
# and the crab_id used in data/clean/crabdata_clean_1.csv.
#
# Filters used throughout:
#   - non-ovigerous crabs only (ovigerous_original == 0; blank records are dropped)
#   - immature crabs and crabs with no maturity record are dropped
#   - CR306 (C. sapidus female, 14 Oct 2021, AP): recorded spine width of 160.89 mm is
#     implausible for a crab with 52.2 mm carapace length and 76.1 g weight (base width
#     88.15 mm), so its spine width is set to NA. Its length and weight are kept.
# Weight and size are log10-transformed; tests are Type II ANOVA (car::Anova).
#
# Run from the Bluecrab_PuertoRico_submission_Final_Revised folder.

library(dplyr)
library(car)

wl <- read.csv("data/analysis_ready/weight_length_crabs.csv", na.strings = c("", "NA"))

wl$cw_spines_mm[wl$crab_id %in% "CR306"] <- NA

wl <- wl %>%
  filter(ovigerous_original == 0,
         !is.na(maturity_original), maturity_original != "IM") %>%
  mutate(logW        = log10(weight_g),
         logL        = log10(carapace_length_mm),
         logCWspines = log10(cw_spines_mm),
         logCWbase   = log10(cw_base_mm),
         sex     = factor(sex),
         station = factor(station),
         survey  = factor(survey))

sap <- filter(wl, species == "Callinectes sapidus")
boc <- filter(wl, species == "Callinectes bocourti")

# Power curve W = a * X^b from a log10-log10 fit
power_params <- function(fit) {
  c(a = unname(10^coef(fit)[1]), b = unname(coef(fit)[2]),
    R2 = summary(fit)$r.squared, n = nobs(fit))
}

# 1. Which size measurement best predicts weight? (both sexes)
#    Paper: C. sapidus F(1,68) = 42.24; C. bocourti F(1,104) = 57.03 (logL rows)
best_pred <- logW ~ logCWspines + logCWbase + logL + logL:logCWspines
Anova(lm(best_pred, data = sap))
Anova(lm(best_pred, data = boc))

# 2. Does the length-weight relationship differ by sex? (sex rows)
#    Paper: C. sapidus F(1,56) = 6.07, p = 0.02; C. bocourti F(1,91) = 8.11, p = 0.005
sex_model <- logW ~ logL * sex + logL * station + logL * survey
Anova(lm(sex_model, data = sap))
Anova(lm(sex_model, data = boc))

# 3. Length-weight relationship of mature non-ovigerous female C. sapidus (Fig. 8A)
#    Paper: a = 0.0063, b = 2.39, R2 = 0.74
sap_f <- filter(sap, sex == "F")
lw_sap_f <- lm(logW ~ logL, data = sap_f)
power_params(lw_sap_f)

#    Stable across surveys? (survey and logL:survey rows not significant)
Anova(lm(logW ~ logL * survey, data = sap_f))

# 4. Carapace width-weight relationship of mature non-ovigerous females (Fig. 9)
#    Paper: C. sapidus logCW x survey F(4,51) = 0.36, p = 0.83; a = 0.00020, b = 2.73
#           C. bocourti logCW x survey F(4,77) = 3.72, p = 0.01
boc_f <- filter(boc, sex == "F")
Anova(lm(logW ~ logCWspines * survey, data = sap_f))
power_params(lm(logW ~ logCWspines, data = sap_f))
Anova(lm(logW ~ logCWspines * survey, data = boc_f))
