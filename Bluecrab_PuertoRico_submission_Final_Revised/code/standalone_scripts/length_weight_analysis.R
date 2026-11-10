# Weight-length relationships of C. sapidus and C. bocourti (Section 3.6; Figs. 8-9)
#
# Cleaned-up version of the weight-length models in the original working script
# (codeforpublication_1stpaper.R). Input: data/analysis_ready/weight_length_crabs.csv,
# the size records exactly as analysed, with the original maturity and ovigerous codes
# and the crab_id used in data/clean/crabdata_clean_1.csv.
#
# Filters used throughout:
#   - non-ovigerous crabs only (ovigerous_original == 0; blank records are dropped)
#   - immature crabs are dropped. Crabs with no maturity stage recorded are kept
#     (they are adult-sized and were mature crabs left unlabelled), except those
#     under 75 mm CW, the cut-off used in the size analysis (one 43-mm C. sapidus).
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
         !maturity_original %in% "IM",                       # drop immature
         !is.na(maturity_original) | cw_spines_mm >= 75) %>%  # unlabelled: keep if >= 75 mm
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
#    Paper: C. sapidus F(1,74) = 43.37; C. bocourti F(1,111) = 66.54 (logL rows)
best_pred <- logW ~ logCWspines + logCWbase + logL + logL:logCWspines
Anova(lm(best_pred, data = sap))
Anova(lm(best_pred, data = boc))

# 2. Does the length-weight relationship differ by sex? (sex rows)
#    Paper: C. sapidus F(1,62) = 6.06, p = 0.02; C. bocourti F(1,98) = 8.96, p = 0.003
sex_model <- logW ~ logL * sex + logL * station + logL * survey
Anova(lm(sex_model, data = sap))
Anova(lm(sex_model, data = boc))

# 3. Length-weight relationship of mature non-ovigerous female C. sapidus (Fig. 8A)
#    Paper: a = 0.0060, b = 2.40, R2 = 0.73, n = 68
sap_f <- filter(sap, sex == "F")
lw_sap_f <- lm(logW ~ logL, data = sap_f)
power_params(lw_sap_f)

#    Within each sex: does the length-weight relationship differ among
#    surveys or stations? (intercepts: survey, station; slopes: logL:survey,
#    logL:station). Too few C. sapidus males (n = 12) to fit this model.
boc_m <- filter(boc, sex == "M")
Anova(lm(logW ~ logL * survey + logL * station, data = sap_f))
Anova(lm(logW ~ logL * survey + logL * station, data = filter(boc, sex == "F")))
Anova(lm(logW ~ logL * survey + logL * station, data = boc_m))

# 4. Carapace width-weight relationship of mature non-ovigerous females (Fig. 9)
#    Paper: C. sapidus logCW x survey F(4,57) = 0.35, p = 0.84; a = 0.00018, b = 2.75
#           C. bocourti logCW x survey F(4,81) = 3.85, p = 0.006
boc_f <- filter(boc, sex == "F")
Anova(lm(logW ~ logCWspines * survey, data = sap_f))
power_params(lm(logW ~ logCWspines, data = sap_f))
Anova(lm(logW ~ logCWspines * survey, data = boc_f))
