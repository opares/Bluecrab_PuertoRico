# =============================================================================
# Carapace width of C. sapidus and C. bocourti: does size differ between
# sexes, among stations, and among quarterly surveys?
#
# Design: 4 stations x 5 surveys (20 station-survey combinations), crabs
# trapped over 2-3 consecutive days per combination. Crabs from the same
# combination are not independent, so a linear mixed model (LMM) is used with
# a random intercept for station-survey combination (the u(q,s) term of Eq. 3).
#
# Model (one per species):
#     CW ~ sex + station + survey + (1 | station:survey)
# Sex is included because the sex ratio differs among stations (Section 3.4)
# and the sexes differ in size, so station effects must be adjusted for sex.
# =============================================================================

library(dplyr)
library(lmerTest)     # lmer() with Satterthwaite F tests
library(emmeans)
library(DHARMa)
library(ggplot2)

# --- 1. Data -----------------------------------------------------------------
crabs <- read.csv("../../data/clean/crabdata_clean_1.csv") %>%
  filter(species %in% c("Callinectes sapidus", "Callinectes bocourti"),
         sex %in% c("F", "M"), !is.na(cw_spines_mm)) %>%
  mutate(species = factor(species, levels = c("Callinectes sapidus", "Callinectes bocourti")),
         station = factor(station, levels = c("BC", "AP", "DY", "IM")),
         survey  = factor(survey),
         sex     = factor(sex),
         ss      = interaction(station, survey, drop = TRUE))

# --- 2. Exploratory ----------------------------------------------------------
# 2a. Summary by species and sex
crabs %>% group_by(species, sex) %>%
  summarise(n = n(), mean = mean(cw_spines_mm), sd = sd(cw_spines_mm),
            min = min(cw_spines_mm), max = max(cw_spines_mm), .groups = "drop") %>%
  print(digits = 4)

# 2b. Sample sizes per station x survey (shows the imbalance)
for (sp in levels(crabs$species)) { cat("\n", sp, "\n"); print(with(filter(crabs, species == sp), table(station, survey))) }

# 2c. Small crabs: the lower tail is sparse (trap selectivity) and includes
#     the few immature individuals
print(crabs %>% filter(cw_spines_mm < 80) %>% select(species, crab_id, sex, station, survey, cw_spines_mm, maturity_coded))

# 2d. Plot: CW by station and by survey, by sex
ggplot(crabs, aes(station, cw_spines_mm, fill = sex)) + geom_boxplot(outlier.size = 0.8) +
  facet_wrap(~ species) + labs(y = "Carapace width (mm)") + theme_bw()
ggsave("fig_cw_by_station.png", width = 8, height = 4, dpi = 300)
ggplot(crabs, aes(survey, cw_spines_mm, fill = sex)) + geom_boxplot(outlier.size = 0.8) +
  facet_wrap(~ species) + labs(y = "Carapace width (mm)") + theme_bw()
ggsave("fig_cw_by_survey.png", width = 8, height = 4, dpi = 300)

# --- 3. Fit the LMM on all crabs and check assumptions ----------------------
for (sp in levels(crabs$species)) {
  dat <- filter(crabs, species == sp)
  m   <- lmer(cw_spines_mm ~ sex + station + survey + (1 | ss), data = dat)
  cat("\n=====", sp, "(all crabs) =====\n")
  print(VarCorr(m))
  print(anova(m, type = 2))                              # Satterthwaite F tests
  cat("Shapiro-Wilk on residuals p =", signif(shapiro.test(resid(m))$p.value, 2), "\n")
  sim <- simulateResiduals(m, n = 500)
  cat("DHARMa uniformity p =", round(testUniformity(sim, plot = FALSE)$p.value, 3),
      "| outlier p =", round(testOutliers(sim, plot = FALSE)$p.value, 3), "\n")
}
# -> C. sapidus fails normality / outlier tests because of three small crabs
#    (43, 61, 72 mm; two recorded immature, the 43-mm crab far below L50).

# --- 4. Final models: crabs >= 75 mm CW ---------------------------------------
# Removes the three immature C. sapidus; no C. bocourti are affected.
adults <- filter(crabs, cw_spines_mm >= 75)

for (sp in levels(crabs$species)) {
  dat <- filter(adults, species == sp)
  m   <- lmer(cw_spines_mm ~ sex + station + survey + (1 | ss), data = dat)
  cat("\n=====", sp, "(>= 75 mm, n =", nobs(m), ") =====\n")
  print(VarCorr(m))
  print(anova(m, type = 2))
  cat("Shapiro-Wilk p =", signif(shapiro.test(resid(m))$p.value, 2), "\n")
  print(ranova(m))                                         # LRT for the random intercept
  cat("\nEstimated marginal means (mm):\n")
  print(emmeans(m, ~ sex,     lmer.df = "satterthwaite"))
  print(emmeans(m, ~ station, lmer.df = "satterthwaite"))
  print(emmeans(m, ~ survey,  lmer.df = "satterthwaite"))
  cat("\nTukey pairwise, survey:\n")
  print(pairs(emmeans(m, ~ survey, lmer.df = "satterthwaite"), adjust = "tukey"))
}

# --- 5. Size-frequency distributions (5-mm classes, pooled) -----------------
ggplot(adults, aes(cw_spines_mm, fill = sex)) +
  geom_histogram(binwidth = 5, boundary = 0, colour = "white", position = "stack") +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  labs(x = "Carapace width (mm)", y = "Number of crabs") + theme_bw()
ggsave("fig_size_frequency.png", width = 6, height = 6, dpi = 300)
