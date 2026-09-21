# =============================================================================
# Overall sex ratio of Callinectes spp. in Torrecilla Lagoon,
# accounting for the study design
#
# Design: 4 stations x 5 quarterly surveys (20 station-survey combinations),
# with 2-3 consecutive trap-days per combination (54 trap-day events).
# Crabs from the same station-survey combination are not independent, so a
# simple binomial test on pooled counts overstates precision.
#
# Approach: one intercept-only binomial GLMM per species,
#     male ~ 1 + (1 | station:survey)
# The intercept is the overall log-odds that a crab is male; the random
# intercept for station-survey combination absorbs spatial and temporal
# structure (the same random term used in the catch models, Eq. 3).
# H0 (1:1 sex ratio)  <=>  intercept = 0 (Wald z test).
# =============================================================================

library(dplyr)
library(lme4)
library(DHARMa)
library(ggplot2)

# --- 1. Data -----------------------------------------------------------------
crabs <- read.csv("../../data/clean/crabdata_clean_1.csv") %>%
  filter(sex %in% c("F", "M"),
         species %in% c("Callinectes sapidus", "Callinectes bocourti",
                        "Callinectes similis", "Callinectes ornatus")) %>%
  mutate(male    = as.integer(sex == "M"),
         species = factor(species, levels = c("Callinectes sapidus", "Callinectes bocourti",
                                              "Callinectes similis", "Callinectes ornatus")),
         station = factor(station, levels = c("BC", "AP", "DY", "IM")),
         survey  = factor(survey),
         ss      = interaction(station, survey, drop = TRUE))   # station-survey combination

# --- 2. Exploratory ----------------------------------------------------------
# 2a. Pooled counts and proportion male (naive estimate), with exact binomial CI
raw <- crabs %>%
  group_by(species) %>%
  summarise(n = n(), males = sum(male), females = n - males,
            prop_male   = males / n,
            naive_low   = binom.test(males, n)$conf.int[1],
            naive_high  = binom.test(males, n)$conf.int[2],
            n_ss        = n_distinct(ss),          # station-survey combos represented
            .groups = "drop")
print(raw, digits = 3)

# 2b. Sex counts by station and by survey (cf. Tables 3-4)
print(with(crabs, table(species, station, sex)))
print(with(crabs, table(species, survey, sex)))

# 2c. Proportion male per station-survey combination: is there structure?
ss_sum <- crabs %>%
  group_by(species, station, survey, ss) %>%
  summarise(n = n(), males = sum(male), prop_male = males / n, .groups = "drop")

ggplot(ss_sum, aes(station, prop_male, size = n, colour = survey)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  facet_wrap(~ species) +
  labs(y = "Proportion male (station-survey combination)", size = "Crabs") +
  theme_bw()
ggsave("fig_sexratio_exploratory.png", width = 8, height = 6, dpi = 300)

# 2d. Extra-binomial variation among station-survey combinations
#     (dispersion > 1 means more heterogeneity than binomial sampling alone)
for (sp in levels(crabs$species)) {
  g <- glm(cbind(males, n - males) ~ 1, family = binomial, data = filter(ss_sum, species == sp))
  cat(sp, ": dispersion =", round(sum(residuals(g, "pearson")^2) / g$df.residual, 2), "\n")
}

# --- 3. Fit one GLMM per species --------------------------------------------
fits <- lapply(setNames(nm = levels(crabs$species)), function(sp)
  glmer(male ~ 1 + (1 | ss), family = binomial,
        data = filter(crabs, species == sp),
        control = glmerControl(optimizer = "bobyqa")))

# --- 4. Diagnostics ----------------------------------------------------------
for (sp in names(fits)) {
  cat("\n", sp, "\n"); print(VarCorr(fits[[sp]]))
  sim <- simulateResiduals(fits[[sp]], n = 1000)
  cat(" DHARMa uniformity p =", round(testUniformity(sim, plot = FALSE)$p.value, 2),
      "| dispersion p =",       round(testDispersion(sim, plot = FALSE)$p.value, 2), "\n")
}
# Likelihood-ratio test of the random effect (boundary-corrected p = 0.5 * chi-square p)
for (sp in names(fits)) {
  g0 <- glm(male ~ 1, family = binomial, data = filter(crabs, species == sp))
  lrt <- 2 * (as.numeric(logLik(fits[[sp]])) - as.numeric(logLik(g0)))
  cat(sp, ": LRT chi2 =", round(lrt, 2), " p =", signif(0.5 * pchisq(lrt, 1, lower.tail = FALSE), 3), "\n")
}

# --- 5. Results --------------------------------------------------------------
res <- bind_rows(lapply(names(fits), function(sp) {
  m  <- fits[[sp]]
  b  <- fixef(m)[[1]]; se <- sqrt(vcov(m)[1, 1]); z <- b / se
  data.frame(species = sp,
             glmm_prop_male = plogis(b),
             glmm_low       = plogis(b - 1.96 * se),
             glmm_high      = plogis(b + 1.96 * se),
             odds_male      = exp(b),               # males per female
             z = z, p = 2 * pnorm(-abs(z)),
             sd_station_survey = as.data.frame(VarCorr(m))$sdcor[1])
}))
res <- left_join(raw, res, by = "species")
print(as.data.frame(res), digits = 3)
write.csv(res, "sexratio_results.csv", row.names = FALSE)
