# =============================================================================
# Reproductive state of mature female Callinectes spp. by survey
#
# Question: was there any quarter in which females were not reproducing?
# Response: "spawning-capable" = ovigerous OR full ovaries OR egg remnants
#           (prior spawn), vs. light/medium ovaries with no eggs.
#
# Two complementary tests per species:
#   (a) Fisher's exact test of spawning-capable x survey (exact, no minimum
#       cell size, but treats crabs as independent);
#   (b) binomial GLMM  capable ~ survey + (1 | station:survey), which allows
#       for the non-independence of crabs trapped at the same station in the
#       same survey (same random term as the catch and sex-ratio models),
#       with survey tested by likelihood-ratio test.
# =============================================================================

library(dplyr)
library(lme4)
library(DHARMa)
library(ggplot2)

# --- 1. Data and classification --------------------------------------------
# crabdata_clean_2.csv: CR149 and CR344 (ovary examined, not stageable) are
# coded ovary_stage = "light"; see data_note column.
fem <- read.csv("../../data/clean/crabdata_clean_2.csv") %>%
  filter(sex == "F",
         species %in% c("Callinectes sapidus", "Callinectes bocourti",
                        "Callinectes similis", "Callinectes ornatus"),
         is.na(maturity_coded) | maturity_coded != "immature") %>%   # drop recorded immature
  mutate(species = factor(species, levels = c("Callinectes sapidus", "Callinectes bocourti",
                                              "Callinectes similis", "Callinectes ornatus")),
         survey  = factor(survey),
         ss      = interaction(station, survey, drop = TRUE),
         ovig    = ovigerous   == "Y",
         prior   = prior_spawn == "Y",
         full    = ovary_stage == "full",
         state   = case_when(ovig  ~ "Ovigerous",
                             prior ~ "Prior spawn",
                             full  ~ "Full ovaries",
                             ovary_stage %in% c("light", "medium") ~ "Light/medium ovaries",
                             TRUE ~ NA_character_),
         capable = state %in% c("Ovigerous", "Prior spawn", "Full ovaries"))

# 1a. How many females could not be classified (missing ovary stage and no eggs)?
cat("Unclassified females by species:\n"); print(table(fem$species, is.na(fem$state)))
fem <- filter(fem, !is.na(state))
cat("\nClassified mature females:", nrow(fem), "\n")
fem$state <- factor(fem$state, levels = c("Ovigerous", "Prior spawn", "Full ovaries", "Light/medium ovaries"))

# --- 2. Descriptive: state composition and proportion capable by survey -----
comp <- fem %>% count(species, survey, state) %>%
  group_by(species, survey) %>% mutate(n_total = sum(n), prop = n / n_total) %>% ungroup()
print(as.data.frame(comp))

cap <- fem %>% group_by(species, survey) %>%
  summarise(n = n(), k = sum(capable), p = k / n,
            lo = binom.test(k, n)$conf.int[1], hi = binom.test(k, n)$conf.int[2],   # exact CI
            n_ss = n_distinct(ss), .groups = "drop")
print(as.data.frame(cap), digits = 2)
write.csv(cap, "repro_capable_by_survey.csv", row.names = FALSE)

# Figure: stacked composition by survey, n above bars
ggplot(comp, aes(survey, prop, fill = state)) +
  geom_col(colour = "grey20", linewidth = 0.25, width = 0.75) +
  geom_text(data = distinct(comp, species, survey, n_total), aes(survey, 1.04, label = n_total),
            inherit.aes = FALSE, size = 3) +
  facet_wrap(~ species, ncol = 2) +
  scale_fill_manual(values = c("Ovigerous" = "#0072B2", "Prior spawn" = "#56B4E9",
                               "Full ovaries" = "#E69F00", "Light/medium ovaries" = "white"), name = NULL) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Survey", y = "Proportion of mature females") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "italic", hjust = 0),
        legend.position = "top", legend.justification = "left")
ggsave("fig_repro_state_by_survey.png", width = 7.5, height = 5.5, dpi = 300)

# --- 3. Test (a): Fisher's exact test, per species ---------------------------
cat("\n--- Fisher exact: capable x survey ---\n")
for (sp in levels(fem$species)) {
  tab <- with(filter(fem, species == sp), table(survey, capable))
  tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  if (ncol(tab) < 2) { cat(sp, ": all females capable in every survey - no test possible\n"); next }
  cat(sp, ": p =", signif(fisher.test(tab)$p.value, 2), "\n")
}

# --- 4. Test (b): binomial GLMM, per species ---------------------------------
cat("\n--- GLMM capable ~ survey + (1|station:survey) ---\n")
fits <- list()
for (sp in levels(fem$species)) {
  dat <- filter(fem, species == sp) %>% droplevels()
  if (length(unique(dat$capable)) < 2) { cat("\n", sp, ": no variation in response - GLMM not fitted\n"); next }
  m1  <- glmer(capable ~ survey + (1 | ss), family = binomial, data = dat,
               control = glmerControl(optimizer = "bobyqa"))
  m0  <- glmer(capable ~ 1 + (1 | ss), family = binomial, data = dat,
               control = glmerControl(optimizer = "bobyqa"))
  lrt <- anova(m0, m1)                                   # ML likelihood-ratio test of survey
  cat("\n", sp, ": n =", nrow(dat), "| station-survey SD =", round(attr(VarCorr(m1)$ss, "stddev"), 2),
      "| singular:", isSingular(m1), "\n")
  cat("  LRT survey: chi2(", lrt$Df[2], ") =", round(lrt$Chisq[2], 2), " p =", signif(lrt$`Pr(>Chisq)`[2], 2), "\n")
  # assumption / fit checks
  b <- fixef(m1); cat("  max |coef| =", round(max(abs(b)), 1),
                      ifelse(max(abs(b)) > 10, " <- separation: a survey with 0% or 100% capable", ""), "\n")
  sim <- simulateResiduals(m1, n = 1000)
  cat("  DHARMa uniformity p =", round(testUniformity(sim, plot = FALSE)$p.value, 2),
      "| dispersion p =", round(testDispersion(sim, plot = FALSE)$p.value, 2),
      "| outlier p =", round(testOutliers(sim, plot = FALSE)$p.value, 2), "\n")
  # model-based proportion capable per survey (typical station-survey combination)
  X  <- model.matrix(~ survey, data = data.frame(survey = levels(dat$survey)))
  eta <- as.vector(X %*% b); se <- sqrt(diag(X %*% vcov(m1) %*% t(X)))
  print(data.frame(survey = levels(dat$survey), p = round(plogis(eta), 2),
                   lo = round(plogis(eta - 1.96 * se), 2), hi = round(plogis(eta + 1.96 * se), 2)))
  fits[[sp]] <- m1
}

# --- 5. Sensitivity: does the definition matter? -----------------------------
# (i) Count "prior spawn" as NOT capable (egg remnants show past, not current, activity)
cat("\n--- Sensitivity: capable = ovigerous or full ovaries only ---\n")
for (sp in levels(fem$species)) {
  dat <- filter(fem, species == sp) %>% mutate(cap2 = state %in% c("Ovigerous", "Full ovaries"))
  tab <- with(dat, table(survey, cap2)); tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  if (ncol(tab) < 2) { cat(sp, ": no variation\n"); next }
  cat(sp, ": Fisher p =", signif(fisher.test(tab)$p.value, 2), "\n")
}
# (ii) Females of unrecorded maturity: are results the same on coded-mature females only?
cat("\n--- Sensitivity: recorded-mature females only ---\n")
for (sp in levels(fem$species)) {
  dat <- filter(fem, species == sp, maturity_coded == "mature")
  tab <- with(dat, table(survey, capable)); tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  if (ncol(tab) < 2) { cat(sp, ": n =", nrow(dat), " no variation\n"); next }
  cat(sp, ": n =", nrow(dat), " Fisher p =", signif(fisher.test(tab)$p.value, 2), "\n")
}

# --- 6. Ovigerous-only proportion alongside spawning-capable ----------------
both <- fem %>% group_by(species, survey) %>%
  summarise(n = n(),
            k_ovig = sum(state == "Ovigerous"), k_cap = sum(capable), .groups = "drop") %>%
  rowwise() %>%
  mutate(p_ovig = k_ovig / n, ovig_lo = binom.test(k_ovig, n)$conf.int[1], ovig_hi = binom.test(k_ovig, n)$conf.int[2],
         p_cap  = k_cap  / n, cap_lo  = binom.test(k_cap,  n)$conf.int[1], cap_hi  = binom.test(k_cap,  n)$conf.int[2]) %>%
  ungroup()
print(as.data.frame(both), digits = 2)
write.csv(both, "repro_ovig_and_capable_by_survey.csv", row.names = FALSE)

cat("\n--- Fisher exact: ovigerous x survey (same 242 females) ---\n")
for (sp in levels(fem$species)) {
  tab <- with(filter(fem, species == sp), table(survey, state == "Ovigerous"))
  tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  cat(sp, ": p =", signif(fisher.test(tab)$p.value, 2), "\n")
}

# Figure: proportion (exact 95% CI) per survey, ovigerous vs spawning-capable
long <- bind_rows(
  both %>% transmute(species, survey, n, measure = "Ovigerous",        p = p_ovig, lo = ovig_lo, hi = ovig_hi),
  both %>% transmute(species, survey, n, measure = "Spawning-capable", p = p_cap,  lo = cap_lo,  hi = cap_hi)) %>%
  mutate(measure = factor(measure, levels = c("Spawning-capable", "Ovigerous")))

ggplot(long, aes(survey, p, colour = measure, group = measure)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, position = position_dodge(0.5), linewidth = 0.5) +
  geom_point(aes(size = n), position = position_dodge(0.5)) +
  facet_wrap(~ species, ncol = 2) +
  scale_colour_manual(values = c("Spawning-capable" = "#E69F00", "Ovigerous" = "#0072B2"), name = NULL) +
  scale_size_area(max_size = 5, breaks = c(5, 10, 20, 30), name = "Females") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Survey", y = "Proportion of mature females (exact 95% CI)") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "italic", hjust = 0),
        legend.position = "top", legend.justification = "left", panel.grid.major.y = element_line(colour = "grey92"))
ggsave("fig_repro_proportions_ci.png", width = 7.5, height = 5.5, dpi = 300)
