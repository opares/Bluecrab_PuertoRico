# Figures 8 and 9: weight-length and width-weight relationships (Section 3.6)
#
# Uses the same input file and filters as length_weight_analysis.R:
#   non-ovigerous crabs; immature crabs dropped; crabs with no maturity stage
#   recorded kept if >= 75 mm CW; CR306's spine width is blank in the data.
# Fitted curves are W = a * X^b from log10-log10 regressions.
#
# Run from the Bluecrab_PuertoRico_submission_Final_Revised folder. Writes:
#   figures/Figure_8_L-W_relationship_M_F.png   (panels A and B side by side)
#   figures/Figure_8A_sapidus.png, figures/Figure_8B_bocourti.png (single panels)
#   figures/Figure_9_CW_Weight_Female_sap.png

library(dplyr)
library(ggplot2)

wl <- read.csv("data/analysis_ready/weight_length_crabs.csv", na.strings = c("", "NA"))
wl$cw_spines_mm[wl$crab_id %in% "CR306"] <- NA   # recording error (see README)

wl <- wl %>%
  filter(ovigerous_original == 0,
         !maturity_original %in% "IM",                       # drop immature
         !is.na(maturity_original) | cw_spines_mm >= 75) %>%  # unlabelled: keep if >= 75 mm
  filter(!is.na(weight_g), !is.na(carapace_length_mm))

sap <- filter(wl, species == "Callinectes sapidus")
boc <- filter(wl, species == "Callinectes bocourti")

# ---- Figure 8: weight vs carapace length, by sex ---------------------------
col_sex <- c(F = "blue", M = "red")
pch_sex <- c(F = 1, M = 16)

panel_lw <- function(d, label = NULL) {
  plot(weight_g ~ carapace_length_mm, data = d,
       pch = pch_sex[d$sex], col = col_sex[d$sex], cex = 0.8,
       xlab = "Carapace Length (mm)", ylab = "Weight (g)", las = 0)
  for (s in c("F", "M")) {
    ds  <- filter(d, sex == s)
    fit <- lm(log10(weight_g) ~ log10(carapace_length_mm), data = ds)
    x   <- seq(min(ds$carapace_length_mm), max(ds$carapace_length_mm), length.out = 200)
    lines(x, 10^predict(fit, data.frame(carapace_length_mm = x)), col = col_sex[s], lwd = 2)
  }
  legend("topleft", legend = c("Female", "Male"), col = col_sex, pch = pch_sex,
         lwd = 2, bty = "n", cex = 0.9)
  if (!is.null(label)) mtext(label, side = 3, adj = 0, line = 0.5, font = 2, cex = 1.2)
}

png("figures/Figure_8_L-W_relationship_M_F.png", width = 7, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
panel_lw(sap, "A"); panel_lw(boc, "B")
dev.off()

# single panels, without letters (the manuscript adds A and B as labels)
for (p in list(list(sap, "figures/Figure_8A_sapidus.png"), list(boc, "figures/Figure_8B_bocourti.png"))) {
  png(p[[2]], width = 2.8, height = 4.06, units = "in", res = 300, pointsize = 9)
  par(mar = c(4, 4, 1, 0.5)); panel_lw(p[[1]]); dev.off()
}

# ---- Figure 9: weight vs carapace width, female C. sapidus -----------------
sap_f  <- filter(sap, sex == "F", !is.na(cw_spines_mm))
fit_cw <- lm(log10(weight_g) ~ log10(cw_spines_mm), data = sap_f)
a <- 10^coef(fit_cw)[1]; b <- coef(fit_cw)[2]
cat(sprintf("Fig. 9 curve: W = %.5f CW^%.2f, R2 = %.2f, n = %d\n", a, b, summary(fit_cw)$r.squared, nobs(fit_cw)))

curve_df <- data.frame(cw = seq(min(sap_f$cw_spines_mm), max(sap_f$cw_spines_mm), length.out = 200))
curve_df$w <- a * curve_df$cw^b

p9 <- ggplot(sap_f, aes(cw_spines_mm, weight_g)) +
  geom_point(colour = "blue", alpha = 0.8, size = 2) +
  geom_line(data = curve_df, aes(cw, w), colour = "blue", linewidth = 1) +
  labs(x = "Carapace width (lateral spines) (mm)", y = "Weight (g)") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), axis.title = element_text(face = "bold"))
ggsave("figures/Figure_9_CW_Weight_Female_sap.png", p9, width = 4.1, height = 4.56, dpi = 300)
