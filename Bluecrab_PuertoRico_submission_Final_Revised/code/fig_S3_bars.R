# Fig. S3 (bar version): proportion of mature females ovigerous / spawning-capable
# per survey, with exact binomial 95% CI. Reads the analysis-ready Table 5 file.
library(dplyr); library(tidyr); library(ggplot2)

tab <- read.csv("data/analysis_ready/table5_reproductive_state_by_survey.csv")

survey_lab <- c("2020-10" = "Oct\n2020", "2021-01" = "Jan\n2021", "2021-04" = "Apr\n2021",
                "2021-07" = "Jul\n2021", "2021-10" = "Oct\n2021")

long <- bind_rows(
  tab %>% transmute(species, survey, n, measure = "Spawning-capable", p = p_cap,  lo = cap_lo,  hi = cap_hi),
  tab %>% transmute(species, survey, n, measure = "Ovigerous",        p = p_ovig, lo = ovig_lo, hi = ovig_hi)
) %>%
  mutate(species = factor(species, levels = c("Callinectes sapidus", "Callinectes bocourti",
                                              "Callinectes similis", "Callinectes ornatus"),
                          labels = c("C. sapidus", "C. bocourti", "C. similis", "C. ornatus")),
         survey  = factor(survey, levels = names(survey_lab), labels = survey_lab),
         measure = factor(measure, levels = c("Spawning-capable", "Ovigerous")))

pal <- c("Spawning-capable" = "#A52A2A", "Ovigerous" = "#FFA501")
nlab <- long %>% distinct(species, survey, n)

p <- ggplot(long, aes(survey, p, fill = measure)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "grey20", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.75),
                width = 0.25, linewidth = 0.4, colour = "grey20") +
  geom_text(data = nlab, aes(survey, y = -0.08, label = paste0("n = ", n)),
            inherit.aes = FALSE, colour = "grey30", size = 2.6) +
  facet_wrap(~ species, nrow = 1, drop = FALSE) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(limits = c(-0.12, 1.02), breaks = seq(0, 1, 0.25),
                     labels = c("0", "0.25", "0.50", "0.75", "1")) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Survey", y = "Proportion of mature females") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "italic"),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggsave("figures/Figure_S3_reproductive_proportions_CI.png", p, width = 9.5, height = 4.2, dpi = 300)
