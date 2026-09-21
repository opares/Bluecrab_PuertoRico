# Size-frequency distributions (5-mm CW classes) of the four well-sampled
# Callinectes species, by sex, with individuals recorded as immature shown in white.
library(dplyr)
library(ggplot2)

crabs <- read.csv("../../data/clean/crabdata_clean_1.csv") %>%
  filter(species %in% c("Callinectes sapidus", "Callinectes bocourti",
                        "Callinectes similis", "Callinectes ornatus"),
         sex %in% c("F", "M"), !is.na(cw_spines_mm)) %>%
  mutate(species = factor(species, levels = c("Callinectes sapidus", "Callinectes bocourti",
                                              "Callinectes similis", "Callinectes ornatus"),
                          labels = c("C. sapidus", "C. bocourti", "C. similis", "C. ornatus")),
         group = case_when(maturity_coded == "immature" ~ "Immature",
                           sex == "F" ~ "Female", TRUE ~ "Male"),
         group = factor(group, levels = c("Immature", "Male", "Female")))

# panel label carries n and number immature
crabs <- crabs %>% group_by(species) %>%
  mutate(panel = sprintf("%s  (n = %d; %d immature)", species, n(), sum(group == "Immature"))) %>%
  ungroup() %>% mutate(panel = factor(panel, levels = unique(panel[order(species)])))

p <- ggplot(crabs, aes(cw_spines_mm, fill = group)) +
  geom_histogram(binwidth = 5, boundary = 0, colour = "grey20", linewidth = 0.25, position = "stack") +
  scale_fill_manual(values = c(Female = "#E69F00", Male = "#0072B2", Immature = "white"),
                    breaks = c("Female", "Male", "Immature"), name = NULL) +
  scale_x_continuous(breaks = seq(40, 160, 20), limits = c(40, 165)) +
  scale_y_continuous(breaks = function(l) seq(0, floor(l[2]), by = ifelse(l[2] > 15, 5, 2)), expand = expansion(mult = c(0, 0.08))) +
  facet_wrap(~ panel, ncol = 2, scales = "free_y") +
  labs(x = "Carapace width (mm)", y = "Number of crabs") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "italic", hjust = 0),
        legend.position = "top", legend.justification = "left", panel.grid.major.y = element_line(colour = "grey92"))

ggsave("Fig6_size_frequency.png", p, width = 7, height = 5.5, dpi = 300)
ggsave("Fig6_size_frequency.pdf", p, width = 7, height = 5.5)
