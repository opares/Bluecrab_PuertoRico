# Builds the analysis-ready tables in data/analysis_ready/ from the clean files.
# Run from the package root:  Rscript code/make_analysis_ready.R
library(dplyr)
crabs  <- read.csv("data/clean/crabdata_clean_2.csv", stringsAsFactors = FALSE)
events <- read.csv("data/clean/sampling_events_final_clean_2.csv", stringsAsFactors = FALSE)
sp4 <- c("Callinectes sapidus", "Callinectes bocourti", "Callinectes similis", "Callinectes ornatus")
out <- "data/analysis_ready"

# 1. Environmental parameters: one value per station per survey (Section 2.5.1, Table 1)
events %>% group_by(station, station_number, survey) %>%
  summarise(n_trap_days = n(), across(c(temp_c, do_pct_sat, salinity_ppt, ph), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  write.csv(file.path(out, "env_station_survey_means.csv"), row.names = FALSE)

# 2. Sex ratio: one row per sexed crab of the four species (Section 2.5.4)
crabs %>% filter(species %in% sp4, sex %in% c("F", "M")) %>%
  transmute(crab_id, species, station, survey, station_survey = paste(station, survey, sep = "_"), sex, male = as.integer(sex == "M")) %>%
  write.csv(file.path(out, "sexratio_crabs.csv"), row.names = FALSE)

# 3. Size: crabs >= 75 mm CW of C. sapidus and C. bocourti (Section 2.5.3), plus all measured crabs for Fig. 5
crabs %>% filter(species %in% sp4, sex %in% c("F", "M"), !is.na(cw_spines_mm)) %>%
  transmute(crab_id, species, station, survey, station_survey = paste(station, survey, sep = "_"), sex, maturity_coded, cw_spines_mm) %>%
  write.csv(file.path(out, "size_all_measured_crabs.csv"), row.names = FALSE)
crabs %>% filter(species %in% sp4[1:2], sex %in% c("F", "M"), !is.na(cw_spines_mm), cw_spines_mm >= 75) %>%
  transmute(crab_id, species, station, survey, station_survey = paste(station, survey, sep = "_"), sex, cw_spines_mm) %>%
  write.csv(file.path(out, "size_lmm_crabs_ge75mm.csv"), row.names = FALSE)

# 4. Reproductive state: one row per classified mature female (Section 2.5.5)
crabs %>% filter(species %in% sp4, sex == "F", repro_state != "") %>%
  transmute(crab_id, species, station, survey, station_survey = paste(station, survey, sep = "_"),
            maturity_coded, ovigerous, prior_spawn, ovary_stage, repro_state, spawning_capable) %>%
  write.csv(file.path(out, "females_reproductive_state.csv"), row.names = FALSE)

# 5. Count tables behind Tables 3 and 4
sx <- crabs %>% filter(species %in% sp4, sex %in% c("F", "M"))
sx %>% count(species, survey, sex) %>% tidyr::pivot_wider(names_from = sex, values_from = n, values_fill = 0) %>%
  mutate(pct_male = round(100 * M / (M + F), 1)) %>% write.csv(file.path(out, "table3_sex_by_survey_counts.csv"), row.names = FALSE)
sx %>% count(species, station, sex) %>% tidyr::pivot_wider(names_from = sex, values_from = n, values_fill = 0) %>%
  mutate(pct_male = round(100 * M / (M + F), 1)) %>% write.csv(file.path(out, "table4_sex_by_station_counts.csv"), row.names = FALSE)
cat("analysis-ready files written to", out, "\n")
