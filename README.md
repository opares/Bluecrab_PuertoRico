# A Caribbean coastal lagoon serves as spawning habitat for multiple Callinectes spp.

## Overviewn
This repository contains data and code associated with the study:

Pares, O., Stevens, B., & Schott, E. (2026). *A Caribbean coastal lagoon serves as spawning habitat for multiple Callinectes spp.: implications for fishery management.*

This study investigates the relative abundance, reproductive state, size structure, and life history of Callinectes species in Torrecilla Lagoon, Puerto Rico. The goal is to provide baseline data for a data-poor artisanal fishery.

## Key Findings
- 349 crabs collected across 7 species
- C. sapidus and C. bocourti comprised 65% of total catch
- Evidence supports Torrecilla Lagoon as a spawning habitat
- Size at maturity:
  - C. sapidus: 113 mm CW
  - C. bocourti: 109 mm CW
- Estimated spawning recovery time: ~21 days

## Study Design
- Quarterly sampling (Oct 2020 – Oct 2021)
- 4 stations along estuarine gradient
- 5 traps per station (mesh: 3.8 cm)
- Environmental variables recorded (YSI ProDSS)

## Repository Structure
Everything is inside Bluecrab_PuertoRico_submission_Final_Revised/:
code/
  bluecrab_puertoRico_revised.Rmd     all analyses in the paper, in manuscript order
  make_analysis_ready.R               builds data/analysis_ready/ from data/clean/
  fig_S3_bars.R                       Fig. S3 (bar version)
  standalone_scripts/
    catch_model_check.R               catch GAMLSS refitted on the clean data (Table 2)
    sex_ratio_analysis.R              binomial GLMM + Fisher tests (Tables 3–4)
    size_analysis.R                   size-frequency and carapace-width LMM (Fig. 5)
    reproductive_state_analysis.R     reproductive state by survey (Fig. 6, Table 5, Table S3)
    fig_size_frequency.R              Fig. 5
data/
  clean/                              one row per crab / one row per sampling event
  analysis_ready/                     tables derived from data/clean by make_analysis_ready.R
figures/                              Figures 1–8 and S1–S3 as submitted


Data
data/clean/
File	Rows	Description
crabdata_clean_2.csv	
**Use this file.**
349 crabs	One row per crab: crab_id, event_id, survey, date_catch, station, species, sex, maturity, carapace width spine-to-spine (cw_spines_mm) and base-to-base (cw_base_mm), carapace_length_mm, abdomen_width_mm, weight_g, ovary and egg data (ovary_weight_g, ovary_stage, egg_stage, prior_spawn), derived repro_state and spawning_capable, and a data_note column recording, for example, the six females that could not be dissected in the field. 

sampling_events_final_clean_2.csv	
** Use this file.**
54 events	One row per station × trap-day: survey, date_catch, station (BC = Station 1, AP = 2, DY = 3, IM = 4), number of traps, YSI sonde readings (temp_c, do_pct_sat, salinity_ppt, ph, depth_m, with the cast date), catch by species (n_sapidus … n_unidentified), n_female, n_male, total_catch and qc_flags (missing sonde readings are flagged rather than filled).


crabdata_clean_1.csv, sampling_events_final_clean.csv		Earlier versions kept for traceability. _1 lacks the ovary stages recovered from the field datasheets for three crabs and the repro_state columns; the events file without _2 carries a sonde reading for BC on 15 Oct 2021 that was later found to have no cast.

Station codes: BC Boca de Cangrejo (Station 1, inlet), AP airport (Station 2), DY Suárez channel (Station 3), IM Isla Mosquito (Station 4).

data/analysis_ready/

Built from data/clean/ by code/make_analysis_ready.R; each file feeds one analysis or table in the paper.

File	Used for
env_station_survey_means.csv	Station × survey means of the sonde readings (n = 20) — Table 1, Section 3.1, Table S1
sexratio_crabs.csv	One row per sexed crab of the four well-sampled species — Section 3.4 GLMM
table3_sex_by_survey_counts.csv, table4_sex_by_station_counts.csv	Tables 3 and 4
size_all_measured_crabs.csv	Carapace width of every measured crab — Fig. 5, Table S2
size_lmm_crabs_ge75mm.csv	Crabs ≥ 75 mm CW used in the carapace-width mixed model — Section 3.3
females_reproductive_state.csv	245 mature females of the four species with reproductive state — Fig. 6, Fig. S3, Table 5, Table S3
table5_reproductive_state_by_survey.csv	Proportion ovigerous / spawning-capable per species and survey with exact 95% CI — Table 5, Fig. S3
table2_catch_model_refit.csv	Wald F tests from the catch GAMLSS fitted to the 54 clean events (clean events (54 rows); Table 2) alongside the fit to a superseded 53-row dataset, kept to document the change
Code

code/bluecrab_puertoRico_revised.Rmd 
reproduces every analysis in the order of the paper: 
1. environmental parameters (two-way ANOVA, Friedman check)
2. catch (GAMLSS negative binomial with AR1 correlation within station–survey)
3. sex ratio (binomial GLMM and Fisher's exact tests)
4. size (size-frequency distributions and linear mixed model)
5. reproductive state (exact binomial CIs, Fisher's exact tests, binomial GLMM, sensitivity table)
6.  weight–length relationships.

Figures are written to figures/.

Paths in the Rmd are relative to the package folder, so set the working directory to it first:

r
setwd("path/to/Bluecrab_PuertoRico_submission_Final_Revised")
rmarkdown::render("code/bluecrab_puertoRico_revised.Rmd")

or, if knitting from RStudio with the Rmd open, add knitr::opts_knit$set(root.dir = "..") to the first chunk.

R ≥ 4.3 with dplyr, tidyr, ggplot2, lme4, lmerTest, emmeans, DHARMa, gamlss and nlme.

Figures

figures/ contains the eight main-text figures (survey map, catch by survey and station, sex ratio by station, size-frequency distributions, reproductive state by survey, length–weight relationships, carapace width–weight relationship of female C. sapidus) and the three supplementary figures (catch per trap, occurrence, reproductive proportions with confidence intervals).

Notes on the data
Six females were not dissected because field time ran out; they are retained in crabdata_clean_2.csv with a data_note and excluded from the reproductive-state analyses.
Seven crabs were recorded as immature; maturity was not recorded for 63 crabs.
Sonde readings are missing for AP and IM on 20 Jan 2021, pH for DY on 20 Oct 2020, and the whole cast for BC on 15 Oct 2021 (see qc_flags).
Raw field datasheets are available from the corresponding author on request.


