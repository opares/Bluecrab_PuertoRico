# Re-run the catch model exactly as in bluecrab_puertoRico.Rmd (GitHub), on
# (a) the GitHub d1.csv (53 rows) and (b) the clean sampling-events file (54 rows).
suppressPackageStartupMessages({library(gamlss); library(nlme); library(dplyr)})

fit_one <- function(dat, resp) {
  f <- as.formula(paste0(resp, " ~ re(fixed = ~ Year_Month + Station, random = ~1|lagoon, correlation = corAR1(form = ~ Date | lagoon/sq), method = 'REML')"))
  assign("dat", dat, envir = .GlobalEnv)
  m <- gamlss(f, family = NBI(), data = dat, trace = FALSE)
  a <- anova(getSmo(m))
  data.frame(resp = resp, term = rownames(a), numDF = a$numDF, denDF = a$denDF, F = round(a[["F-value"]], 2), p = signif(a[["p-value"]], 2))
}

prep <- function(d) {
  d$sq <- factor(d$sq); d$Year_Month <- factor(d$Year_Month); d$Station <- factor(d$Station); d$lagoon <- factor(d$lagoon)
  d$Date <- as.Date(d$Date); d
}

# (a) GitHub data
d1 <- read.csv("d1.csv", check.names = FALSE)

d1
d1$Date <- as.Date(d1$Date, format = "%m/%d/%y")
d1 <- prep(d1)
d1$total_catch <- with(d1, cpue.sap + cpue.boc + cpue.sim + cpue.orn + cpue.lar + cpue.exa)
resA <- bind_rows(lapply(c("total_catch", "cpue.sap", "cpue.boc", "cpue.sim", "cpue.orn"), function(r) fit_one(d1, r)))
resA$data <- "github_d1 (53 rows)"

# (b) clean events file
s <- read.csv("sampling_events_final_clean_2.csv")
#or 
#s<- sampling_events_final_clean_2

d2 <- data.frame(Date = as.Date(s$date_catch, format = "%m/%d/%y"), Station = s$station, Year_Month = s$survey,
                 lagoon = "lagoon", sq = paste(s$survey, s$station, sep = "."),
                 cpue.sap = s$n_sapidus, cpue.boc = s$n_bocourti, cpue.sim = s$n_similis, cpue.orn = s$n_ornatus,
                 cpue.lar = s$n_larvatus, cpue.exa = s$n_exasperatus)
d2 <- prep(d2)
d2$total_catch <- with(d2, cpue.sap + cpue.boc + cpue.sim + cpue.orn + cpue.lar + cpue.exa)
resB <- bind_rows(lapply(c("total_catch", "cpue.sap", "cpue.boc", "cpue.sim", "cpue.orn"), function(r) fit_one(d2, r)))
resB$data <- "clean events (54 rows)"

out <- bind_rows(resA, resB)
print(out, row.names = FALSE)
#write.csv(out, "catch_model_check.csv", row.names = FALSE)

# random-effect variance and AR1 parameter for the pooled model (clean data)
m <- gamlss(total_catch ~ re(fixed = ~ Year_Month + Station, random = ~1|lagoon, correlation = corAR1(form = ~ Date | lagoon/sq), method = "REML"), family = NBI(), data = d2, trace = FALSE)
print(VarCorr(getSmo(m))); print(getSmo(m)$modelStruct$corStruct)
