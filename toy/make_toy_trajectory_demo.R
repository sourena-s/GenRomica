# =============================================================================
# Trajectory figure on simulated data with a known GPS x age effect
# =============================================================================
# For designing the step-7 figure only. Simulates one repeated outcome
# (attention problems) at four waves for 3,000 individuals, where the GPS
# association grows with age, and fits it with the unchanged step-5 model
# (R/models.R, fit_lmm_one_score).
#
# Usage (from the pipeline folder): Rscript toy/make_toy_trajectory_demo.R

suppressPackageStartupMessages(library(data.table))
source("R/utils.R")
source("config/project.R")
source("R/plot_style.R")
source("R/plot_trajectory.R")
set.seed(7)

n <- 3000
gpcs <- cfg$genetic_pcs
genr <- data.table(idc = 1:n,
                   sex = factor(sample(c("girl", "boy"), n, TRUE), levels = c("girl", "boy")),
                   genotype_cohort = factor(sample(1:2, n, TRUE, prob = c(.7, .3))),
                   C3 = rnorm(n))
for (v in gpcs) genr[, (v) := rnorm(n)]
source("R/models.R")

waves <- c(f05 = 6, f09 = 9.8, f13 = 13.5, f17 = 17.4)
u <- rnorm(n, 0, 0.7)                                   # individual level
long <- rbindlist(lapply(names(waves), function(w) {
  age <- waves[[w]] + rnorm(n, 0, 0.4)
  outcome <- 2 + 0.25 * (age - 6) - 0.02 * (age - 6)^2 + u +
             genr$C3 * (0.05 + 0.03 * (age - 6)) +       # GPS association grows with age
             rnorm(n, 0, 0.8)
  data.table(genr[, c("idc", "sex", "genotype_cohort", gpcs), with = FALSE],
             age = age, outcome = outcome, domain = "CBCL", outcome_name = "att", wave = w)
}))
long <- long[runif(.N) > 0.2]                           # 20% missing visits

fit <- fit_lmm_one_score(long, "C3", return_model = TRUE)
p_lines <- fit$main[, sprintf("%s: p = %s", hypothesis, formatC(p_value, format = "e", digits = 1))]
p_text <- paste(paste(p_lines[1:2], collapse = "     "),
                paste(c(p_lines[3:4], "(simulated data)"), collapse = "     "), sep = "\n")
dir.create("toy/data/viz/figures", recursive = TRUE, showWarnings = FALSE)
save_png(plot_trajectory(fit$fit, "C3", "att", p_text), "trajectory_demo_GPS3_att.png",
         width = 8.5, height = 6, dir = "toy/data/viz/figures")
