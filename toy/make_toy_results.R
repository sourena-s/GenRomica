# =============================================================================
# Simulated association results, for designing the result figures
# =============================================================================
# Writes toy/data/viz/results.csv in the format of step 5 (results.csv), for
# 20 GenomICA components (C1-C20) and every outcome in config/project.R:
#   H1 for all 17 outcomes, H2 for the 12 repeated measures.
# No models are fitted. Most associations are null; a few planted effects make
# the figures show significant and FDR-significant cells.
#
# Usage (from the pipeline folder): Rscript toy/make_toy_results.R

suppressPackageStartupMessages(library(data.table))
source("config/project.R")
set.seed(20)

components <- paste0("C", 1:20)
long  <- unique(cfg$longitudinal[, .(domain, outcome)])
single <- cfg$single_outcomes[, .(domain, outcome)]
outcomes <- rbind(long[, type := "repeated"], single[, type := "single"])

# Planted effects (standardised beta) for H1, and planted age interactions for H2
h1_effects <- data.table(
  predictor = c("C3", "C3", "C3", "C7", "C7", "C12", "C15", "C18", "C9"),
  outcome   = c("att", "agg", "rul", "EDUCM", "INCOME5", "TotalGrayVol", "onlineweek_child_F17", "anx", "mean_FA_genr"),
  beta      = c(0.09, 0.07, 0.06, -0.11, -0.08, 0.06, 0.07, 0.05, -0.05)
)
h2_effects <- data.table(predictor = c("C3", "C12"), outcome = c("att", "cortical_thickness"))

h1 <- CJ(predictor = components, outcome = outcomes$outcome, sorted = FALSE)
h1 <- merge(h1, outcomes, by = "outcome", sort = FALSE)
h1[, n_id := ifelse(type == "repeated", round(runif(.N, 2500, 5500)), round(runif(.N, 3000, 7000)))]
h1[, std_error := 1 / sqrt(n_id)]
h1 <- merge(h1, h1_effects, by = c("predictor", "outcome"), all.x = TRUE, sort = FALSE)
h1[is.na(beta), beta := 0]
h1[, estimate := beta + rnorm(.N, 0, std_error)]
h1[, p_value := 2 * pnorm(-abs(estimate / std_error))]
h1[, incremental_r2 := estimate^2]
h1[, hypothesis := "H1"]

h2 <- CJ(predictor = components, outcome = long$outcome, sorted = FALSE)
h2 <- merge(h2, outcomes, by = "outcome", sort = FALSE)
h2[, n_id := round(runif(.N, 2500, 5500))]
h2[, f := rchisq(.N, df = 2) / 2]
h2[h2_effects, on = .(predictor, outcome), f := f + 9]
h2[, p_value := pf(f, 2, n_id * 2, lower.tail = FALSE)]
h2[, `:=`(estimate = NA_real_, std_error = NA_real_, incremental_r2 = f * 2 / (n_id * 2), hypothesis = "H2")]

results <- rbind(h1, h2, fill = TRUE)
results[, `:=`(score_group = "PCS", singular = FALSE, n_obs = n_id * fifelse(type == "repeated", 3, 1))]
results[, p_fdr := p.adjust(p_value, method = "BH"), by = .(score_group, hypothesis)]
results <- results[, .(hypothesis, estimate, std_error, p_value, incremental_r2, predictor,
                       n_id, n_obs, singular, score_group, domain, outcome, p_fdr)]

dir.create("toy/data/viz", recursive = TRUE, showWarnings = FALSE)
fwrite(results, "toy/data/viz/results.csv")
message(sprintf("Wrote toy/data/viz/results.csv: %d rows (%d H1, %d H2); %d nominal and %d FDR-significant",
                nrow(results), sum(results$hypothesis == "H1"), sum(results$hypothesis == "H2"),
                sum(results$p_value < 0.05), sum(results$p_fdr < 0.05)))
