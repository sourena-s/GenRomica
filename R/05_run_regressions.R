# =============================================================================
# Step 5 - Association of genetic scores with individual outcomes in GenR
# =============================================================================
# Each genetic score (ICA component scores, polygenic scores) is tested against
#   (a) repeated measures (CBCL, MRI, DTI) with linear mixed models, and
#   (b) outcomes measured once (SES, diet, online) with linear regression.
# Every score is tested alone; each set in cfg$combined_sets is also tested
# jointly (all scores of the set in one model).
#
# Hypotheses
#   H1  the score(s) are associated with the outcome
#   H2  (single scores, repeated measures) the association changes with age
#
# The statistical models are described in R/models.R. p-values
# are FDR-corrected (Benjamini-Hochberg) within each score set and hypothesis.
#
# Output (cfg$results_dir)
#   results.csv             one row per score set x score x outcome x hypothesis
#   h2_details.csv          age-interaction estimates behind H2
#   covariate_results.csv   covariate estimates of the H2 mixed models
#   failed_models.csv       models that stopped with an error
#   model_warnings.csv      warnings raised while fitting (e.g. convergence)

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)     # Satterthwaite t- and F-tests for lmer models
  library(performance)  # Nakagawa R2
})
source("R/utils.R")
cfg <- load_config()

gpcs <- cfg$genetic_pcs
all_scores <- unlist(cfg$score_sets, use.names = FALSE)


# =============================================================================
# 1. Load data
# =============================================================================

genr <- load_genr_with_scores(cfg)
message(sprintf("Loaded %d individuals, %d scores", nrow(genr), length(all_scores)))


# =============================================================================
# 2. Repeated measures in long format
# =============================================================================
# One row per individual x wave x outcome, keeping observations with both an outcome
# value and an age. Columns: idc, covariates, age, outcome, domain,
# outcome_name, wave.

long_data <- make_long_data(cfg, genr)
long_outcomes <- unique(cfg$longitudinal[, .(domain, outcome)])


# =============================================================================
# 3-6. Models (R/models.R): helpers, mixed models for repeated measures (one
# score at a time, or a score set jointly) and linear regression for outcomes
# measured once. The statistical models are described there.
# =============================================================================

source("R/models.R")


# =============================================================================
# 7. Run all models
# =============================================================================

all_results    <- list()
all_h2         <- list()
all_covariates <- list()

add_labels <- function(x, score_group, domain, outcome) {
  x[, `:=`(score_group = score_group, domain = domain, outcome = outcome)]
}

# ---- Each score on its own ---------------------------------------------------
for (score_group in names(cfg$score_sets)) {
  for (score in cfg$score_sets[[score_group]]) {
    message("Running ", score_group, ": ", score)

    for (k in seq_len(nrow(long_outcomes))) {
      dom <- long_outcomes$domain[k]
      out <- long_outcomes$outcome[k]

      fit <- safe_fit(fit_lmm_one_score(long_data[domain == dom & outcome_name == out], score),
                      paste(score_group, score, dom, out, sep = " | "))
      if (is.null(fit)) next

      all_results[[length(all_results) + 1]]       <- add_labels(fit$main, score_group, dom, out)
      all_h2[[length(all_h2) + 1]]                 <- add_labels(fit$h2_details, score_group, dom, out)
      all_covariates[[length(all_covariates) + 1]] <- add_labels(fit$covariates, score_group, dom, out)
    }

    for (k in seq_len(nrow(cfg$single_outcomes))) {
      o <- cfg$single_outcomes[k]

      fit <- safe_fit(fit_lm(o$outcome, score, score, combined = FALSE, age_var = o$age_var),
                      paste(score_group, score, o$domain, o$outcome, sep = " | "))
      if (is.null(fit)) next

      all_results[[length(all_results) + 1]] <- add_labels(fit, score_group, o$domain, o$outcome)
    }
  }
}

# ---- All scores of a set jointly ---------------------------------------------
for (score_group in names(cfg$combined_sets)) {
  scores <- cfg$combined_sets[[score_group]]
  message("Running ", score_group)

  for (k in seq_len(nrow(long_outcomes))) {
    dom <- long_outcomes$domain[k]
    out <- long_outcomes$outcome[k]

    fit <- safe_fit(fit_lmm_score_set(long_data[domain == dom & outcome_name == out], scores, score_group),
                    paste(score_group, dom, out, sep = " | "))
    if (is.null(fit)) next

    all_results[[length(all_results) + 1]] <- add_labels(fit, score_group, dom, out)
  }

  for (k in seq_len(nrow(cfg$single_outcomes))) {
    o <- cfg$single_outcomes[k]

    fit <- safe_fit(fit_lm(o$outcome, scores, score_group, combined = TRUE, age_var = o$age_var),
                    paste(score_group, o$domain, o$outcome, sep = " | "))
    if (is.null(fit)) next

    all_results[[length(all_results) + 1]] <- add_labels(fit, score_group, o$domain, o$outcome)
  }
}


# =============================================================================
# 8. Multiple-testing correction and output
# =============================================================================

if (!length(all_results)) stop("No models ran successfully.")

results           <- rbindlist(all_results, fill = TRUE)
h2_details        <- rbindlist(all_h2, fill = TRUE)
covariate_results <- rbindlist(all_covariates, fill = TRUE)
failed_models     <- rbind(data.table(model = character(), error = character()),
                           rbindlist(failures, fill = TRUE))
warnings_table    <- rbind(data.table(model = character(), warning = character()),
                           rbindlist(model_warnings, fill = TRUE))

# Benjamini-Hochberg FDR within each score set and hypothesis (over all
# predictors and outcomes of the set)
results[, p_fdr := p.adjust(p_value, method = "BH"), by = .(score_group, hypothesis)]

if (nrow(h2_details)) {
  h2_fdr <- results[hypothesis == "H2",
                    .(predictor, score_group, domain, outcome, joint_p_fdr = p_fdr)]
  h2_details <- merge(h2_details, h2_fdr,
                      by = c("predictor", "score_group", "domain", "outcome"), all.x = TRUE)
}

dir.create(cfg$results_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(results,           file.path(cfg$results_dir, "results.csv"))
fwrite(h2_details,        file.path(cfg$results_dir, "h2_details.csv"))
fwrite(covariate_results, file.path(cfg$results_dir, "covariate_results.csv"))
fwrite(failed_models,     file.path(cfg$results_dir, "failed_models.csv"))
fwrite(warnings_table,    file.path(cfg$results_dir, "model_warnings.csv"))

message(sprintf("Done: %d results, %d failed models, %d warnings. Output in %s",
                nrow(results), nrow(failed_models), nrow(warnings_table), cfg$results_dir))
