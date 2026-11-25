# =============================================================================
# Step 5 - Association of genetic scores with child outcomes in GenR
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
# The statistical models are described above each fitting function. p-values
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

# ---- GenR variables: covariates first, then each outcome domain --------------
tables <- genr_tables(cfg)
genr <- read_table(genr_file(cfg, "covariates"))
stop_if_duplicated_ids(genr, "covariates.csv.gz")
for (table in setdiff(names(tables), "covariates")) {
  add_columns_by_idc(genr, read_table(genr_file(cfg, table)), paste0(table, ".csv.gz"))
}
require_columns(genr, unique(unlist(tables)), "The exported GenR data",
                load_command = sprintf('d <- data.table::fread("%s")', genr_file(cfg, "<table>")))

# Categorical covariates, with the level order of the original RDS file so the
# reference category is unchanged (first level = reference)
factor_levels <- read_table(genr_file(cfg, "factor_levels"))
for (v in cfg$categorical_covariates) {
  lv <- factor_levels[variable == v][order(order), level]
  genr[, (v) := if (length(lv)) droplevels(factor(get(v), levels = lv)) else factor(get(v))]
}

# ---- Genetic scores (children without genotype data get NA) -------------------
for (path in cfg$score_files) add_columns_by_idc(genr, read_table(path), path)
require_columns(genr, all_scores, "The genetic score file(s)",
                load_command = sprintf('d <- data.table::fread("%s")', cfg$score_files[1]))

message(sprintf("Loaded %d children, %d scores", nrow(genr), length(all_scores)))


# =============================================================================
# 2. Repeated measures in long format
# =============================================================================
# One row per child x wave x outcome, keeping observations with both an outcome
# value and an age. Columns: idc, covariates, age, outcome, domain,
# outcome_name, wave.

id_and_covariates <- c("idc", cfg$covariates)

long_data <- rbindlist(lapply(seq_len(nrow(cfg$longitudinal)), function(k) {
  m <- cfg$longitudinal[k]
  age_col <- paste0(m$age_prefix, "_", m$wave)
  outcome_col <- paste0(m$column, "_", m$wave)

  d <- genr[, c(id_and_covariates, age_col, outcome_col), with = FALSE]
  setnames(d, c(age_col, outcome_col), c("age", "outcome"))
  d[, `:=`(age = as.numeric(age), outcome = as.numeric(outcome),
           domain = m$domain, outcome_name = m$outcome, wave = m$wave)]
  d[!is.na(outcome) & !is.na(age)]
}))

long_outcomes <- unique(cfg$longitudinal[, .(domain, outcome)])


# =============================================================================
# 3. Model helpers
# =============================================================================

# Covariates that vary in the analysed sample; a constant covariate (e.g. a
# single genotyping cohort) cannot be estimated and is left out.
varying_covariates <- function(d) {
  covars <- c("sex", "genotype_cohort", gpcs)
  covars[sapply(covars, function(v) uniqueN(d[[v]]) > 1)]
}

# Standardise the genetic PCs over children
standardise_gpcs <- function(d) {
  for (v in gpcs) d[, (v) := z_by_id(get(v), idc)]
  d
}

# Marginal R2 (Nakagawa): variance explained by the fixed effects.
# NA for singular fits (random-intercept variance estimated as ~0).
marginal_r2 <- function(model) {
  if (isSingular(model, tol = 1e-4)) return(NA_real_)
  r2 <- tryCatch(
    suppressWarnings(performance::r2_nakagawa(model, tolerance = 1e-12)$R2_marginal),
    error = function(e) NA_real_
  )
  if (!length(r2)) return(NA_real_)
  as.numeric(r2[1])
}

# Name of the fixed effect for the interaction a:b (either order)
interaction_term <- function(model, a, b) {
  terms <- names(fixef(model))
  hit <- terms[vapply(strsplit(terms, ":", fixed = TRUE),
                      function(x) length(x) == 2 && setequal(x, c(a, b)), logical(1))]
  if (length(hit) != 1) stop("Could not identify interaction: ", a, " x ", b)
  hit
}

# Runs one model fit. Errors are logged in `failures` and return NULL (the
# loop continues); warnings are logged in `model_warnings`.
failures <- list()
model_warnings <- list()

safe_fit <- function(expr, label) {
  withCallingHandlers(
    tryCatch(expr, error = function(e) {
      failures[[length(failures) + 1]] <<- data.table(model = label, error = conditionMessage(e))
      NULL
    }),
    warning = function(w) {
      model_warnings[[length(model_warnings) + 1]] <<- data.table(model = label, warning = conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
}


# =============================================================================
# 4. Linear mixed model: repeated measures, one score at a time
# =============================================================================
# Data: one row per observation j (wave) of child i.
#
#   M0: y_ij = b0 + b1 a_ij + b2 a_ij^2 + g'X_i                              + u_i + e_ij
#   M1: y_ij = b0 + b1 a_ij + b2 a_ij^2 + g'X_i + b3 S_i                     + u_i + e_ij
#   M2: y_ij = b0 + b1 a_ij + b2 a_ij^2 + g'X_i + b3 S_i
#                 + b4 S_i a_ij + b5 S_i a_ij^2                              + u_i + e_ij
#
#   y_ij  outcome, z-scored over all observations (all waves pooled)
#   a_ij  age - 17 (cfg$reference_age), in years; a_ij^2 its square
#   S_i   genetic score, z-scored over children (one value per child)
#   X_i   sex, genotyping cohort, gPC1-4 (gPCs z-scored over children);
#         covariates that are constant in the sample are dropped
#   u_i   random intercept of child i, u_i ~ N(0, s_u^2)
#   e_ij  residual, e_ij ~ N(0, s_e^2)
#
# Fixed effects
#   - age and age^2: average (quadratic) developmental trajectory of the outcome
#   - sex: sex differences in the outcome
#   - genotyping cohort: differences between genotyping batches/arrays
#   - gPC1-4: population stratification (ancestry differences that could
#     correlate with both the score and the outcome)
#   - S (M1, M2): the genetic score being tested
#   - S x age, S x age^2 (M2): lets the score effect change with age
# Random effects
#   - a child-specific intercept u_i. Repeated observations of the same child
#     are correlated; u_i absorbs stable between-child differences so the
#     standard errors account for this non-independence. There are no random
#     slopes: apart from the score terms, all children follow the same trajectory.
# Estimation: maximum likelihood (REML = FALSE).
#
# Tests
#   H1: b3 = 0 in M1 (t-test, Satterthwaite df). M1 has no interaction, so b3 is
#       the association averaged over ages, in SD of outcome per SD of score.
#   H2: b4 = b5 = 0 in M2 (joint F-test, Satterthwaite df): does the association
#       change with age? In M2, b3 is the association at age 17.
#   Incremental R2 (Nakagawa marginal): M1 - M0 for H1, M2 - M1 for H2.

fit_lmm_one_score <- function(data, score_name) {

  d <- copy(data)
  d[, score := genr[[score_name]][match(idc, genr$idc)]]

  model_vars <- c("idc", "outcome", "age", "sex", "genotype_cohort", gpcs, "score")
  d <- d[complete.cases(d[, ..model_vars])]

  # Too little variation to fit the model
  if (!nrow(d) || uniqueN(d$outcome) < 2 || uniqueN(d$score) < 2 || uniqueN(d$age) < 3) {
    return(NULL)
  }

  d[, outcome := z(outcome)]
  d[, score := z_by_id(score, idc)]
  d <- standardise_gpcs(d)
  d[, `:=`(age_c = age - cfg$reference_age, age2 = (age - cfg$reference_age)^2)]

  covars <- varying_covariates(d)
  random_intercept <- "(1 | idc)"

  m0 <- lmer(reformulate(c("age_c", "age2", covars, random_intercept), "outcome"),
             data = d, REML = FALSE)
  m1 <- lmer(reformulate(c("age_c", "age2", "score", covars, random_intercept), "outcome"),
             data = d, REML = FALSE)
  m2 <- lmer(reformulate(c("age_c", "age2", "score", "score:age_c", "score:age2", covars, random_intercept),
                         "outcome"),
             data = d, REML = FALSE)

  r2_m0 <- marginal_r2(m0)
  r2_m1 <- marginal_r2(m1)
  r2_m2 <- marginal_r2(m2)

  coef_m1 <- summary(m1)$coefficients
  coef_m2 <- summary(m2)$coefficients

  # H2: joint F-test of the two score x age terms
  age_term  <- interaction_term(m2, "score", "age_c")
  age2_term <- interaction_term(m2, "score", "age2")
  contrast <- matrix(0, 2, length(fixef(m2)), dimnames = list(NULL, names(fixef(m2))))
  contrast[1, age_term]  <- 1
  contrast[2, age2_term] <- 1
  h2_test <- contestMD(m2, contrast, ddf = "Satterthwaite")

  n_id <- uniqueN(d$idc)

  main <- data.table(
    hypothesis     = c("H1", "H2"),
    estimate       = c(coef_m1["score", "Estimate"], NA_real_),
    std_error      = c(coef_m1["score", "Std. Error"], NA_real_),
    p_value        = c(coef_m1["score", "Pr(>|t|)"], h2_test[1, "Pr(>F)"]),
    incremental_r2 = c(r2_m1 - r2_m0, r2_m2 - r2_m1),
    predictor      = score_name,
    n_id           = n_id,
    n_obs          = nrow(d),
    singular       = c(isSingular(m1), isSingular(m2))
  )

  h2_details <- data.table(
    predictor          = score_name,
    age_interaction    = coef_m2[age_term, "Estimate"],
    age_interaction_se = coef_m2[age_term, "Std. Error"],
    age_interaction_p  = coef_m2[age_term, "Pr(>|t|)"],
    age2_interaction    = coef_m2[age2_term, "Estimate"],
    age2_interaction_se = coef_m2[age2_term, "Std. Error"],
    age2_interaction_p  = coef_m2[age2_term, "Pr(>|t|)"],
    joint_F        = h2_test[1, "F value"],
    joint_p        = h2_test[1, "Pr(>F)"],
    incremental_r2 = r2_m2 - r2_m1,
    n_id           = n_id,
    n_obs          = nrow(d)
  )

  covariate_terms <- grep("^(sex|genotype_cohort|gPC)", rownames(coef_m2), value = TRUE)
  covariates <- data.table(
    term      = covariate_terms,
    estimate  = as.numeric(coef_m2[covariate_terms, "Estimate"]),
    p_value   = as.numeric(coef_m2[covariate_terms, "Pr(>|t|)"]),
    predictor = score_name
  )

  list(main = main, h2_details = h2_details, covariates = covariates)
}


# =============================================================================
# 5. Linear mixed model: repeated measures, all scores of a set jointly
# =============================================================================
# Same data preparation, covariates and random intercept as section 4, with
# all K scores of the set entered together (each z-scored over children):
#
#   M0: y_ij = b0 + b1 a_ij + b2 a_ij^2 + g'X_i                      + u_i + e_ij
#   M1: y_ij = b0 + b1 a_ij + b2 a_ij^2 + g'X_i + sum_k c_k S_ik     + u_i + e_ij
#
# H1: c_1 = ... = c_K = 0, likelihood-ratio test of M1 vs M0 (chi-square, K df;
#     this is why both models are fitted with maximum likelihood).
# No per-score estimate is reported. Incremental R2 = R2(M1) - R2(M0).

fit_lmm_score_set <- function(data, scores, label) {

  d <- copy(data)
  rows <- match(d$idc, genr$idc)
  d[, (scores) := genr[rows, scores, with = FALSE]]

  model_vars <- c("idc", "outcome", "age", "sex", "genotype_cohort", gpcs, scores)
  d <- d[complete.cases(d[, ..model_vars])]

  if (!nrow(d) || uniqueN(d$outcome) < 2 || uniqueN(d$age) < 3) return(NULL)
  if (any(sapply(scores, function(v) uniqueN(d[[v]]) < 2))) return(NULL)

  d[, outcome := z(outcome)]
  for (s in scores) d[, (s) := z_by_id(get(s), idc)]
  d <- standardise_gpcs(d)
  d[, `:=`(age_c = age - cfg$reference_age, age2 = (age - cfg$reference_age)^2)]

  covars <- varying_covariates(d)
  random_intercept <- "(1 | idc)"

  m0 <- lmer(reformulate(c("age_c", "age2", covars, random_intercept), "outcome"),
             data = d, REML = FALSE)
  m1 <- lmer(reformulate(c("age_c", "age2", scores, covars, random_intercept), "outcome"),
             data = d, REML = FALSE)

  lr_test <- anova(m0, m1)

  data.table(
    hypothesis     = "H1",
    estimate       = NA_real_,
    std_error      = NA_real_,
    p_value        = lr_test[2, "Pr(>Chisq)"],
    incremental_r2 = marginal_r2(m1) - marginal_r2(m0),
    predictor      = label,
    n_id           = uniqueN(d$idc),
    n_obs          = nrow(d),
    singular       = isSingular(m1)
  )
}


# =============================================================================
# 6. Linear regression: outcomes measured once
# =============================================================================
# One row per child i:
#
#   M0: y_i = b0 + g'X_i                  + e_i
#   M1: y_i = b0 + g'X_i + sum_k c_k S_ik + e_i,   e_i ~ N(0, s^2)
#
#   y_i   outcome, z-scored
#   S_ik  genetic score(s), z-scored (one score, or all scores of a set)
#   X_i   sex, genotyping cohort, gPC1-4 (z-scored), plus the age at
#         measurement when the outcome has one (cfg$single_outcomes$age_var);
#         constant covariates are dropped
#
# Fitted with ordinary least squares. No random effects: one observation per child.
# H1: all c_k = 0, F-test of M1 vs M0. For a single score, the estimate is c_1
#     (SD of outcome per SD of score). Incremental R2 = R2(M1) - R2(M0).

fit_lm <- function(outcome, scores, label, combined = FALSE, age_var = NA_character_) {

  model_vars <- c("idc", outcome, "sex", "genotype_cohort", gpcs, scores)
  if (!is.na(age_var)) model_vars <- c(model_vars, age_var)

  d <- copy(genr[, ..model_vars])
  d <- d[complete.cases(d)]

  if (!nrow(d) || uniqueN(d[[outcome]]) < 2) return(NULL)
  if (any(sapply(scores, function(v) uniqueN(d[[v]]) < 2))) return(NULL)

  d[, (outcome) := z(get(outcome))]
  for (s in scores) d[, (s) := z(get(s))]
  d <- standardise_gpcs(d)

  covars <- varying_covariates(d)
  if (!is.na(age_var) && uniqueN(d[[age_var]]) > 1) covars <- c(age_var, covars)

  m0 <- lm(reformulate(covars, outcome), data = d)
  m1 <- lm(reformulate(c(covars, scores), outcome), data = d)

  f_test <- anova(m0, m1)
  coefs <- summary(m1)$coefficients

  data.table(
    hypothesis     = "H1",
    estimate       = if (combined) NA_real_ else coefs[label, "Estimate"],
    std_error      = if (combined) NA_real_ else coefs[label, "Std. Error"],
    p_value        = f_test[2, "Pr(>F)"],
    incremental_r2 = summary(m1)$r.squared - summary(m0)$r.squared,
    predictor      = label,
    n_id           = nrow(d),
    n_obs          = nrow(d),
    singular       = NA
  )
}


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
