# =============================================================================
# Model definitions shared by step 5 (all models) and step 7 (trajectory figures)
# =============================================================================
# Sourced after `cfg`, `genr` (individual-level data with scores) and `gpcs`
# exist. The functions below use these objects.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)     # Satterthwaite t- and F-tests for lmer models
  library(performance)  # Nakagawa R2
})

# =============================================================================
# 3. Model helpers
# =============================================================================

# Covariates that vary in the analysed sample; a constant covariate (e.g. a
# single genotyping cohort) cannot be estimated and is left out.
varying_covariates <- function(d) {
  covars <- c("sex", "genotype_cohort", gpcs)
  covars[sapply(covars, function(v) uniqueN(d[[v]]) > 1)]
}

# Standardise the genetic PCs over individuals
#
# DANGER! This line is a GUESS. Line 134 of the original regression script
# (inside scale_gpcs()) was not visible in the screenshots we received. We
# assumed it is  d[, (v) := z_by_id(get(v), idc)] , which standardises each gPC
# over individuals, matching how the scores are standardised. If the real line is
# plain  z(get(v))  (standardised over all rows), only the gPC rows of
# covariate_results.csv from the mixed models would change. Rescaling a
# covariate does not change any score estimate, p-value or R2.
# Check against the original script and update if it differs.
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

# Joint F-test that the fixed effects `terms` are all zero, with Satterthwaite
# denominator df (lmerTest::contestMD). Used for every mixed-model hypothesis
# test, so all mixed-model p-values come from the same kind of test. With one
# term, F = t^2 and the p-value equals that of the Satterthwaite t-test.
joint_f_test <- function(model, terms) {
  contrast <- matrix(0, length(terms), length(fixef(model)),
                     dimnames = list(terms, names(fixef(model))))
  contrast[cbind(terms, terms)] <- 1
  contestMD(model, contrast, ddf = "Satterthwaite")
}

# Adds centred age (age_c) and its square (age2). The original code centres at
# cfg$reference_age (17). Centring does not change H1, H2 or R2; it changes the
# meaning and precision of the single M2 coefficients (see README, concern C1).
centre_age <- function(d) {
  centre <- cfg$reference_age
  # Uncomment to centre on the average age of the model's analysis sample instead:
  # centre <- mean(d$age)
  d[, `:=`(age_c = age - centre, age2 = (age - centre)^2)]
  d
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
# Data: one row per observation j (wave) of individual i.
#
#   M0: y_ij = b0 + b1 a_ij + b2 a_ij^2 + g'X_i                              + u_i + e_ij
#   M1: y_ij = b0 + b1 a_ij + b2 a_ij^2 + g'X_i + b3 S_i                     + u_i + e_ij
#   M2_lin: M1 + b4 S_i a_ij                                                 (linear interaction)
#   M2: y_ij = b0 + b1 a_ij + b2 a_ij^2 + g'X_i + b3 S_i
#                 + b4 S_i a_ij + b5 S_i a_ij^2                              + u_i + e_ij
#
#   y_ij  outcome, z-scored over all observations (all waves pooled)
#   a_ij  age - 17 (cfg$reference_age), in years; a_ij^2 its square
#   S_i   genetic score: a GPS (GenomICA polygenic score) or PGS, z-scored over
#         individuals (one value per individual)
#   X_i   sex, genotyping cohort, gPC1-4 (gPCs z-scored over individuals);
#         covariates that are constant in the sample are dropped
#   u_i   random intercept of individual i, u_i ~ N(0, s_u^2)
#   e_ij  residual, e_ij ~ N(0, s_e^2)
#
# Fixed effects
#   - age and age^2: average (quadratic) developmental trajectory of the outcome
#   - sex: sex differences in the outcome
#   - genotyping cohort: differences between genotyping batches/arrays
#   - gPC1-4: population stratification (ancestry differences that could
#     correlate with both the score and the outcome)
#   - S (M1, M2_lin, M2): the genetic score being tested
#   - S x age (M2_lin, M2), S x age^2 (M2): lets the score effect change with age
# Random effects
#   - an individual-specific intercept u_i. Repeated observations of the same individual
#     are correlated; u_i absorbs stable between-individual differences so the
#     standard errors account for this non-independence. There are no random
#     slopes: apart from the score terms, all individuals follow the same trajectory.
# Estimation: maximum likelihood (REML = FALSE).
#
# Tests
#   H1: b3 = 0 in M1 (F-test with 1 df, Satterthwaite; identical p-value to the
#       Satterthwaite t-test). M1 has no interaction, so b3 is the association
#       averaged over ages, in SD of outcome per SD of score.
#   H2_lin_interact: b4 = 0 in M2_lin (F-test with 1 df, Satterthwaite): does the
#       association change linearly with age? b4 does not depend on where age
#       is centred.
#   H2_quad_interact: b5 = 0 in M2 (F-test with 1 df, Satterthwaite): does the
#       change with age curve, on top of the linear change? b5 does not depend
#       on where age is centred.
#   H2_lin_and_quad_interact: b4 = b5 = 0 in M2 (joint F-test, 2 df, Satterthwaite):
#       does the association change with age, linearly or curved? This is the
#       original code's H2. In M2, b3 is the association at age 17.
#   Incremental R2 (Nakagawa marginal): M1 - M0 for H1, M2_lin - M1 for
#       H2_lin_interact, M2 - M2_lin for H2_quad_interact, M2 - M1 for
#       H2_lin_and_quad_interact.

fit_lmm_one_score <- function(data, score_name, return_model = FALSE) {

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
  d <- centre_age(d)

  covars <- varying_covariates(d)
  random_intercept <- "(1 | idc)"

  m0 <- lmer(reformulate(c("age_c", "age2", covars, random_intercept), "outcome"),
             data = d, REML = FALSE)
  m1 <- lmer(reformulate(c("age_c", "age2", "score", covars, random_intercept), "outcome"),
             data = d, REML = FALSE)
  m2_lin <- lmer(reformulate(c("age_c", "age2", "score", "score:age_c", covars, random_intercept), "outcome"),
                 data = d, REML = FALSE)
  m2 <- lmer(reformulate(c("age_c", "age2", "score", "score:age_c", "score:age2", covars, random_intercept),
                         "outcome"),
             data = d, REML = FALSE)

  r2_m0 <- marginal_r2(m0)
  r2_m1 <- marginal_r2(m1)
  r2_m2_lin <- marginal_r2(m2_lin)
  r2_m2 <- marginal_r2(m2)

  coef_m1 <- summary(m1)$coefficients
  coef_m2_lin <- summary(m2_lin)$coefficients
  coef_m2 <- summary(m2)$coefficients

  # H1: F-test of the score effect in M1
  h1_test <- joint_f_test(m1, "score")
  h1_p <- h1_test[1, "Pr(>F)"]
  # The old code used the t-test p-value from summary(m1) (same value):
  # h1_p <- coef_m1["score", "Pr(>|t|)"]

  # H2_lin_interact: F-test (1 df) of the linear score x age term in M2_lin
  lin_term <- interaction_term(m2_lin, "score", "age_c")
  h2_lin_test <- joint_f_test(m2_lin, lin_term)

  # H2_quad_interact: F-test (1 df) of the quadratic score x age^2 term in M2,
  # on top of the linear term
  age_term  <- interaction_term(m2, "score", "age_c")
  age2_term <- interaction_term(m2, "score", "age2")
  h2_quad_test <- joint_f_test(m2, age2_term)

  # H2_lin_and_quad_interact: joint F-test (2 df) of both score x age terms in M2
  # (the original code's H2)
  h2_test <- joint_f_test(m2, c(age_term, age2_term))

  n_id <- uniqueN(d$idc)

  main <- data.table(
    hypothesis     = c("H1", "H2_lin_interact", "H2_quad_interact", "H2_lin_and_quad_interact"),
    estimate       = c(coef_m1["score", "Estimate"], coef_m2_lin[lin_term, "Estimate"],
                       coef_m2[age2_term, "Estimate"], NA_real_),
    std_error      = c(coef_m1["score", "Std. Error"], coef_m2_lin[lin_term, "Std. Error"],
                       coef_m2[age2_term, "Std. Error"], NA_real_),
    p_value        = c(h1_p, h2_lin_test[1, "Pr(>F)"], h2_quad_test[1, "Pr(>F)"], h2_test[1, "Pr(>F)"]),
    incremental_r2 = c(r2_m1 - r2_m0, r2_m2_lin - r2_m1, r2_m2 - r2_m2_lin, r2_m2 - r2_m1),
    predictor      = score_name,
    n_id           = n_id,
    n_obs          = nrow(d),
    singular       = c(isSingular(m1), isSingular(m2_lin), isSingular(m2), isSingular(m2))
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
    lin_interaction    = coef_m2_lin[lin_term, "Estimate"],     # from M2_lin
    lin_interaction_se = coef_m2_lin[lin_term, "Std. Error"],
    lin_interaction_p  = h2_lin_test[1, "Pr(>F)"],
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

  result <- list(main = main, h2_details = h2_details, covariates = covariates)
  # Step 7 (figures) uses the fitted models and their analysis data
  if (return_model) result$fit <- list(m1 = m1, m2_lin = m2_lin, m2 = m2, data = d)
  result
}


# =============================================================================
# 5. Linear mixed model: repeated measures, all scores of a set jointly
# =============================================================================
# Same data preparation, covariates and random intercept as section 4, with
# all K scores of the set entered together (each z-scored over individuals):
#
#   M0: y_ij = b0 + b1 a_ij + b2 a_ij^2 + g'X_i                      + u_i + e_ij
#   M1: y_ij = b0 + b1 a_ij + b2 a_ij^2 + g'X_i + sum_k c_k S_ik     + u_i + e_ij
#
# H1: c_1 = ... = c_K = 0, joint F-test in M1 (K numerator df, Satterthwaite
#     denominator df), the same test as H2 in section 4. The original code used
#     a likelihood-ratio test of M1 vs M0 instead (see below).
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
  d <- centre_age(d)

  covars <- varying_covariates(d)
  random_intercept <- "(1 | idc)"

  m0 <- lmer(reformulate(c("age_c", "age2", covars, random_intercept), "outcome"),
             data = d, REML = FALSE)
  m1 <- lmer(reformulate(c("age_c", "age2", scores, covars, random_intercept), "outcome"),
             data = d, REML = FALSE)

  h1_p <- joint_f_test(m1, scores)[1, "Pr(>F)"]
  # The old code used a likelihood-ratio chi-square test of M1 vs M0 (K df);
  # to reproduce the original results, use this line instead:
  # h1_p <- anova(m0, m1)[2, "Pr(>Chisq)"]

  data.table(
    hypothesis     = "H1",
    estimate       = NA_real_,
    std_error      = NA_real_,
    p_value        = h1_p,
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
# One row per individual i:
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
# Fitted with ordinary least squares. No random effects: one observation per individual.
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


