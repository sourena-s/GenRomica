# Original "Regression code", transcribed from Code_Hailey_Components.docx
# (screenshots of lines 1-476). Used only by toy/compare_with_original.R.
# Edits, all marked "# EDIT":
#   - data path and score-set sizes, so it runs on the toy data
#   - line 134 (body of scale_gpcs) was not visible: assumed z_by_id per gPC
#   - lines after 476 (end of merge + writing output) were not visible: completed

library(data.table)
library(lme4)
library(lmerTest)
library(performance)

# ============================================================
# Loading data
# ============================================================

data_dir <- Sys.getenv("ORIGINAL_DATA_DIR")                                  # EDIT (was "/mnt/data/Hailey/data")
full_data <- as.data.table(readRDS(file.path(data_dir, "genr_full_file.rds")))

waves <- c("f05", "f09", "f13", "f17")
reference_age <- 17
gpcs <- paste0("gPC", 1:4)

score_sets <- list(
  PCS = paste0("C", 1:as.integer(Sys.getenv("TOY_N_PCS", "20"))),            # EDIT (was 1:20)
  UK_PCS = paste0("UK_PCS_", 1:as.integer(Sys.getenv("TOY_N_UKPCS", "20"))), # EDIT (was 1:20)
  PGS = paste0("PGS_", 1:as.integer(Sys.getenv("TOY_N_PGS", "18")))          # EDIT (was 1:18)
)

combined_sets <- list(
  PCS_combined = score_sets$PCS,
  UK_PCS_combined = score_sets$UK_PCS
)

single_outcomes <- data.table(
  outcome = c("INCOME5", "EDUCM", "EDUCP", "DietScore_8y", "onlineweek_child_F17"),
  domain = c("SES", "SES", "SES", "Diet", "Online"),
  age_var = c(NA_character_, NA_character_, NA_character_, NA_character_, "AgeChild_GR1203")
)

all_scores <- unlist(score_sets, use.names = FALSE)
base_vars <- c("idc", "sex", "genotype_cohort", gpcs)

full_data[, sex := factor(sex)]
full_data[, genotype_cohort := factor(genotype_cohort)]

if ("sum_wid_f05" %in% names(full_data) && !"sum_wit_f05" %in% names(full_data)) {
  setnames(full_data, "sum_wid_f05", "sum_wit_f05")
}

required <- c(
  base_vars, all_scores, single_outcomes$outcome,
  paste0("age_cbcl_", waves), paste0("age_child_mri_", waves),
  "AgeChild_GR1203"
)

# ============================================================
# Longitudinal data
# ============================================================

make_long <- function(outcomes, age_prefix, domain) {

  rbindlist(lapply(waves, function(wave) {

    age_var <- paste0(age_prefix, "_", wave)
    outcome_vars <- paste0(outcomes, "_", wave)
    outcome_vars <- outcome_vars[outcome_vars %in% names(full_data)]

    if (!length(outcome_vars)) return(NULL)

    x <- copy(full_data[, c(base_vars, age_var, outcome_vars), with = FALSE])
    setnames(x, age_var, "age")

    x[, age := as.numeric(age)]
    x[, (outcome_vars) := lapply(.SD, as.numeric), .SDcols = outcome_vars]

    x <- melt(
      x,
      id.vars = c(base_vars, "age"),
      measure.vars = outcome_vars,
      variable.name = "outcome_name",
      value.name = "outcome",
      variable.factor = FALSE,
      na.rm = TRUE
    )

    x <- x[!is.na(age)]
    x[, outcome_name := sub(paste0("_", wave, "$"), "", outcome_name)]
    if (domain == "CBCL") x[, outcome_name := sub("^sum_", "", outcome_name)]
    x[, `:=`(wave = wave, domain = domain)]
    x
  }), fill = TRUE)
}

cbcl_long <- make_long(
  paste0("sum_", c("anx", "som", "wit", "sop", "tho", "att", "rul", "agg")),
  "age_cbcl",
  "CBCL"
)

mri_long <- make_long(
  c("TotalGrayVol", "cortical_thickness"),
  "age_child_mri",
  "MRI"
)

dti_long <- make_long(
  c("mean_FA_genr", "mean_MD_genr"),
  "age_child_mri",
  "DTI"
)

long_data <- rbindlist(list(cbcl_long, mri_long, dti_long), fill = TRUE)
long_outcomes <- unique(long_data[, .(domain, outcome_name)])

# ============================================================
# HELPERS
# ============================================================

z <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / s
}

z_by_id <- function(x, id) {
  keep <- !duplicated(id)
  m <- mean(x[keep])
  s <- sd(x[keep])
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

get_covars <- function(d) {
  covars <- c("sex", "genotype_cohort", gpcs)
  covars[sapply(covars, function(v) uniqueN(d[[v]]) > 1)]
}

scale_gpcs <- function(d) {
  for (v in gpcs) {
    d[, (v) := z_by_id(get(v), idc)]                                           # EDIT (line 134 not visible; assumed)
  }
  d
}

get_r2 <- function(model) {
  if (isSingular(model, tol = 1e-4)) return(NA_real_)
  x <- tryCatch(
    suppressWarnings(performance::r2_nakagawa(model, tolerance = 1e-12)$R2_marginal),
    error = function(e) NA_real_
  )
  if (!length(x)) return(NA_real_)
  as.numeric(x[1])
}

interaction_term <- function(model, a, b) {
  nms <- names(fixef(model))
  hit <- nms[vapply(strsplit(nms, ":", fixed = TRUE), function(x) {
    length(x) == 2 && setequal(x, c(a, b))
  }, logical(1))]
  if (length(hit) != 1) stop("Could not identify interaction: ", a, " x ", b)
  hit
}

failures <- list()

safe_fit <- function(expr, label) {
  tryCatch(
    expr,
    error = function(e) {
      failures[[length(failures) + 1]] <<- data.table(
        model = label,
        error = conditionMessage(e)
      )
      NULL
    }
  )
}

# ============================================================
# Longitudinal model
# ============================================================

fit_lmm <- function(data, predictor) {

  d <- copy(data)
  d[, score := full_data[[predictor]][match(idc, full_data$idc)]]

  vars <- c("idc", "outcome", "age", "sex", "genotype_cohort", gpcs, "score")
  d <- d[complete.cases(d[, ..vars])]

  if (!nrow(d) || uniqueN(d$outcome) < 2 || uniqueN(d$score) < 2 || uniqueN(d$age) < 3) {
    return(NULL)
  }

  d[, outcome := z(outcome)]
  d[, score := z_by_id(score, idc)]
  d <- scale_gpcs(d)
  d[, `:=`(age_c = age - reference_age, age2 = (age - reference_age)^2)]

  covars <- get_covars(d)

  m0 <- lmer(
    reformulate(c("age_c", "age2", covars, "(1 | idc)"), "outcome"),
    data = d, REML = FALSE
  )

  m1 <- lmer(
    reformulate(c("age_c", "age2", "score", covars, "(1 | idc)"), "outcome"),
    data = d, REML = FALSE
  )

  m2 <- lmer(
    reformulate(
      c("age_c", "age2", "score", "score:age_c", "score:age2", covars, "(1 | idc)"),
      "outcome"
    ),
    data = d,
    REML = FALSE
  )

  r0 <- get_r2(m0)
  r1 <- get_r2(m1)
  r2 <- get_r2(m2)

  co1 <- summary(m1)$coefficients
  co2 <- summary(m2)$coefficients

  age_term <- interaction_term(m2, "score", "age_c")
  age2_term <- interaction_term(m2, "score", "age2")

  L <- matrix(0, 2, length(fixef(m2)), dimnames = list(NULL, names(fixef(m2))))
  L[1, age_term] <- 1
  L[2, age2_term] <- 1

  h2 <- contestMD(m2, L, ddf = "Satterthwaite")

  main <- data.table(
    hypothesis = c("H1", "H2"),
    estimate = c(co1["score", "Estimate"], NA_real_),
    std_error = c(co1["score", "Std. Error"], NA_real_),
    p_value = c(co1["score", "Pr(>|t|)"], h2[1, "Pr(>F)"]),
    incremental_r2 = c(r1 - r0, r2 - r1),
    predictor = predictor,
    n_id = uniqueN(d$idc),
    n_obs = nrow(d),
    singular = c(isSingular(m1), isSingular(m2))
  )

  h2_details <- data.table(
    predictor = predictor,
    age_interaction = co2[age_term, "Estimate"],
    age_interaction_se = co2[age_term, "Std. Error"],
    age_interaction_p = co2[age_term, "Pr(>|t|)"],
    age2_interaction = co2[age2_term, "Estimate"],
    age2_interaction_se = co2[age2_term, "Std. Error"],
    age2_interaction_p = co2[age2_term, "Pr(>|t|)"],
    joint_F = h2[1, "F value"],
    joint_p = h2[1, "Pr(>F)"],
    incremental_r2 = r2 - r1,
    n_id = uniqueN(d$idc),
    n_obs = nrow(d)
  )

  cov_names <- grep("^(sex|genotype_cohort|gPC)", rownames(co2), value = TRUE)

  covariates <- data.table(
    term = cov_names,
    estimate = as.numeric(co2[cov_names, "Estimate"]),
    p_value = as.numeric(co2[cov_names, "Pr(>|t|)"]),
    predictor = predictor
  )

  list(main = main, h2_details = h2_details, covariates = covariates)
}

# ============================================================
# COMBINED LONGITUDINAL MODEL
# ============================================================

fit_combined_lmm <- function(data, genes, label) {

  d <- copy(data)
  idx <- match(d$idc, full_data$idc)
  d[, (genes) := full_data[idx, ..genes]]

  vars <- c("idc", "outcome", "age", "sex", "genotype_cohort", gpcs, genes)
  d <- d[complete.cases(d[, ..vars])]

  if (!nrow(d) || uniqueN(d$outcome) < 2 || uniqueN(d$age) < 3) return(NULL)
  if (any(sapply(genes, function(v) uniqueN(d[[v]]) < 2))) return(NULL)

  d[, outcome := z(outcome)]
  for (gene in genes) d[, (gene) := z_by_id(get(gene), idc)]
  d <- scale_gpcs(d)
  d[, `:=`(age_c = age - reference_age, age2 = (age - reference_age)^2)]

  covars <- get_covars(d)

  m0 <- lmer(
    reformulate(c("age_c", "age2", covars, "(1 | idc)"), "outcome"),
    data = d, REML = FALSE
  )

  m1 <- lmer(
    reformulate(c("age_c", "age2", genes, covars, "(1 | idc)"), "outcome"),
    data = d, REML = FALSE
  )

  test <- anova(m0, m1)

  data.table(
    hypothesis = "H1",
    estimate = NA_real_,
    std_error = NA_real_,
    p_value = test[2, "Pr(>Chisq)"],
    incremental_r2 = get_r2(m1) - get_r2(m0),
    predictor = label,
    n_id = uniqueN(d$idc),
    n_obs = nrow(d),
    singular = isSingular(m1)
  )
}

# ============================================================
# SINGLE-TIME MODEL
# ============================================================

fit_single <- function(outcome, genes, label, combined = FALSE, age_var = NA_character_) {

  vars <- c("idc", outcome, "sex", "genotype_cohort", gpcs, genes)
  if (!is.na(age_var)) vars <- c(vars, age_var)

  d <- copy(full_data[, ..vars])
  d <- d[complete.cases(d)]

  if (!nrow(d) || uniqueN(d[[outcome]]) < 2) return(NULL)
  if (any(sapply(genes, function(v) uniqueN(d[[v]]) < 2))) return(NULL)

  d[, (outcome) := z(get(outcome))]
  for (gene in genes) d[, (gene) := z(get(gene))]
  d <- scale_gpcs(d)

  covars <- get_covars(d)
  if (!is.na(age_var) && uniqueN(d[[age_var]]) > 1) covars <- c(age_var, covars)

  m0 <- lm(reformulate(covars, outcome), data = d)
  m1 <- lm(reformulate(c(covars, genes), outcome), data = d)

  test <- anova(m0, m1)
  co <- summary(m1)$coefficients

  data.table(
    hypothesis = "H1",
    estimate = if (combined) NA_real_ else co[label, "Estimate"],
    std_error = if (combined) NA_real_ else co[label, "Std. Error"],
    p_value = test[2, "Pr(>F)"],
    incremental_r2 = summary(m1)$r.squared - summary(m0)$r.squared,
    predictor = label,
    n_id = nrow(d),
    n_obs = nrow(d),
    singular = NA
  )
}

# ============================================================
# Run models
# ============================================================

all_results <- list()
all_h2 <- list()
all_covariates <- list()

for (score_group in names(score_sets)) {

  for (predictor in score_sets[[score_group]]) {

    message("Running ", score_group, ": ", predictor)

    for (i in seq_len(nrow(long_outcomes))) {

      domain_i <- long_outcomes$domain[i]
      outcome_i <- long_outcomes$outcome_name[i]
      x <- long_data[domain == domain_i & outcome_name == outcome_i]

      fit <- safe_fit(
        fit_lmm(x, predictor),
        paste(score_group, predictor, domain_i, outcome_i, sep = " | ")
      )

      if (is.null(fit)) next

      fit$main[, `:=`(score_group = score_group, domain = domain_i, outcome = outcome_i)]
      fit$h2_details[, `:=`(score_group = score_group, domain = domain_i, outcome = outcome_i)]
      fit$covariates[, `:=`(score_group = score_group, domain = domain_i, outcome = outcome_i)]

      all_results[[length(all_results) + 1]] <- fit$main
      all_h2[[length(all_h2) + 1]] <- fit$h2_details
      all_covariates[[length(all_covariates) + 1]] <- fit$covariates
    }

    for (i in seq_len(nrow(single_outcomes))) {

      outcome_i <- single_outcomes$outcome[i]
      domain_i <- single_outcomes$domain[i]
      age_i <- single_outcomes$age_var[i]

      fit <- safe_fit(
        fit_single(outcome_i, predictor, predictor, FALSE, age_i),
        paste(score_group, predictor, domain_i, outcome_i, sep = " | ")
      )

      if (is.null(fit)) next

      fit[, `:=`(score_group = score_group, domain = domain_i, outcome = outcome_i)]
      all_results[[length(all_results) + 1]] <- fit
    }
  }
}

for (score_group in names(combined_sets)) {

  genes <- combined_sets[[score_group]]
  message("Running ", score_group)

  for (i in seq_len(nrow(long_outcomes))) {

    domain_i <- long_outcomes$domain[i]
    outcome_i <- long_outcomes$outcome_name[i]
    x <- long_data[domain == domain_i & outcome_name == outcome_i]

    fit <- safe_fit(
      fit_combined_lmm(x, genes, score_group),
      paste(score_group, domain_i, outcome_i, sep = " | ")
    )

    if (is.null(fit)) next

    fit[, `:=`(score_group = score_group, domain = domain_i, outcome = outcome_i)]
    all_results[[length(all_results) + 1]] <- fit
  }

  for (i in seq_len(nrow(single_outcomes))) {

    outcome_i <- single_outcomes$outcome[i]
    domain_i <- single_outcomes$domain[i]
    age_i <- single_outcomes$age_var[i]

    fit <- safe_fit(
      fit_single(outcome_i, genes, score_group, TRUE, age_i),
      paste(score_group, domain_i, outcome_i, sep = " | ")
    )

    if (is.null(fit)) next

    fit[, `:=`(score_group = score_group, domain = domain_i, outcome = outcome_i)]
    all_results[[length(all_results) + 1]] <- fit
  }
}

# ============================================================
# Results + FDR
# ============================================================

if (!length(all_results)) stop("No models ran successfully.")

results <- rbindlist(all_results, fill = TRUE)
h2_details <- rbindlist(all_h2, fill = TRUE)
covariate_results <- rbindlist(all_covariates, fill = TRUE)
failed_models <- rbindlist(failures, fill = TRUE)

results[, p_fdr := p.adjust(p_value, method = "BH"), by = .(score_group, hypothesis)]

if (nrow(h2_details)) {
  h2_fdr <- results[
    hypothesis == "H2",
    .(predictor, score_group, domain, outcome, joint_p_fdr = p_fdr)
  ]

  h2_details <- merge(
    h2_details,
    h2_fdr,
    by = c("predictor", "score_group", "domain", "outcome"),
    all.x = TRUE                                                               # EDIT (lines after 476 not visible)
  )
}

out_dir <- Sys.getenv("ORIGINAL_RESULTS_DIR")                                  # EDIT (output step not visible)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)                    # EDIT
fwrite(results, file.path(out_dir, "results.csv"))                             # EDIT
fwrite(h2_details, file.path(out_dir, "h2_details.csv"))                       # EDIT
fwrite(covariate_results, file.path(out_dir, "covariate_results.csv"))         # EDIT
fwrite(failed_models, file.path(out_dir, "failed_models.csv"))                 # EDIT
