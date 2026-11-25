# =============================================================================
# Step 1 - Export the GenR RDS file to one csv.gz file per data type
# =============================================================================
# Run once. Later steps read only these plain-text files (with data.table::fread),
# so the variables used by the analysis are explicit and do not depend on R's
# binary RDS format. Output (in cfg$genr_csv_dir, each with a header row):
#
#   covariates.csv.gz      idc, sex, genotype_cohort, gPC1-gPC4
#   cbcl.csv.gz            idc, age_cbcl_<wave>, sum_<scale>_<wave>
#   mri.csv.gz, dti.csv.gz idc, age_child_mri_<wave>, imaging measures
#   ses.csv.gz, diet.csv.gz, online.csv.gz   single-time-point outcomes
#   factor_levels.csv.gz   level order of factor covariates (keeps the
#                          reference category used in the regressions)
#
# The columns come from config (see genr_tables() in R/utils.R). Genetic scores
# are deliberately NOT exported: they are based on weights from another study
# and are computed from the genotypes in step 4.

source("R/utils.R")
cfg <- load_config()

require_files(cfg$genr_rds)
genr <- as.data.table(readRDS(cfg$genr_rds))
message("Read ", nrow(genr), " individuals x ", ncol(genr), " variables from ", cfg$genr_rds)

# ---- Fix known naming inconsistencies --------------------------------------
for (old in names(cfg$genr_renames)) {
  new <- cfg$genr_renames[[old]]
  if (old %in% names(genr) && !new %in% names(genr)) {
    setnames(genr, old, new)
    message("Renamed ", old, " -> ", new)
  }
}

# ---- Check that every needed variable exists and has the expected type ------
tables <- genr_tables(cfg)
needed <- unique(unlist(tables))

require_columns(genr, needed, "The GenR RDS file",
                load_command = sprintf('d <- readRDS("%s")', cfg$genr_rds))
stop_if_duplicated_ids(genr, "The GenR RDS file")

# All variables except idc and the categorical covariates enter the models as
# numbers. A factor here would be turned into level codes, so stop instead.
numeric_vars <- setdiff(needed, c("idc", cfg$categorical_covariates))
not_numeric  <- numeric_vars[!vapply(genr[, ..numeric_vars], is.numeric, logical(1))]
if (length(not_numeric)) {
  stop("These variables should be numeric but are not:\n  ",
       paste(sprintf("%s (%s)", not_numeric,
                     vapply(genr[, ..not_numeric], function(x) class(x)[1], "")), collapse = "\n  "),
       "\n\nInspect them in R with:\n",
       sprintf('  d <- readRDS("%s")\n', cfg$genr_rds),
       sprintf('  str(d[, c(%s)])', paste0('"', not_numeric, '"', collapse = ", ")),
       "\nand decide on a numeric coding before exporting.", call. = FALSE)
}

# ---- Record factor level order (the first level is the reference category) --
factor_levels <- rbindlist(lapply(cfg$categorical_covariates, function(v) {
  if (!is.factor(genr[[v]])) return(NULL)
  data.table(variable = v, level = levels(genr[[v]]), order = seq_along(levels(genr[[v]])))
}))
if (!nrow(factor_levels)) factor_levels <- data.table(variable = character(), level = character(), order = integer())

# ---- Write -----------------------------------------------------------------
dir.create(cfg$genr_csv_dir, recursive = TRUE, showWarnings = FALSE)

# fwrite() rounds decimals to 15 significant digits. To reproduce results from
# the RDS exactly, decimal numbers are written with 17 digits, which read back
# into exactly the same value.
full_precision <- function(d) {
  for (col in names(d)) {
    x <- d[[col]]
    if (is.double(x) && !is.factor(x)) set(d, j = col, value = ifelse(is.na(x), NA_character_, sprintf("%.17g", x)))
  }
  d
}

for (table in names(tables)) {
  cols <- tables[[table]]
  fwrite(full_precision(genr[, ..cols]), genr_file(cfg, table))
  message(sprintf("Wrote %-12s %d variables", paste0(table, ".csv.gz"), length(cols) - 1))
}
fwrite(factor_levels, genr_file(cfg, "factor_levels"))

# ---- Report what was left out ----------------------------------------------
not_exported <- setdiff(names(genr), needed)
message(length(not_exported), " variable(s) in the RDS are not used and were not exported.")

score_like <- grep("^(C[0-9]+|UK_PCS_[0-9]+|PGS_[0-9]+)$", not_exported, value = TRUE)
if (length(score_like)) {
  message("Note: the RDS contains genetic score columns (", paste(head(score_like, 5), collapse = ", "),
          if (length(score_like) > 5) ", ..." else "",
          "). They were not exported; scores are computed in step 4.")
}
