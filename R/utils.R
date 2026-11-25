# =============================================================================
# Shared helpers used by every pipeline step
# =============================================================================

suppressPackageStartupMessages(library(data.table))


# ---- Configuration ----------------------------------------------------------

# Every step is run from the pipeline folder as
#   Rscript R/<step>.R <config file>
# This reads the config file and returns its `cfg` list.
load_config <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1) {
    stop("Usage: Rscript R/<step>.R <config file>\n",
         "  e.g. Rscript R/05_run_regressions.R config/project.R", call. = FALSE)
  }
  require_files(args)
  env <- new.env()
  sys.source(args, envir = env)
  if (!exists("cfg", envir = env)) stop("Config file ", args, " does not define `cfg`.", call. = FALSE)
  env$cfg
}


# ---- Checks that stop the pipeline with an actionable message ---------------

require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing input file(s):\n  ", paste(missing, collapse = "\n  "),
         "\nCheck the paths in your config file.", call. = FALSE)
  }
}

# Stops if `data` lacks any of `columns`. The message lists the missing columns
# and, for each, an R command that searches for similarly named columns.
# `load_command` is R code that loads the data as `d`, e.g. 'd <- readRDS("x.rds")'.
require_columns <- function(data, columns, data_name, load_command) {
  missing <- setdiff(columns, names(data))
  if (!length(missing)) return(invisible(TRUE))

  # Search stem: the name without its last "_suffix" (usually the wave)
  stems <- unique(ifelse(grepl("_", missing), sub("_[^_]*$", "", missing), missing))
  search_commands <- sprintf('  grep("%s", names(d), value = TRUE)', stems)

  stop(sprintf(paste0(
    "%s is missing %d required column(s):\n  %s\n\n",
    "To look for similarly named columns, run in R:\n  %s\n%s\n\n",
    "If a column was renamed, rename it (see `genr_renames` in the config).\n",
    "If it does not exist by design (e.g. a measure not collected at a wave),\n",
    "remove it from the config file."),
    data_name, length(missing), paste(missing, collapse = ", "),
    load_command, paste(search_commands, collapse = "\n")), call. = FALSE)
}

stop_if_duplicated_ids <- function(data, data_name) {
  n_dup <- sum(duplicated(data$idc))
  if (n_dup) stop(data_name, " has ", n_dup, " duplicated idc value(s); expected one row per individual.", call. = FALSE)
}


# ---- Reading and joining ----------------------------------------------------

read_table <- function(path, ...) {
  require_files(path)
  fread(path, ...)
}

# Adds the columns of `extra` to `base` (in place), matching rows on idc and
# keeping the row order of `base` (a left join). Columns present in both tables
# (e.g. MRI age, stored with both MRI and DTI data) must hold identical values.
add_columns_by_idc <- function(base, extra, extra_name) {
  stop_if_duplicated_ids(extra, extra_name)
  rows <- match(base$idc, extra$idc)

  shared <- setdiff(intersect(names(base), names(extra)), "idc")
  for (col in shared) {
    if (!isTRUE(all.equal(base[[col]], extra[[col]][rows], check.attributes = FALSE))) {
      stop("Column ", col, " differs between ", extra_name, " and earlier tables.", call. = FALSE)
    }
  }

  new_cols <- setdiff(names(extra), names(base))
  base[, (new_cols) := extra[rows, new_cols, with = FALSE]]
  invisible(base)
}


# ---- Standardisation --------------------------------------------------------

# z-score over all values
z <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / s
}

# z-score of an individual-level variable in long data (several rows per individual):
# mean and SD are computed over individuals (first row per idc), so individuals with
# more observations do not get more weight.
z_by_id <- function(x, id) {
  first_row <- !duplicated(id)
  s <- sd(x[first_row])
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x[first_row])) / s
}


# ---- GenR export layout (shared by steps 1 and 5) ----------------------------

# Returns the GenR data files and their columns as a named list:
#   covariates = idc + covariates; one entry per outcome domain (lower case).
genr_tables <- function(cfg) {
  tables <- list(covariates = c("idc", cfg$covariates))

  long <- cfg$longitudinal
  for (dom in unique(long$domain)) {
    m <- long[domain == dom]
    tables[[tolower(dom)]] <- unique(c("idc",
                                       paste0(m$age_prefix, "_", m$wave),
                                       paste0(m$column, "_", m$wave)))
  }

  single <- cfg$single_outcomes
  for (dom in unique(single$domain)) {
    m <- single[domain == dom]
    name <- tolower(dom)
    tables[[name]] <- unique(c("idc", tables[[name]], m$outcome, na.omit(m$age_var)))
  }
  tables
}

genr_file <- function(cfg, table) file.path(cfg$genr_csv_dir, paste0(table, ".csv.gz"))


# ---- Analysis data (shared by steps 5 and 7) ----------------------------------

# GenR variables from the step-1 export plus the genetic scores, one row per
# individual. Categorical covariates keep the level order of the original RDS
# file, so the reference category is unchanged (first level = reference).
load_genr_with_scores <- function(cfg) {
  tables <- genr_tables(cfg)
  genr <- read_table(genr_file(cfg, "covariates"))
  stop_if_duplicated_ids(genr, "covariates.csv.gz")
  for (table in setdiff(names(tables), "covariates")) {
    add_columns_by_idc(genr, read_table(genr_file(cfg, table)), paste0(table, ".csv.gz"))
  }
  require_columns(genr, unique(unlist(tables)), "The exported GenR data",
                  load_command = sprintf('d <- data.table::fread("%s")', genr_file(cfg, "<table>")))

  factor_levels <- read_table(genr_file(cfg, "factor_levels"))
  for (v in cfg$categorical_covariates) {
    lv <- factor_levels[variable == v][order(order), level]
    genr[, (v) := if (length(lv)) droplevels(factor(get(v), levels = lv)) else factor(get(v))]
  }

  # Genetic scores; individuals without genotype data get NA
  for (path in cfg$score_files) add_columns_by_idc(genr, read_table(path), path)
  require_columns(genr, unlist(cfg$score_sets, use.names = FALSE), "The genetic score file(s)",
                  load_command = sprintf('d <- data.table::fread("%s")', cfg$score_files[1]))
  genr
}

# Repeated measures in long format: one row per individual x wave x outcome,
# keeping observations with both an outcome value and an age.
# Columns: idc, covariates, age, outcome, domain, outcome_name, wave.
make_long_data <- function(cfg, genr) {
  id_and_covariates <- c("idc", cfg$covariates)
  rbindlist(lapply(seq_len(nrow(cfg$longitudinal)), function(k) {
    m <- cfg$longitudinal[k]
    age_col <- paste0(m$age_prefix, "_", m$wave)
    outcome_col <- paste0(m$column, "_", m$wave)

    d <- genr[, c(id_and_covariates, age_col, outcome_col), with = FALSE]
    setnames(d, c(age_col, outcome_col), c("age", "outcome"))
    d[, `:=`(age = as.numeric(age), outcome = as.numeric(outcome),
             domain = m$domain, outcome_name = m$outcome, wave = m$wave)]
    d[!is.na(outcome) & !is.na(age)]
  }))
}
