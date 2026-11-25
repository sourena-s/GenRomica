# =============================================================================
# Project configuration: file paths and analysis settings
# =============================================================================
# Every step takes this file as its only argument, e.g.
#   Rscript R/05_run_regressions.R config/project.R
# Paths are those of the original scripts; adjust them to the machine you use.
# File name patterns use sprintf() syntax: %d is replaced by the component number.

cfg <- list()


# ---- Step 1: export of the GenR RDS file ------------------------------------
cfg$genr_rds     <- "/mnt/data/Hailey/data/genr_full_file.rds"
cfg$genr_csv_dir <- "/mnt/data/Hailey/data/genr_csv"   # output: one csv.gz per data type

# Columns renamed during export (old = new). Applied only when the new name is
# absent. sum_wid_f05 is the CBCL "withdrawn" scale at wave f05, named
# sum_wit_<wave> at the other waves.
cfg$genr_renames <- c(sum_wid_f05 = "sum_wit_f05")


# ---- Step 2: attach ICA component Z-scores to the variant catalogue ---------
cfg$n_components      <- 20
cfg$variant_catalogue <- "/gpfs/home4/hdavis/genomica/variant_catalogue.txt"
cfg$zscore_pattern    <- "/gpfs/home4/hdavis/genomica/melodic_oIc_%d.txt"
cfg$merged_pattern    <- "/scratch-shared/hdavis/component%d.txt"
cfg$qc_plot_dir       <- "/scratch-shared/hdavis/component_qc"   # Z-score histograms
cfg$z_hist_binwidth   <- 0.05                                    # histogram bin width (Z units)


# ---- Step 3: QC of the component summary statistics -------------------------
# NOTE: the original QC script reads "*_corrected.txt" files, produced by a step
# that is not part of this pipeline (see README). Point this at the merged
# files (cfg$merged_pattern) if no correction is needed.
cfg$qc_input_pattern  <- "/mnt/data/Hailey/sumstats/components/component%d_corrected.txt"
cfg$qc_output_pattern <- "/mnt/data/Hailey/sumstats/components_clean/component_%d_QC.txt"


# ---- Step 4: subject-level scores from GenR genotypes (bigsnpr) -------------
cfg$genotype_rds       <- "/mnt/data/Hailey/genotypes/genr.rds"  # made with bigsnpr::snp_readBed()
cfg$genotype_id_column <- "sample.ID"   # column of the .fam file that holds the GenR idc
cfg$ncores             <- 4

# snp_match() is run with strand_flip = TRUE first; if fewer than this many
# matched SNPs needed a strand flip, it is rerun with strand_flip = FALSE to
# keep strand-ambiguous (A/T, C/G) SNPs.
cfg$min_strand_flips <- 10

# Which summary statistics make which score. Every file must have the step-3
# output format (CHR, BP, A1, A2, BETA, ...; A2 = effect allele).
cfg$score_weights <- data.table(
  score    = paste0("C", 1:20),
  sumstats = sprintf(cfg$qc_output_pattern, 1:20)
)
# TODO: add rows for UK_PCS_1..20 and PGS_1..18 once their summary statistics
# are located, or list a file with these precomputed scores in cfg$score_files.
cfg$scores_file <- "/mnt/data/Hailey/scores/genr_scores.csv.gz"


# ---- Step 5: regression analyses ----------------------------------------------
cfg$score_files <- c(cfg$scores_file)   # all files with idc + score columns
cfg$results_dir <- "/mnt/data/Hailey/results"

cfg$genetic_pcs <- paste0("gPC", 1:4)
cfg$covariates  <- c("sex", "genotype_cohort", cfg$genetic_pcs)
cfg$categorical_covariates <- c("sex", "genotype_cohort")

# Each score is tested alone; each combined set is also tested jointly.
cfg$score_sets <- list(
  PCS    = paste0("C", 1:20),
  UK_PCS = paste0("UK_PCS_", 1:20),
  PGS    = paste0("PGS_", 1:18)
)
cfg$combined_sets <- list(
  PCS_combined    = cfg$score_sets$PCS,
  UK_PCS_combined = cfg$score_sets$UK_PCS
)

# Age (years) at which age is centred in the mixed models
cfg$reference_age <- 17

# Repeated measures: one row per (column, wave). The data column is
# <column>_<wave>; the age at that measurement is <age_prefix>_<wave>.
# Every listed column must exist. Measures not collected at a wave must be
# removed from this table, e.g.
#   cfg$longitudinal <- cfg$longitudinal[!(column == "sum_rul" & wave == "f05")]
repeated_measure <- function(domain, columns, age_prefix, waves, strip_prefix = "") {
  CJ(column = columns, wave = waves, sorted = FALSE)[, `:=`(
    domain     = domain,
    outcome    = sub(paste0("^", strip_prefix), "", column),  # label used in the results
    age_prefix = age_prefix
  )][]
}

cfg$waves <- c("f05", "f09", "f13", "f17")
cfg$longitudinal <- rbind(
  repeated_measure("CBCL", paste0("sum_", c("anx", "som", "wit", "sop", "tho", "att", "rul", "agg")),
                   "age_cbcl", cfg$waves, strip_prefix = "sum_"),
  repeated_measure("MRI", c("TotalGrayVol", "cortical_thickness"), "age_child_mri", cfg$waves),
  repeated_measure("DTI", c("mean_FA_genr", "mean_MD_genr"), "age_child_mri", cfg$waves)
)

# Uncomment to remove the CBCL f05 timepoint from all analyses. At f05 the
# CBCL/1½-5 was used, which is not comparable with the CBCL/6-18 used at later
# waves (README, concern D). The CBCL models then use waves f09, f13 and f17 only.
# cfg$longitudinal <- cfg$longitudinal[!(domain == "CBCL" & wave == "f05")]

# Outcomes measured once. age_var (optional) is added as a covariate.
cfg$single_outcomes <- data.table(
  outcome = c("INCOME5", "EDUCM", "EDUCP", "DietScore_8y", "onlineweek_child_F17"),
  domain  = c("SES",     "SES",   "SES",   "Diet",         "Online"),
  age_var = c(NA,        NA,      NA,      NA,             "AgeChild_GR1203")
)
