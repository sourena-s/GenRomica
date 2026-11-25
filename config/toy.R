# =============================================================================
# Toy configuration: the project settings, pointed at toy/data and scaled down
# =============================================================================
# Create the toy data first:  Rscript toy/make_toy_data.R

source("config/project.R", local = TRUE)

toy <- "toy/data"

cfg$genr_rds     <- file.path(toy, "genr", "genr_full_file.rds")
cfg$genr_csv_dir <- file.path(toy, "genr_csv")

cfg$n_components      <- 4
cfg$variant_catalogue <- file.path(toy, "ica", "variant_catalogue.txt")
cfg$zscore_pattern    <- file.path(toy, "ica", "melodic_oIc_%d.txt")
cfg$pvalue_pattern    <- file.path(toy, "ica", "p_val_component%d.txt")
cfg$merged_pattern    <- file.path(toy, "pipeline_out", "merged", "component%d.txt")

cfg$qc_input_pattern  <- cfg$merged_pattern   # no "corrected" step in the toy
cfg$qc_output_pattern <- file.path(toy, "pipeline_out", "clean", "component_%d_QC.txt")

cfg$genotype_rds <- file.path(toy, "genotypes", "toy.rds")
cfg$ncores       <- 1
cfg$score_weights <- data.table(
  score    = c(paste0("C", 1:4), paste0("UK_PCS_", 1:3), paste0("PGS_", 1:2)),
  sumstats = c(sprintf(cfg$qc_output_pattern, 1:4),
               file.path(toy, "sumstats_other", paste0(c(paste0("UK_PCS_", 1:3), paste0("PGS_", 1:2)), ".txt")))
)
cfg$scores_file <- file.path(toy, "pipeline_out", "scores", "genr_scores.csv.gz")
cfg$score_files <- cfg$scores_file

cfg$results_dir <- file.path(toy, "pipeline_out", "results")
cfg$score_sets <- list(
  PCS    = paste0("C", 1:4),
  UK_PCS = paste0("UK_PCS_", 1:3),
  PGS    = paste0("PGS_", 1:2)
)
cfg$combined_sets <- list(
  PCS_combined    = cfg$score_sets$PCS,
  UK_PCS_combined = cfg$score_sets$UK_PCS
)
