# =============================================================================
# Step 4 - Subject-level genetic scores from GenR genotypes (bigsnpr)
# =============================================================================
# For every row of cfg$score_weights, the score of individual n is
#
#     score_n = sum over SNPs k of  dosage_nk * w_k
#
# where dosage_nk is the number of effect alleles individual n carries at SNP k and
# w_k is the weight from the summary statistics (BETA, here a Z-score). The
# same computation gives ICA component scores and polygenic scores.
#
# Alleles are aligned with bigsnpr::snp_match() before scoring. It
#   - joins summary statistics and genotypes on chromosome + position,
#   - flips the sign of w_k when the effect allele is the genotype's other allele,
#   - matches alleles reported on the opposite DNA strand (strand_flip = TRUE),
#     dropping strand-ambiguous A/T and C/G SNPs,
#   - drops SNPs that are duplicated after matching.
# Matching is first done with strand_flip = TRUE. If fewer than
# cfg$min_strand_flips matched SNPs needed a strand flip, it is redone with
# strand_flip = FALSE so that strand-ambiguous SNPs are kept (those few
# strand-flipped SNPs are then dropped). The choice is logged per score.
# Positions must be on the same genome build in both files.
#
# Input:  a bigSNP .rds (bigsnpr::snp_readBed()); QC'd summary statistics
#         (step 3 format, A2 = effect allele).
# Output: cfg$scores_file (idc + one column per score) and, next to it,
#         <name>_snp_matching.csv with the number of SNPs used per score.

suppressPackageStartupMessages(library(bigsnpr))
source("R/utils.R")
cfg <- load_config()

# ---- Genotypes ------------------------------------------------------------
require_files(cfg$genotype_rds)
genotypes <- snp_attach(cfg$genotype_rds)
G <- genotypes$genotypes

# Column names expected by snp_match(); in bigsnpr, allele1 is the counted
# allele (a1) and allele2 the other allele (a0).
genotype_map <- setNames(genotypes$map[c("chromosome", "marker.ID", "physical.pos", "allele1", "allele2")],
                         c("chr", "rsid", "pos", "a1", "a0"))
require_columns(genotypes$fam, cfg$genotype_id_column, "The genotype .fam table",
                load_command = sprintf('d <- bigsnpr::snp_attach("%s")$fam', cfg$genotype_rds))
idc <- genotypes$fam[[cfg$genotype_id_column]]
message(sprintf("Genotypes: %d individuals x %d SNPs", nrow(G), ncol(G)))

# Missing genotype calls are replaced by twice the allele frequency (mean
# dosage), so a missing call adds the average contribution to the score.
n_missing <- big_apply(G, function(X, ind) sum(is.na(X[, ind])), a.combine = "+", ncores = cfg$ncores)
if (n_missing > 0) {
  message(n_missing, " missing genotype calls: imputed with the mean dosage")
  G <- snp_fastImputeSimple(G, method = "mean2", ncores = cfg$ncores)
}

# ---- One score per row of cfg$score_weights ---------------------------------
require_files(cfg$score_weights$sumstats)
scores <- data.table(idc = idc)
matching <- list()

for (k in seq_len(nrow(cfg$score_weights))) {
  score_name <- cfg$score_weights$score[k]
  path <- cfg$score_weights$sumstats[k]

  stats <- read_table(path)
  require_columns(stats, c("CHR", "BP", "A1", "A2", "BETA"), path,
                  load_command = sprintf('d <- data.table::fread("%s")', path))
  weights <- data.frame(chr = stats$CHR, pos = stats$BP,
                        a0 = stats$A1, a1 = stats$A2,   # A2 is the effect allele
                        beta = stats$BETA)

  # First match allowing strand flips. If fewer than cfg$min_strand_flips SNPs
  # needed one, match again without strand flipping (keeps ambiguous SNPs).
  match_alleles <- function(strand_flip) {
    snp_match(weights, genotype_map, strand_flip = strand_flip,
              join_by_pos = TRUE, remove_dups = TRUE, return_flip_and_rev = TRUE)
  }
  matched <- match_alleles(strand_flip = TRUE)
  strand_flip_used <- sum(matched$`_FLIP_`) >= cfg$min_strand_flips
  if (!strand_flip_used) matched <- match_alleles(strand_flip = FALSE)

  scores[, (score_name) := big_prodVec(G, matched$beta, ind.col = matched$`_NUM_ID_`,
                                       ncores = cfg$ncores)]
  matching[[k]] <- data.table(score = score_name, sumstats = path,
                              n_sumstats = nrow(stats), n_used = nrow(matched),
                              strand_flip = strand_flip_used,
                              n_sign_flipped = sum(matched$`_REV_`),
                              n_strand_flipped = sum(matched$`_FLIP_`))
  message(sprintf("%-10s %d of %d SNPs used (strand flip %s)", score_name, nrow(matched), nrow(stats),
                  if (strand_flip_used) "on" else "off"))
}

# ---- Write -------------------------------------------------------------------
dir.create(dirname(cfg$scores_file), recursive = TRUE, showWarnings = FALSE)
fwrite(scores, cfg$scores_file)
fwrite(rbindlist(matching), sub("(\\.csv)?(\\.gz)?$", "_snp_matching.csv", cfg$scores_file))
message("Wrote ", ncol(scores) - 1, " scores for ", nrow(scores), " individuals to ", cfg$scores_file)
