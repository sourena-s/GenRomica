# =============================================================================
# Step 3 - Quality control of the component summary statistics
# =============================================================================
# The input files are never modified: cleaned copies are written to a separate
# location (cfg$qc_output_pattern).
#
# SNPs are removed, in this order, when they have
#   1. a missing position, a non-finite Z-score (BETA) or a p-value outside (0, 1]
#   2. a chromosome other than 1-22
#   3. an allele other than A/C/G/T, or two identical alleles
#   4. strand-ambiguous (palindromic) alleles: A/T or C/G
#   5. a position (CHR:BP) that occurs more than once; all copies are removed
#
# Output: component_<i>_QC.txt (tab-separated) with columns
#   SNP (CHR:BP), RSID, CHR, BP, A1 (other allele), A2 (effect allele), BETA (Z), P

source("R/utils.R")
cfg <- load_config()

dir.create(dirname(sprintf(cfg$qc_output_pattern, 1)), recursive = TRUE, showWarnings = FALSE)
nucleotides <- c("A", "C", "G", "T")

for (i in seq_len(cfg$n_components)) {
  input <- sprintf(cfg$qc_input_pattern, i)
  stats <- read_table(input)
  require_columns(stats, c("SNP", "CHR", "BP", "A1", "A2", "BETA", "P"), input,
                  load_command = sprintf('d <- data.table::fread("%s")', input))
  if ("RSID" %in% names(stats)) stop(input, " already has an RSID column.", call. = FALSE)

  # The SNP column holds rsIDs; the new SNP identifier is CHR:BP
  setnames(stats, "SNP", "RSID")
  stats[, SNP := paste0(CHR, ":", BP)]

  n <- c(input = nrow(stats))

  stats <- stats[!is.na(CHR) & !is.na(BP) & is.finite(BETA) & is.finite(P) & P > 0 & P <= 1]
  n["valid values"] <- nrow(stats)

  stats <- stats[CHR %in% as.character(1:22)]
  n["autosomes"] <- nrow(stats)

  # TODO: this also removes indels (alleles longer than one base). Check whether
  # indels can be aligned with the GenR genotypes and kept.
  stats <- stats[A1 %in% nucleotides & A2 %in% nucleotides & A1 != A2]
  n["ACGT alleles"] <- nrow(stats)

  stats <- stats[!paste0(A1, A2) %in% c("AT", "TA", "CG", "GC")]
  n["not palindromic"] <- nrow(stats)

  stats <- stats[!duplicated(SNP) & !duplicated(SNP, fromLast = TRUE)]
  n["unique position"] <- nrow(stats)

  stats <- stats[, .(SNP, RSID, CHR, BP, A1, A2, BETA, P)]
  fwrite(stats, sprintf(cfg$qc_output_pattern, i), sep = "\t")

  removed <- -diff(n)
  message(sprintf("Component %2d: %d -> %d SNPs (removed: %s)", i, n[1], n[length(n)],
                  paste(sprintf("%s %d", names(removed), removed), collapse = ", ")))
}
