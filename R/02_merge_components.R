# =============================================================================
# Step 2 - Attach ICA component Z-scores and p-values to the variant catalogue
# =============================================================================
# Inputs (from the ICA study):
#   variant_catalogue.txt    no header, 7 columns:
#                              CHR, BP, A1 (other allele), A2 (effect allele),
#                              SNP (rsID), CHR:BP:A, include flag (1 = used in ICA)
#   melodic_oIc_<i>.txt      one column, no header: Z-score of each SNP on component i
#   p_val_component<i>.txt   one column, no header: p-value of each SNP on component i
#
# The Z-score and p-value files have no SNP identifier: row k belongs to the
# k-th included SNP of the catalogue. This script checks that the row counts
# match before combining them.
#
# Output: component<i>.txt (comma-separated) with columns
#   CHR, BP, A1, A2, SNP, CHR:BP:A, BETA, P
# BETA holds the Z-score, whose sign refers to the effect allele A2. It is
# called BETA by convention so standard summary-statistics tools recognise it.

source("R/utils.R")
cfg <- load_config()

# Reads a one-column, header-less numeric file.
read_one_numeric_column <- function(path, what) {
  x <- read_table(path, header = FALSE)
  if (ncol(x) != 1) stop(path, ": expected 1 column of ", what, ", found ", ncol(x), ".", call. = FALSE)
  if (!is.numeric(x[[1]])) stop(path, ": ", what, " are not numeric (does the file have a header line?).", call. = FALSE)
  x[[1]]
}

# ---- Variant catalogue -----------------------------------------------------
catalogue <- read_table(cfg$variant_catalogue, header = FALSE)
if (ncol(catalogue) != 7) stop("Variant catalogue: expected 7 columns, found ", ncol(catalogue), ".", call. = FALSE)
setnames(catalogue, c("CHR", "BP", "A1", "A2", "SNP", "CHR:BP:A", "include"))

catalogue <- catalogue[include == 1]
catalogue[, include := NULL]
message(nrow(catalogue), " SNPs flagged for inclusion in the variant catalogue")

# ---- One output file per component -------------------------------------------
dir.create(dirname(sprintf(cfg$merged_pattern, 1)), recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(cfg$n_components)) {
  zscores <- read_one_numeric_column(sprintf(cfg$zscore_pattern, i), "Z-scores")
  pvalues <- read_one_numeric_column(sprintf(cfg$pvalue_pattern, i), "p-values")

  # Safeguard: rows are matched by position, so the counts must be identical
  if (length(zscores) != nrow(catalogue) || length(pvalues) != nrow(catalogue)) {
    stop(sprintf("Component %d: %d Z-scores and %d p-values, but %d included SNPs in the catalogue.",
                 i, length(zscores), length(pvalues), nrow(catalogue)), call. = FALSE)
  }

  component <- copy(catalogue)
  component[, `:=`(BETA = zscores, P = pvalues)]
  fwrite(component, sprintf(cfg$merged_pattern, i))
}
message("Wrote ", cfg$n_components, " component files: ", sprintf(cfg$merged_pattern, 1), ", ...")
