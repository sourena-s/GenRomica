# =============================================================================
# Checks the new pipeline against the original code on the toy data
# =============================================================================
# 1. Runs the transcribed original merge + cleaning code and compares its
#    output files with those of steps 2-3.
# 2. Recomputes score C1 by hand (plain R allele matching) and compares it with
#    the bigsnpr score of step 4.
# 3. Builds the original single RDS (GenR data + step-4 scores), runs the
#    transcribed original regression code and compares all result tables with
#    step 5.
# Usage (from the pipeline folder, after toy/run_toy.sh has run steps 1-5):
#   Rscript toy/compare_with_original.R

suppressPackageStartupMessages({ library(data.table); library(bigsnpr) })
toy <- "toy/data"
ok <- TRUE
check <- function(label, result) {
  passed <- isTRUE(result)
  cat(sprintf("[%s] %s%s\n", if (passed) "PASS" else "FAIL", label, if (passed) "" else paste0(": ", paste(result, collapse = "; "))))
  if (!passed) ok <<- FALSE
}

# ---- 1. Merge and cleaning ---------------------------------------------------
Sys.setenv(TOY_DIR = toy, TOY_N_COMPONENTS = "4")
system2("Rscript", "reference/original_merge_and_clean.R", stdout = FALSE)

for (i in 1:4) {
  check(sprintf("merged component %d identical", i),
        all.equal(fread(sprintf("%s/reference_out/merged/component%d.txt", toy, i)),
                  fread(sprintf("%s/pipeline_out/merged/component%d.txt", toy, i))))
  check(sprintf("cleaned component %d identical", i),
        all.equal(fread(sprintf("%s/reference_out/clean/component_%d_QC.txt", toy, i)),
                  fread(sprintf("%s/pipeline_out/clean/component_%d_QC.txt", toy, i))))
}

# ---- 2. Score C1 by hand -----------------------------------------------------
geno <- snp_attach(file.path(toy, "genotypes", "toy.rds"))
G <- geno$genotypes[]
# mean imputation, rounded to 2 decimals as in snp_fastImputeSimple(method = "mean2")
G <- apply(G, 2, function(g) { g[is.na(g)] <- round(mean(g, na.rm = TRUE), 2); g })
map <- as.data.table(geno$map)[, k := .I]

stats <- fread(file.path(toy, "pipeline_out", "clean", "component_1_QC.txt"))
complement <- c(A = "T", T = "A", C = "G", G = "C")
m <- merge(stats, map, by.x = c("CHR", "BP"), by.y = c("chromosome", "physical.pos"))
m[, weight := fcase(
  A2 == allele1 & A1 == allele2, BETA,                                  # effect allele counted
  A2 == allele2 & A1 == allele1, -BETA,                                 # alleles swapped
  complement[A2] == allele1 & complement[A1] == allele2, BETA,          # other strand
  complement[A2] == allele2 & complement[A1] == allele1, -BETA,
  default = NA_real_)]
m <- m[!is.na(weight)]
manual <- data.table(idc = geno$fam$sample.ID, C1 = as.vector(G[, m$k] %*% m$weight))
scores <- fread(file.path(toy, "pipeline_out", "scores", "genr_scores.csv.gz"))
check("score C1 equals hand computation", all.equal(manual$C1, scores$C1[match(manual$idc, scores$idc)]))

# ---- 3. Regressions ------------------------------------------------------------
original_dir <- file.path(toy, "reference_out", "rds")
dir.create(original_dir, recursive = TRUE, showWarnings = FALSE)
genr <- as.data.table(readRDS(file.path(toy, "genr", "genr_full_file.rds")))
genr[, c("C1", "C2") := NULL]                                    # stale columns: use the step-4 scores
genr <- merge(genr, scores, by = "idc", all.x = TRUE, sort = FALSE)
saveRDS(as.data.frame(genr), file.path(original_dir, "genr_full_file.rds"))

Sys.setenv(ORIGINAL_DATA_DIR = original_dir,
           ORIGINAL_RESULTS_DIR = file.path(toy, "reference_out", "results"),
           TOY_N_PCS = "4", TOY_N_UKPCS = "3", TOY_N_PGS = "2")
message("Running the original regression code on the toy data ...")
system2("Rscript", "reference/original_regression.R", stdout = FALSE, stderr = FALSE)

compare_table <- function(name, keys) {
  a <- fread(file.path(toy, "reference_out", "results", name))
  b <- fread(file.path(toy, "pipeline_out", "results", name))
  if (!nrow(a) && !nrow(b)) return(TRUE)
  if (nrow(a) != nrow(b)) return(sprintf("%d vs %d rows", nrow(a), nrow(b)))
  setcolorder(b, names(a))
  setkeyv(a, keys); setkeyv(b, keys)
  all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE)
}
check("results.csv matches original (tolerance 1e-10)",
      compare_table("results.csv", c("score_group", "predictor", "domain", "outcome", "hypothesis")))
check("h2_details.csv matches original",
      compare_table("h2_details.csv", c("score_group", "predictor", "domain", "outcome")))
check("covariate_results.csv matches original",
      compare_table("covariate_results.csv", c("score_group", "predictor", "domain", "outcome", "term")))
check("failed_models.csv matches original", compare_table("failed_models.csv", "model"))

cat(if (ok) "\nAll checks passed.\n" else "\nSome checks FAILED.\n")
if (!ok) quit(status = 1)
