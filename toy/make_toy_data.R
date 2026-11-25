# =============================================================================
# Toy data for testing the pipeline on a laptop
# =============================================================================
# Creates, in toy/data/:
#   genr/genr_full_file.rds  simulated GenR data (500 children), including the
#                            problems the pipeline must handle: sum_wid_f05
#                            naming, missing values, stale score columns
#   genotypes/toy.{bed,bim,fam,rds,bk}  480 genotyped children (some not in
#                            GenR, some GenR children not genotyped), 3000 SNPs
#   ica/                     variant catalogue + per-component Z and p-value files,
#                            with SNPs that should be removed by QC or flipped by
#                            snp_match (swapped alleles, opposite strand)
#   sumstats_other/          QC'd summary statistics for toy UK_PCS and PGS scores
#
# Usage (from the pipeline folder): Rscript toy/make_toy_data.R

suppressPackageStartupMessages({ library(data.table); library(bigsnpr) })
set.seed(2026)

out <- "toy/data"
unlink(out, recursive = TRUE)
for (d in c("genr", "genotypes", "ica", "sumstats_other")) dir.create(file.path(out, d), recursive = TRUE)

n_components <- 4
complement <- c(A = "T", T = "A", C = "G", G = "C")


# ---- Genotypes ---------------------------------------------------------------
n_geno <- 480
m_snps <- 3000

genr_ids <- sample(100001:100500)                               # GenR idc, shuffled order
geno_ids <- c(sample(genr_ids, 460), 900001:900020)             # 40 GenR children not genotyped, 20 extra

non_palindromic <- list(c("A","C"), c("A","G"), c("T","C"), c("T","G"), c("C","A"), c("G","A"), c("C","T"), c("G","T"))
alleles <- do.call(rbind, sample(non_palindromic, m_snps, replace = TRUE))

map <- data.table(
  chromosome   = sort(sample(1:22, m_snps, replace = TRUE)),
  marker.ID    = paste0("rs", 1:m_snps),
  genetic.dist = 0,
  physical.pos = NA_integer_,
  allele1      = alleles[, 1],   # counted allele
  allele2      = alleles[, 2]
)
map[, physical.pos := sort(sample(1e4:1e8, .N)), by = chromosome]

maf <- runif(m_snps, 0.05, 0.5)
geno <- sapply(maf, function(p) rbinom(n_geno, 2, p))
geno[sample(length(geno), 0.01 * length(geno))] <- NA           # 1% missing calls

fake <- snp_fake(n_geno, m_snps)
geno_codes <- geno
geno_codes[is.na(geno_codes)] <- 3L                            # bigsnpr code for a missing call
fake$genotypes[] <- as.raw(geno_codes)
fake$map <- as.data.frame(map)
fake$fam$sample.ID <- geno_ids
fake$fam$family.ID <- geno_ids
bed <- snp_writeBed(fake, file.path(out, "genotypes", "toy.bed"))
invisible(snp_readBed(bed, backingfile = file.path(out, "genotypes", "toy")))


# ---- Variant catalogue (ICA study) ------------------------------------------
# A1 = other allele, A2 = effect allele. 2500 genotyped SNPs, reported as
#   - same alleles: A2 = genotype allele1
#   - swapped:      A2 = genotype allele2 (snp_match flips the sign)
#   - other strand: complement of allele1/allele2
# plus SNPs QC must remove and SNPs absent from the genotypes.
idx <- sort(sample(m_snps, 2500))
orientation <- sample(c("same", "swapped", "strand"), length(idx), replace = TRUE, prob = c(.6, .3, .1))
a1 <- ifelse(orientation == "same", map$allele2[idx], map$allele1[idx])
a2 <- ifelse(orientation == "same", map$allele1[idx], map$allele2[idx])
a1[orientation == "strand"] <- complement[map$allele2[idx][orientation == "strand"]]
a2[orientation == "strand"] <- complement[map$allele1[idx][orientation == "strand"]]

catalogue <- data.table(CHR = map$chromosome[idx], BP = map$physical.pos[idx], A1 = a1, A2 = a2,
                        SNP = map$marker.ID[idx])
opposite_strand_snps <- catalogue$SNP[orientation == "strand"]

problems <- data.table(
  CHR = c(23L, 23L, 5L, 6L, 7L, 8L, 9L, rep(10L, 2), 11L, sample(1:22, 300, replace = TRUE)),
  BP  = as.integer(c(5e5, 6e5, 1.5e8, 1.5e8, 1.5e8, 1.5e8, 1.5e8, rep(1.6e8, 2), 1.5e8, sample(1.2e8:1.3e8, 300))),
  A1  = c("A", "C", "A", "C", "a", "A", "A", "A", "G", "G", sample(c("A", "G"), 300, replace = TRUE)),
  A2  = c("G", "T", "T", "G", "g", "A", "N", "G", "T", "C", sample(c("C", "T"), 300, replace = TRUE)),
  SNP = paste0("rsP", 1:310)
)
# rows: chrX x2, palindromic A/T and C/G, lowercase, identical alleles, non-ACGT,
# duplicated position x2, palindromic G/C, 300 SNPs not genotyped
catalogue <- rbind(catalogue, problems)
catalogue[, `CHR:BP:A` := paste0(CHR, ":", BP, ":", A2)]
catalogue[, include := rbinom(.N, 1, 0.95)]                    # 5% not used in the ICA
fwrite(catalogue, file.path(out, "ica", "variant_catalogue.txt"), col.names = FALSE, sep = " ")

n_included <- sum(catalogue$include == 1)
for (i in seq_len(n_components)) {
  zs <- rnorm(n_included)
  zs[sample(n_included, 3)] <- NA                               # removed by QC
  ps <- 2 * pnorm(-abs(zs))
  ps[sample(which(!is.na(zs)), 2)] <- 0                         # removed by QC (P > 0)
  fwrite(data.table(zs), file.path(out, "ica", sprintf("melodic_oIc_%d.txt", i)), col.names = FALSE, na = "NA")
  fwrite(data.table(ps), file.path(out, "ica", sprintf("p_val_component%d.txt", i)), col.names = FALSE, na = "NA")
}


# ---- Other summary statistics (already QC'd, step-3 format) -------------------
clean <- catalogue[include == 1 & CHR <= 22 & A1 %in% c("A","C","G","T") & A2 %in% c("A","C","G","T") & A1 != A2 &
                   !paste0(A1, A2) %in% c("AT","TA","CG","GC")]
clean <- clean[!duplicated(paste(CHR, BP)) & !duplicated(paste(CHR, BP), fromLast = TRUE)]
# The PGS weights contain no opposite-strand SNPs, so step 4 matches them
# without strand flipping (tests both branches of that step).
for (name in c(paste0("UK_PCS_", 1:3), paste0("PGS_", 1:2))) {
  pool <- if (startsWith(name, "PGS")) clean[!SNP %in% opposite_strand_snps] else clean
  s <- pool[sample(.N, 1500), .(SNP = paste0(CHR, ":", BP), RSID = SNP, CHR, BP, A1, A2,
                                 BETA = rnorm(1500), P = runif(1500))]
  fwrite(s, file.path(out, "sumstats_other", paste0(name, ".txt")), sep = "\t")
}


# ---- GenR data ---------------------------------------------------------------
n <- length(genr_ids)
waves <- c(f05 = 5, f09 = 9.5, f13 = 13.5, f17 = 17)
with_missing <- function(x, p = 0.15) { x[runif(length(x)) < p] <- NA; x }

genr <- data.table(
  idc = genr_ids,
  sex = factor(sample(c("girl", "boy"), n, replace = TRUE), levels = c("girl", "boy")),  # non-alphabetical reference
  genotype_cohort = sample(1:2, n, replace = TRUE, prob = c(.7, .3)),
  gPC1 = rnorm(n, 0, 0.01), gPC2 = rnorm(n, 0, 0.01), gPC3 = rnorm(n, 0, 0.01), gPC4 = rnorm(n, 0, 0.01)
)

child_level <- rnorm(n)                                          # stable child effect
for (w in names(waves)) {
  genr[, (paste0("age_cbcl_", w))      := with_missing(waves[[w]] + rnorm(n, 0, 0.4), 0.05)]
  genr[, (paste0("age_child_mri_", w)) := with_missing(waves[[w]] + 0.8 + rnorm(n, 0, 0.4), 0.05)]
  for (s in c("anx", "som", "wit", "sop", "tho", "att", "rul", "agg")) {
    name <- if (s == "wit" && w == "f05") "sum_wid_f05" else paste0("sum_", s, "_", w)
    genr[, (name) := with_missing(rpois(n, exp(1 + 0.3 * child_level)))]
  }
  genr[, (paste0("TotalGrayVol_", w))       := with_missing(7e5 + 2e4 * waves[[w]] + 3e4 * child_level + rnorm(n, 0, 2e4), 0.4)]
  genr[, (paste0("cortical_thickness_", w)) := with_missing(3 - 0.03 * waves[[w]] + 0.1 * child_level + rnorm(n, 0, 0.1), 0.4)]
  genr[, (paste0("mean_FA_genr_", w))       := with_missing(0.4 + 0.005 * waves[[w]] + rnorm(n, 0, 0.02), 0.45)]
  genr[, (paste0("mean_MD_genr_", w))       := with_missing(8e-4 - 1e-5 * waves[[w]] + rnorm(n, 0, 3e-5), 0.45)]
}
genr[, `:=`(
  INCOME5 = with_missing(sample(1:5, n, replace = TRUE)),
  EDUCM   = with_missing(sample(1:3, n, replace = TRUE)),
  EDUCP   = with_missing(sample(1:3, n, replace = TRUE), 0.3),
  DietScore_8y = with_missing(rnorm(n, 5, 1.5)),
  onlineweek_child_F17 = with_missing(rgamma(n, 2, 0.2), 0.4),
  AgeChild_GR1203 = with_missing(rnorm(n, 17.5, 0.5), 0.05),
  # Columns the pipeline must ignore
  C1 = rnorm(n), C2 = rnorm(n), unrelated_variable = rnorm(n)
)]
saveRDS(as.data.frame(genr), file.path(out, "genr", "genr_full_file.rds"))

message("Toy data written to ", out)
