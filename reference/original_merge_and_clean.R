# Original "Merging components" and "Cleaning" code, transcribed from
# Code_Hailey_Components.docx. Used only by toy/compare_with_original.R.
# Only edits: paths (marked "# EDIT") so it runs on the toy data, and the
# number of components (toy has 4).

library(data.table)

n_comp <- as.integer(Sys.getenv("TOY_N_COMPONENTS", "20"))            # EDIT
toy_dir <- normalizePath(Sys.getenv("TOY_DIR"))                         # EDIT
out_merge <- file.path(toy_dir, "reference_out", "merged")              # EDIT
dir.create(out_merge, recursive = TRUE, showWarnings = FALSE)           # EDIT

# ---- Merging components ----
setwd(file.path(toy_dir, "ica"))                                         # EDIT (was /gpfs/home4/hdavis/genomica)

snp_info <- fread("variant_catalogue.txt")

snp_info_filtered <- snp_info[V7 == 1]

snp_info_filtered<-snp_info_filtered[,-7]

colnames(snp_info_filtered)<-c("CHR","BP","A1","A2","SNP","CHR:BP:A")

for(i in 1:n_comp) {                                                     # EDIT (was 1:20)

loadings <- fread(paste0("melodic_oIc_", i,".txt"))

pvals <- fread(paste0("p_val_component", i, ".txt"))

snp_info_filtered[, BETA := loadings[[1]]]

snp_info_filtered[, P := pvals[[1]]]

fwrite(snp_info_filtered,paste0(out_merge, "/","component",i,".txt"))   # EDIT (was /scratch-shared/hdavis/)

}

# ---- Cleaning ----
input_dir <- out_merge                                                   # EDIT (was .../components, *_corrected.txt)

output_dir <- file.path(toy_dir, "reference_out", "clean")               # EDIT
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)           # EDIT

for (i in 1:n_comp) {                                                    # EDIT (was 1:20)

  df <- fread(file.path(input_dir,paste0("component",i,".txt")))         # EDIT (was "_corrected.txt")

  n_before <- nrow(df)

  setnames(df,"SNP","RSID")

  df[,SNP := paste0(CHR,":",BP)]

  df<- df[!is.na(CHR) &
          !is.na(BP) &
          is.finite(BETA) &
          is.finite(P) &
            P > 0 & P <= 1]

  df <- df[CHR %in% as.character(1:22)]

  df <- df[
    A1 %in% c("A","C","G","T") &
    A2 %in% c("A","C","G","T") &
    A1 != A2
  ]

  df <- df[!paste0(A1,A2) %in% c("AT","TA","CG","GC")]

  df <- df[!duplicated(SNP) & !duplicated(SNP, fromLast = TRUE)]

  df <- df[, .(SNP, RSID, CHR, BP, A1, A2, BETA, P)]

  fwrite(df, file.path(output_dir, paste0("component_",i,"_QC.txt")),sep = "\t")

  cat("Component",i,":",
      n_before, "->", nrow(df),"SNPs\n")
}
