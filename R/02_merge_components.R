# =============================================================================
# Step 2 - Attach ICA component Z-scores to the variant catalogue
# =============================================================================
# Inputs (from the ICA study):
#   variant_catalogue.txt    no header, 7 columns:
#                              CHR, BP, A1 (other allele), A2 (effect allele),
#                              SNP (rsID), CHR:BP:A, include flag (1 = used in ICA)
#   melodic_oIc_<i>.txt      one column, no header: Z-score of each SNP on component i
#
# The Z-score files have no SNP identifier: row k belongs to the k-th included
# SNP of the catalogue. This script checks that the row counts match before
# combining them. Two-sided p-values are computed from the Z-scores,
# P = 2 * Phi(-|Z|), where Phi is the standard normal distribution function.
#
# Output
#   component<i>.txt (comma-separated) with columns
#     CHR, BP, A1, A2, SNP, CHR:BP:A, BETA, P
#   BETA holds the Z-score, whose sign refers to the effect allele A2. It is
#   called BETA by convention so standard summary-statistics tools recognise it.
#
#   QC, in cfg$qc_plot_dir:
#     component<i>_z_histogram.png  histogram of the Z-scores with the standard
#                                   normal density (dashed), on a linear and a
#                                   log scale (the log scale shows the tails)
#     component_z_summary.csv       per component: n, mean, SD, genomic inflation
#                                   factor lambda, number of |Z| > 5.45 (P < 5e-8)

source("R/utils.R")
cfg <- load_config()

# Reads a one-column, header-less numeric file.
read_one_numeric_column <- function(path, what) {
  x <- read_table(path, header = FALSE)
  if (ncol(x) != 1) stop(path, ": expected 1 column of ", what, ", found ", ncol(x), ".", call. = FALSE)
  if (!is.numeric(x[[1]])) stop(path, ": ", what, " are not numeric (does the file have a header line?).", call. = FALSE)
  x[[1]]
}

# Histogram of Z-scores against the standard normal density.
# Components have millions of SNPs, so a fixed bin width (cfg$z_hist_binwidth,
# default 0.05) is used rather than a rule that grows the number of bins with
# n: with 1 million SNPs a 0.05-wide bin near 0 holds ~20,000 SNPs, so the
# histogram itself is smooth and any departure from N(0, 1) is real. The x-axis
# shows +-(99.99th percentile of |Z|), kept between 5 and 10; SNPs beyond it
# are counted below the plot. Binning is done once (hist() is linear in n), so
# millions of SNPs take about a second. In the log-scale panel each non-empty
# bin is a point; the y-axis floor is the density of a single SNP in a bin.
plot_z_histogram <- function(z, title, path, binwidth) {
  z <- z[is.finite(z)]
  breaks <- seq(floor(min(z) / binwidth) * binwidth,
                ceiling(max(z) / binwidth) * binwidth + binwidth / 2, by = binwidth)
  h <- hist(z, breaks = breaks, plot = FALSE)
  x_max <- min(10, max(5, ceiling(quantile(abs(z), 0.9999))))
  n_outside <- sum(abs(z) > x_max)
  lambda <- median(z^2) / qchisq(0.5, df = 1)

  # Style
  col_data   <- "#3A3F44"   # data outline / points
  col_fill   <- "#C9CED6"   # histogram fill
  col_normal <- "#1F5FBF"   # N(0, 1) reference, blue dashed
  shown <- abs(h$mids) <= x_max

  png(path, width = 10, height = 4.3, units = "in", res = 300, type = "cairo")
  on.exit(dev.off())
  par(mfrow = c(1, 2), mar = c(3.8, 4.6, 2, 1), oma = c(2.4, 0, 2.2, 0), las = 1,
      mgp = c(2.5, 0.6, 0), tcl = -0.3, bty = "l", cex.axis = 0.85, cex.lab = 0.95)

  # (a) Linear scale: histogram drawn as a filled step outline
  x_step <- rep(h$breaks[c(shown, FALSE) | c(FALSE, shown)], each = 2)
  y_step <- c(0, rep(h$density[shown], each = 2), 0)
  plot(NA, xlim = c(-x_max, x_max), ylim = c(0, max(h$density[shown], dnorm(0)) * 1.05),
       xlab = "Z-score", ylab = "Density", xaxs = "i", yaxs = "i")
  polygon(x_step, y_step, col = col_fill, border = col_data, lwd = 0.6)
  curve(dnorm(x), from = -x_max, to = x_max, n = 1001, lty = 2, lwd = 2, col = col_normal, add = TRUE)
  legend("topleft", legend = c("Z-scores", "N(0, 1)"), inset = c(0.02, 0),
         fill = c(col_fill, NA), border = c(col_data, NA), lty = c(NA, 2), lwd = c(NA, 2),
         col = c(NA, col_normal), merge = TRUE, bty = "n", cex = 0.8)
  mtext("a  Linear scale", side = 3, adj = 0, line = 0.5, font = 2)

  # (b) Log scale: one point per non-empty bin; floor = density of one SNP in a bin
  nonzero <- shown & h$density > 0
  one_snp <- 1 / (length(z) * binwidth)
  y_lim <- c(one_snp / 2, max(h$density, dnorm(0)) * 1.5)
  plot(h$mids[nonzero], h$density[nonzero], log = "y", pch = 16, cex = 0.35, col = col_data,
       xlim = c(-x_max, x_max), ylim = y_lim, xaxs = "i", yaxt = "n",
       xlab = "Z-score", ylab = "Density (log scale)")
  powers <- seq(ceiling(log10(y_lim[1])), floor(log10(y_lim[2])))
  axis(2, at = 10^powers, labels = parse(text = paste0("10^", powers)))
  curve(dnorm(x), from = -x_max, to = x_max, n = 1001, lty = 2, lwd = 2, col = col_normal, add = TRUE)
  mtext("b  Log scale (tails)", side = 3, adj = 0, line = 0.5, font = 2)

  # Title and summary statistics
  mtext(title, outer = TRUE, line = 0.6, font = 2, cex = 1.15)
  stats <- sprintf("n = %s   |   mean = %.3f   |   SD = %.3f   |   \u03bb = %.3f   |   bin width = %g",
                   format(length(z), big.mark = ","), round(mean(z), 3) + 0,  # + 0 avoids "-0.000"
                   sd(z), lambda, binwidth)
  if (n_outside) stats <- paste0(stats, sprintf("   |   %s SNPs with |Z| > %g not shown",
                                                format(n_outside, big.mark = ","), x_max))
  mtext(stats, side = 1, outer = TRUE, line = 0.8, cex = 0.8, col = "grey25")

  data.table(n = length(z), mean = mean(z), sd = sd(z), lambda = lambda,
             n_abs_z_gt_5.45 = sum(abs(z) > 5.45))
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
dir.create(cfg$qc_plot_dir, recursive = TRUE, showWarnings = FALSE)
z_summary <- list()

for (i in seq_len(cfg$n_components)) {
  zscores <- read_one_numeric_column(sprintf(cfg$zscore_pattern, i), "Z-scores")

  # Safeguard: rows are matched by position, so the counts must be identical
  if (length(zscores) != nrow(catalogue)) {
    stop(sprintf("Component %d: %d Z-scores, but %d included SNPs in the catalogue.",
                 i, length(zscores), nrow(catalogue)), call. = FALSE)
  }

  component <- copy(catalogue)
  component[, `:=`(BETA = zscores, P = 2 * pnorm(-abs(zscores)))]
  fwrite(component, sprintf(cfg$merged_pattern, i))

  z_summary[[i]] <- data.table(
    component = i,
    n_missing = sum(!is.finite(zscores)),
    plot_z_histogram(zscores, sprintf("Component %d: Z-scores vs N(0, 1)", i),
                     file.path(cfg$qc_plot_dir, sprintf("component%d_z_histogram.png", i)),
                     cfg$z_hist_binwidth)
  )
}
fwrite(rbindlist(z_summary), file.path(cfg$qc_plot_dir, "component_z_summary.csv"))
message("Wrote ", cfg$n_components, " component files: ", sprintf(cfg$merged_pattern, 1), ", ...")
message("QC histograms and summary in ", cfg$qc_plot_dir)
