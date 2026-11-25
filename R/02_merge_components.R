# =============================================================================
# Step 2 - Attach ICA component Z-scores to the variant catalogue
# =============================================================================
# Inputs (from the ICA study):
#   variant_catalogue.txt    no header, 7 columns:
#                              CHR, BP, A1 (other allele), A2 (effect allele),
#                              SNP (rsID), CHR:BP:A, include flag (1 = used in ICA)
#   melodic_oIc_<i>.txt      one column, no header: weight of each SNP on component i,
#                            from MELODIC melodic_oIC. Treated as a Z-score here, but
#                            not Z-transformed (open issue: README, concern B6)
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
source("R/plot_style.R")
suppressPackageStartupMessages(library(patchwork))
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
  h <- hist(z, breaks = breaks, plot = FALSE)        # binning only; drawn with ggplot2 below
  bins <- data.table(mid = h$mids, density = h$density)
  x_max <- min(10, max(5, ceiling(quantile(abs(z), 0.9999))))
  n_outside <- sum(abs(z) > x_max)
  lambda <- median(z^2) / qchisq(0.5, df = 1)
  bins <- bins[abs(mid) <= x_max]
  normal <- data.table(mid = seq(-x_max, x_max, length.out = 1001))[, density := dnorm(mid)]
  one_snp <- 1 / (length(z) * binwidth)               # density of a bin holding one SNP

  col_fill <- "#c9ced6"; col_data <- "#3a3f44"; col_normal <- "#1f5fbf"   # N(0, 1): blue dashed
  base <- function() list(
    geom_line(data = normal, aes(mid, density, linetype = "N(0, 1)"), colour = col_normal, linewidth = 0.8),
    scale_linetype_manual(values = c("N(0, 1)" = "22"), name = NULL),
    scale_x_continuous(limits = c(-x_max, x_max), expand = c(0, 0)),
    labs(x = "Z-score"),
    theme_figure(),
    theme(panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3),
          plot.subtitle = element_text(face = "bold", size = 10))
  )

  linear <- ggplot(bins, aes(mid, density)) +
    geom_col(width = binwidth, fill = col_fill, colour = NA) +
    geom_step(aes(x = mid - binwidth / 2), colour = col_data, linewidth = 0.25, direction = "hv") +
    base() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(y = "Density", subtitle = "a   Linear scale") +
    theme(legend.position = c(0.85, 0.9))

  powers <- seq(ceiling(log10(one_snp / 2)), 0)
  log_scale <- ggplot(bins[density > 0], aes(mid, density)) +
    geom_point(colour = col_data, size = 0.5) +
    base() +
    scale_y_log10(limits = c(one_snp / 2, 1), breaks = 10^powers,
                  labels = function(b) parse(text = paste0("10^", log10(b)))) +
    labs(y = "Density (log scale)", subtitle = "b   Log scale (tails)") +
    theme(legend.position = "none")

  stats <- sprintf("n = %s   |   mean = %.3f   |   SD = %.3f   |   λ = %.3f   |   bin width = %g",
                   format(length(z), big.mark = ","), round(mean(z), 3) + 0,  # + 0 avoids "-0.000"
                   sd(z), lambda, binwidth)
  if (n_outside) stats <- paste0(stats, sprintf("   |   %s SNPs with |Z| > %g not shown",
                                                format(n_outside, big.mark = ","), x_max))

  figure <- (linear | log_scale) +
    plot_annotation(title = title, caption = stats,
                    theme = theme_figure() + theme(plot.caption = element_text(hjust = 0.5, colour = "grey25")))
  suppressWarnings(ggsave(path, figure, width = 10, height = 4.3, dpi = 300, device = "png",
                          type = "cairo", bg = "white"))

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
