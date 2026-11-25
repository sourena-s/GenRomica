# =============================================================================
# Step 7 - Trajectory figures: the M2 mixed model against the data, by wave
# =============================================================================
# For each selected GPS x repeated outcome, refits the step-5 mixed model M2
# with the same code (R/models.R, fit_lmm_one_score) and draws
# cfg$results_dir/figures/trajectory_<score>_<outcome>.png (see R/plot_trajectory.R).
#
# Pairs drawn: cfg$trajectory_pairs (data.table with columns score, outcome) if
# set; otherwise the pairs with FDR-corrected p < 0.05 for H1 or H2 in
# results.csv, or the cfg$max_trajectory_plots pairs with the smallest p if
# none pass FDR.

source("R/utils.R")
cfg <- load_config()
source("R/plot_style.R")
source("R/plot_trajectory.R")

gpcs <- cfg$genetic_pcs
genr <- load_genr_with_scores(cfg)
long_data <- make_long_data(cfg, genr)
source("R/models.R")

results <- read_table(file.path(cfg$results_dir, "results.csv"))
fig_dir <- file.path(cfg$results_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
max_plots <- if (is.null(cfg$max_trajectory_plots)) 6 else cfg$max_trajectory_plots

# ---- Which pairs to draw ---------------------------------------------------------
single_sets <- names(cfg$score_sets)
long_results <- results[score_group %in% single_sets & domain %in% unique(cfg$longitudinal$domain)]

if (!is.null(cfg$trajectory_pairs)) {
  pairs <- cfg$trajectory_pairs
} else {
  ranked <- long_results[, .(p_min = min(p_value), fdr_min = min(p_fdr)), by = .(score = predictor, outcome)][order(p_min)]
  pairs <- if (any(ranked$fdr_min < 0.05)) ranked[fdr_min < 0.05] else ranked
  pairs <- head(pairs, max_plots)
}
message("Drawing ", nrow(pairs), " trajectory figure(s)")

format_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3))

# ---- Draw ------------------------------------------------------------------------
for (k in seq_len(nrow(pairs))) {
  score_k <- pairs$score[k]
  outcome_k <- pairs$outcome[k]
  x <- long_data[outcome_name == outcome_k]

  fit <- fit_lmm_one_score(x, score_k, return_model = TRUE)
  if (is.null(fit)) { message("Skipped ", score_k, " x ", outcome_k, ": too little data"); next }

  r <- long_results[predictor == score_k & outcome == outcome_k]
  p_text <- paste(vapply(c("H1", "H2"), function(h) {
    row <- r[hypothesis == h]
    if (!nrow(row)) return(sprintf("%s: not in results.csv", h))
    sprintf("%s: p = %s (FDR p = %s)", h, format_p(row$p_value[1]), format_p(row$p_fdr[1]))
  }, ""), collapse = "     ")

  save_png(plot_trajectory(fit$fit, score_k, outcome_k, p_text),
           sprintf("trajectory_%s_%s.png", label_score(score_k), outcome_k), width = 8.5, height = 6, dir = fig_dir)
}
