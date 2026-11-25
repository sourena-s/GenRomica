# =============================================================================
# Step 6 - Figures of the association results
# =============================================================================
# Reads results.csv from step 5 (cfg$results_dir) and writes, for every score
# set tested one score at a time (e.g. PCS = C1-C20), to cfg$results_dir/figures:
#
#   <set>_phewas.png        PheWAS-style plot: -log10(p) of every score x
#                           outcome, outcomes grouped and coloured by domain.
#                           (a) H1 (triangles point up/down with the sign of
#                           the estimate), (b) H2 (age interaction).
#   <set>_heatmap_H1.png    score x outcome heatmap of the standardised
#                           estimate (diverging colour scale)
#   <set>_heatmap_H2.png    score x outcome heatmap of -log10(p) for the
#                           GPS x age interaction (sequential colour scale)
#
# Significance marks: * p < 0.05; ** FDR-corrected p < 0.05.
# Jointly tested sets (*_combined) have one result per outcome and are not
# plotted here.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})
source("R/utils.R")
cfg <- load_config()

results <- read_table(file.path(cfg$results_dir, "results.csv"))
require_columns(results, c("hypothesis", "estimate", "p_value", "p_fdr", "predictor",
                           "score_group", "domain", "outcome"),
                "results.csv", load_command = sprintf('d <- data.table::fread("%s")',
                                                      file.path(cfg$results_dir, "results.csv")))
fig_dir <- file.path(cfg$results_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ---- Labels, order and colours -----------------------------------------------

domain_order <- c("CBCL", "SES", "Diet", "Online", "MRI", "DTI")
domain_colours <- c(CBCL = "#2a78d6", MRI = "#eb6834", DTI = "#1baf7a",    # fixed order of the
                    SES = "#eda100", Diet = "#e87ba4", Online = "#008300")  # reference palette
outcome_labels <- c(
  anx = "Anxious/depressed", som = "Somatic complaints", wit = "Withdrawn",
  sop = "Social problems", tho = "Thought problems", att = "Attention problems",
  rul = "Rule-breaking", agg = "Aggressive behaviour",
  TotalGrayVol = "Grey-matter volume", cortical_thickness = "Cortical thickness",
  mean_FA_genr = "Mean FA", mean_MD_genr = "Mean MD",
  INCOME5 = "Household income", EDUCM = "Maternal education", EDUCP = "Paternal education",
  DietScore_8y = "Diet quality (8y)", onlineweek_child_F17 = "Online activity (17y)"
)
ink <- "#0b0b0b"; ink_muted <- "#52514e"; grid_col <- "#e4e3df"
# threshold lines use palette slots 7-8, which no domain uses
col_fdr <- "#6250d6"; col_nominal <- "#e34948"

label_outcome <- function(x) ifelse(x %in% names(outcome_labels), outcome_labels[x], x)

theme_figure <- function() {
  theme_classic(base_size = 10) %+replace% theme(
    text = element_text(colour = ink),
    axis.text = element_text(colour = ink_muted),
    axis.line = element_line(colour = ink_muted, linewidth = 0.3),
    axis.ticks = element_line(colour = ink_muted, linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 11, hjust = 0, margin = margin(b = 6)),
    plot.caption = element_text(colour = ink_muted, size = 8, hjust = 0),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", colour = ink, size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )
}

save_png <- function(plot, name, width, height) {
  ggsave(file.path(fig_dir, name), plot, width = width, height = height, dpi = 300,
         device = "png", type = "cairo", bg = "white")
  message("Wrote ", file.path(fig_dir, name))
}

stars <- function(p, p_fdr) fifelse(p_fdr < 0.05, "**", fifelse(p < 0.05, "*", ""))


# ---- PheWAS plot ---------------------------------------------------------------

phewas_plot <- function(d, set_name) {
  d <- copy(d)
  d[, domain := factor(domain, levels = intersect(domain_order, unique(domain)))]
  d[, minus_log10_p := -log10(p_value)]
  d[, direction := fcase(hypothesis == "H2", "Age interaction (H2)",
                         estimate >= 0, "Positive estimate",
                         default = "Negative estimate")]
  d[, direction := factor(direction, levels = c("Positive estimate", "Negative estimate", "Age interaction (H2)"))]
  d[, panel := factor(fifelse(hypothesis == "H1", "a   H1: association, one score at a time",
                              "b   H2: change of the association with age (GPS × age)"))]

  # x position: outcomes in domain order (1, 2, ...); the scores of one outcome
  # are spread over a narrow band around that position so they do not overlap
  scores <- unique(d$predictor)
  scores <- scores[order(as.integer(gsub("\\D", "", scores)))]
  d[, outcome_rank := match(outcome, names(outcome_labels))]
  d[, outcome_x := match(outcome, unique(.SD[order(domain, outcome_rank)]$outcome)), by = panel]
  d[, x := outcome_x + (match(predictor, scores) - (length(scores) + 1) / 2) * 0.6 / length(scores)]

  # nominal and FDR thresholds (FDR: largest p with FDR-corrected p < 0.05)
  thresholds <- d[, .(fdr = if (any(p_fdr < 0.05)) max(p_value[p_fdr < 0.05]) else NA_real_), by = panel]

  # alternating background bands mark the domains; domain names sit above them
  bands <- d[, .(xmin = min(outcome_x) - 0.5, xmax = max(outcome_x) + 0.5), by = .(panel, domain)]
  bands[, shade := seq_len(.N) %% 2 == 1, by = panel]
  # narrow domains (one outcome) alternate in height so their names do not overlap
  bands[, label_raise := (xmax - xmin) <= 1 & seq_len(.N) %% 2 == 0, by = panel]
  # each panel has its own y range; domain names sit just above its highest point
  bands[d[, .(y_top = max(c(minus_log10_p, -log10(0.05))) * 1.12), by = panel], on = "panel", y_top := i.y_top]
  axis_labels <- unique(d[, .(panel, outcome, outcome_x)])
  hits <- d[p_fdr < 0.05]

  ggplot(d, aes(x, minus_log10_p)) +
    geom_rect(data = bands[shade == TRUE], aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "#f4f3f0") +
    geom_blank(data = bands, aes(x = xmin, y = y_top * 1.12), inherit.aes = FALSE) +
    geom_text(data = bands, aes(x = (xmin + xmax) / 2, y = y_top * fifelse(label_raise, 1.08, 1), label = domain),
              inherit.aes = FALSE, size = 3, fontface = "bold", colour = ink_muted, vjust = 0) +
    geom_hline(yintercept = -log10(0.05), linetype = "22", colour = col_nominal, linewidth = 0.5) +
    geom_hline(data = thresholds[!is.na(fdr)], aes(yintercept = -log10(fdr)),
               linetype = "22", colour = col_fdr, linewidth = 0.5) +
    annotate("text", x = Inf, y = -log10(0.05), label = "p = 0.05 (single test)",
             colour = col_nominal, size = 2.7, hjust = -0.05, vjust = 0.4, fontface = "bold") +
    geom_text(data = thresholds[!is.na(fdr)], aes(x = Inf, y = -log10(fdr)),
              label = "FDR 5% (all tests)", inherit.aes = FALSE,
              colour = col_fdr, size = 2.7, hjust = -0.05, vjust = 0.4, fontface = "bold") +
    geom_point(aes(fill = domain, shape = direction), colour = "white", stroke = 0.3, size = 2.3) +
    geom_text_repel(data = hits, aes(label = predictor), size = 2.6, colour = ink,
                    min.segment.length = 0, segment.colour = ink_muted, segment.size = 0.25,
                    box.padding = 0.3, seed = 1) +
    geom_text(data = axis_labels, aes(x = outcome_x, y = -Inf, label = label_outcome(outcome)),
              inherit.aes = FALSE, angle = 40, hjust = 1, vjust = 1.2, size = 2.7, colour = ink_muted) +
    facet_wrap(~ panel, ncol = 1, scales = "free") +
    scale_x_continuous(breaks = NULL, expand = expansion(add = 0.1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = domain_colours, name = "Domain") +
    scale_shape_manual(values = c("Positive estimate" = 24, "Negative estimate" = 25, "Age interaction (H2)" = 21),
                       name = NULL, drop = FALSE) +
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 4, colour = NA), order = 1),
           shape = guide_legend(override.aes = list(fill = ink_muted, size = 2.6), order = 2)) +
    labs(x = NULL, y = expression(-log[10](italic(p))),
         title = sprintf("%s: GenomICA polygenic scores and outcomes in Generation R", set_name),
         caption = paste0("Each point is one score × outcome test (one point per component for every outcome). ",
                          "Red dashed line: p = 0.05 for a single test.\nPurple dashed line: FDR 5% threshold across all ",
                          "scores and outcomes of the panel (Benjamini–Hochberg). Labelled points pass FDR.")) +
    coord_cartesian(clip = "off") +
    theme_figure() +
    theme(panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3),
          strip.text = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 14)),
          panel.spacing = unit(5.5, "lines"),
          plot.margin = margin(10, 120, 10, 70),
          legend.position = "bottom", legend.box = "horizontal",
          legend.box.margin = margin(t = 85),
          axis.title.y = element_text(angle = 90, margin = margin(r = 34)),
          plot.caption = element_text(colour = ink_muted, size = 8, hjust = 0, margin = margin(t = 12)),
          axis.line.x = element_blank())
}


# ---- Heatmaps -------------------------------------------------------------------

heatmap_plot <- function(d, hyp, set_name) {
  d <- copy(d[hypothesis == hyp])
  d[, domain := factor(domain, levels = intersect(domain_order, unique(domain)))]
  d[, outcome_label := factor(label_outcome(outcome),
                              levels = unique(label_outcome(d[order(domain, match(outcome, names(outcome_labels)))]$outcome)))]
  scores <- unique(d$predictor)
  d[, predictor := factor(predictor, levels = rev(scores[order(as.integer(gsub("\\D", "", scores)))]))]
  d[, mark := stars(p_value, p_fdr)]

  if (hyp == "H1") {
    d[, value := estimate]
    limit <- max(abs(d$value), na.rm = TRUE)
    fill_scale <- scale_fill_gradient2(low = "#1c5cab", mid = "#f0efec", high = "#c22f2f", midpoint = 0,
                                       limits = c(-limit, limit), name = "Standardised\nestimate")
    title <- sprintf("%s, H1: association of each score with each outcome", set_name)
    dark <- abs(d$value) > 0.55 * limit
  } else {
    d[, value := -log10(p_value)]
    fill_scale <- scale_fill_gradientn(colours = c("#f0efec", "#b7d3f6", "#6da7ec", "#2a78d6", "#1c5cab", "#0d366b"),
                                       limits = c(0, max(d$value)), name = expression(-log[10](italic(p))))
    title <- sprintf("%s, H2: change of the association with age (GPS × age)", set_name)
    dark <- d$value > 0.55 * max(d$value)
  }
  d[, mark_colour := fifelse(dark, "white", ink)]

  ggplot(d, aes(outcome_label, predictor, fill = value)) +
    geom_tile(colour = "#3a3f44", linewidth = 0.25) +
    geom_text(aes(label = mark, colour = mark_colour), size = 6.8, vjust = 0.78, fontface = "bold") +
    scale_colour_identity() +
    fill_scale +
    facet_grid(. ~ domain, scales = "free_x", space = "free_x") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = title,
         caption = "* p < 0.05;  ** FDR-corrected p < 0.05 (Benjamini–Hochberg within score set and hypothesis).") +
    theme_figure() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          plot.margin = margin(8, 8, 8, 40),
          axis.line = element_blank(), axis.ticks = element_blank(),
          panel.spacing = unit(0.35, "lines"),
          strip.text = element_text(face = "bold", size = 9, margin = margin(b = 4)),
          legend.key.height = unit(1.1, "lines"))
}


# ---- Draw every score set ---------------------------------------------------------

for (set_name in setdiff(unique(results$score_group), grep("_combined$", unique(results$score_group), value = TRUE))) {
  d <- results[score_group == set_name]
  n_scores <- uniqueN(d$predictor)
  save_png(phewas_plot(d, set_name), sprintf("%s_phewas.png", set_name), width = 10, height = 8)
  if (any(d$hypothesis == "H1"))
    save_png(heatmap_plot(d, "H1", set_name), sprintf("%s_heatmap_H1.png", set_name),
             width = 9, height = 2 + 0.22 * n_scores)
  if (any(d$hypothesis == "H2"))
    save_png(heatmap_plot(d, "H2", set_name), sprintf("%s_heatmap_H2.png", set_name),
             width = 7.2, height = 2 + 0.22 * n_scores)
}
