# =============================================================================
# Shared figure style: labels, colours and theme for steps 6 and 7
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# ---- Labels, order and colours -----------------------------------------------

domain_order <- c("CBCL", "SES", "Diet", "Online", "MRI", "DTI")
# Domain names shown in the figures (data keep the short names)
domain_display <- c(CBCL = "CBCL", SES = "SES", Diet = "Diet", Online = "Online",
                    MRI = "T1w MRI", DTI = "dMRI")
domain_colours <- c("CBCL" = "#2a78d6", "T1w MRI" = "#eb6834", "dMRI" = "#1baf7a",  # colours follow the
                    "SES" = "#eda100", "Diet" = "#e87ba4", "Online" = "#008300")    # domain, not its position
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
col_nominal <- "#6250d6"; col_fdr <- "#e34948"   # single test: purple; FDR (all tests): red

label_outcome <- function(x) ifelse(x %in% names(outcome_labels), outcome_labels[x], x)

# Score names in the figures: GenomICA component scores C1, C2, ... are shown as
# GPS1, GPS2, ... (GenomICA polygenic scores). Other score names are unchanged.
label_score <- function(x) sub("^C([0-9]+)$", "GPS\\1", x)

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

save_png <- function(plot, name, width, height, dir = fig_dir) {
  ggsave(file.path(dir, name), plot, width = width, height = height, dpi = 300,
         device = "png", type = "cairo", bg = "white")
  message("Wrote ", file.path(dir, name))
}



# Hypotheses, in figure order. `signed`: the test has one coefficient with a
# sign (drawn as up/down triangles and as a diverging heatmap of the estimate).
hypothesis_info <- data.table(
  hypothesis = c("H1", "H2_lin_interact", "H2_quad_interact", "H2_lin_and_quad_interact"),
  title = c("H1: association, one score at a time",
            "H2_lin_interact: linear change with age (GPS \u00d7 age)",
            "H2_quad_interact: curved change with age (GPS \u00d7 age\u00b2, beyond linear)",
            "H2_lin_and_quad_interact: any change with age (GPS \u00d7 age and GPS \u00d7 age\u00b2)"),
  signed = c(TRUE, TRUE, TRUE, FALSE),
  estimate_label = c("Standardised\nestimate", "GPS \u00d7 age\nestimate", "GPS \u00d7 age\u00b2\nestimate", NA)
)
