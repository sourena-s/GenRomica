# =============================================================================
# Trajectory figure for one GPS x repeated outcome (used by step 7)
# =============================================================================
# Shows the existing mixed model M2 of step 5 (R/models.R, fit_lmm_one_score)
# against the data. Nothing is re-modelled here:
#   (a) outcome against age. Faint lines: observed trajectories of a sample of
#       individuals. Points with 95% CI: observed mean per wave in GPS tertiles
#       (descriptive). Curves with 95% CI bands: M2 fixed-effect predictions for
#       GPS at -1 SD, the mean and +1 SD, with covariates at their sample means.
#   (b) M2 residuals (observed - fitted, including the individual's random
#       intercept) per wave. A good fit is centred on 0 at every wave.

suppressPackageStartupMessages(library(patchwork))

gps_levels <- c("-1 SD", "Mean", "+1 SD")
gps_colours <- c("-1 SD" = "#1c5cab", "Mean" = "#52514e", "+1 SD" = "#c22f2f")   # diverging poles + neutral
tertile_to_level <- c("Low" = "-1 SD", "Middle" = "Mean", "High" = "+1 SD")

# M2 fixed-effect predictions over the observed age range, for GPS = -1, 0, +1 SD
predict_m2 <- function(m2, d) {
  X_fit <- getME(m2, "X")
  beta <- fixef(m2)
  V <- as.matrix(vcov(m2))
  centre <- mean(d$age - d$age_c)              # centring age used when the model was fitted

  grid <- CJ(age = seq(min(d$age), max(d$age), length.out = 120), gps = c(-1, 0, 1))
  X <- matrix(colMeans(X_fit), nrow(grid), ncol(X_fit), byrow = TRUE,
              dimnames = list(NULL, colnames(X_fit)))   # covariates at their sample means
  age_c <- grid$age - centre
  X[, "(Intercept)"] <- 1
  X[, "age_c"] <- age_c
  X[, "age2"]  <- age_c^2
  X[, "score"] <- grid$gps
  X[, interaction_term(m2, "score", "age_c")] <- grid$gps * age_c
  X[, interaction_term(m2, "score", "age2")]  <- grid$gps * age_c^2

  grid[, fit := as.vector(X %*% beta)]
  grid[, se := sqrt(rowSums((X %*% V) * X))]
  grid[, `:=`(lower = fit - 1.96 * se, upper = fit + 1.96 * se,
              level = factor(gps_levels[gps + 2], levels = gps_levels))]
  grid
}

plot_trajectory <- function(fit, score_name, outcome_name, p_text = NULL, n_lines = 100, seed = 1) {
  m2 <- fit$m2
  d <- copy(fit$data)
  d[, residual := residuals(m2)]

  # GPS tertiles over individuals (one value per individual)
  per_id <- d[!duplicated(idc), .(idc, score)]
  per_id[, tertile := cut(score, quantile(score, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                          labels = c("Low", "Middle", "High"))]
  d[per_id, on = "idc", tertile := i.tertile]
  d[, level := factor(tertile_to_level[as.character(tertile)], levels = gps_levels)]

  observed <- d[, .(age = mean(age), mean = mean(outcome), se = sd(outcome) / sqrt(.N), n = .N),
                by = .(wave, level)]
  dodge <- c("-1 SD" = -0.12, "Mean" = 0, "+1 SD" = 0.12)       # small horizontal offset per tertile
  observed[, age := age + dodge[as.character(level)]]

  set.seed(seed)
  shown <- sample(unique(d$idc), min(n_lines, uniqueN(d$idc)))
  curves <- predict_m2(m2, d)

  y_label <- paste0(label_outcome(outcome_name), " (z-score)")
  panel_a <- ggplot() +
    geom_line(data = d[idc %in% shown], aes(age, outcome, group = idc),
              colour = "#b5b4af", alpha = 0.25, linewidth = 0.25) +
    geom_ribbon(data = curves, aes(age, ymin = lower, ymax = upper, fill = level), alpha = 0.18) +
    geom_line(data = curves, aes(age, fit, colour = level), linewidth = 0.9) +
    geom_errorbar(data = observed, aes(age, ymin = mean - 1.96 * se, ymax = mean + 1.96 * se, colour = level),
                  width = 0, linewidth = 0.5) +
    geom_point(data = observed, aes(age, mean, fill = level), shape = 21, colour = "white", size = 2.4, stroke = 0.4) +
    scale_colour_manual(values = gps_colours, name = label_score(score_name), drop = FALSE) +
    scale_fill_manual(values = gps_colours, name = label_score(score_name), drop = FALSE) +
    labs(x = "Age (years)", y = y_label, subtitle = "a   Observed data and M2 fit") +
    coord_cartesian(ylim = range(c(curves$lower, curves$upper, observed$mean + c(-2, 2) * observed$se,
                                   quantile(d$outcome, c(0.02, 0.98))))) +
    theme_figure() +
    theme(panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3),
          plot.subtitle = element_text(face = "bold", size = 10),
          legend.position = "bottom")

  wave_ages <- d[, .(age = median(age)), by = wave][order(age)]
  d[, wave := factor(wave, levels = wave_ages$wave,
                     labels = sprintf("%s\n(~%.0f y)", wave_ages$wave, wave_ages$age))]
  panel_b <- ggplot(d, aes(wave, residual)) +
    geom_hline(yintercept = 0, colour = ink_muted, linewidth = 0.4) +
    geom_violin(fill = "#e4e3df", colour = NA, scale = "width", width = 0.8) +
    geom_boxplot(width = 0.18, outlier.shape = NA, colour = ink, fill = "white", linewidth = 0.35) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "#c22f2f", colour = "white") +
    labs(x = "Wave", y = "M2 residual (z-score)", subtitle = "b   M2 residuals by wave") +
    theme_figure() +
    theme(panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3),
          plot.subtitle = element_text(face = "bold", size = 10))

  title <- sprintf("%s and %s", label_score(score_name), label_outcome(outcome_name))
  caption <- paste0(
    "(a) Grey lines: observed trajectories of ", length(shown), " randomly chosen individuals. ",
    "Points: observed mean ± 95% CI per wave in GPS tertiles.\n",
    "Curves: M2 fixed-effect predictions (step 5) for the GPS at −1 SD, the mean and +1 SD, ",
    "covariates at sample means; bands: 95% CI.\n",
    "(b) Residuals of the same M2 fit, including each individual's random intercept. ",
    "Box: median and IQR; red diamond: mean.")

  (panel_a | panel_b) +
    plot_layout(widths = c(2, 1)) +
    plot_annotation(title = title, subtitle = p_text, caption = caption,
                    theme = theme_figure() + theme(plot.subtitle = element_text(colour = ink_muted, size = 9)))
}
