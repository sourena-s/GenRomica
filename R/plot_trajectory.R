# =============================================================================
# Trajectory figure for one GPS x repeated outcome (used by step 7)
# =============================================================================
# Shows the existing mixed model M2 of step 5 (R/models.R, fit_lmm_one_score)
# against the data. Nothing is re-modelled here:
#   Outcome against age. Points: every observation (scatter per wave).
#   Boxplots: outcome per wave, at the wave's median age. Curves with 95% CI
#   bands on top: M2 fixed-effect predictions for GPS at -1 SD, the mean and
#       +1 SD, with covariates at their sample means.

gps_levels <- c("-1 SD", "Mean", "+1 SD")
gps_colours <- c("-1 SD" = "#1c5cab", "Mean" = "#52514e", "+1 SD" = "#c22f2f")   # diverging poles + neutral

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

plot_trajectory <- function(fit, score_name, outcome_name, p_text = NULL) {
  m2 <- fit$m2
  d <- copy(fit$data)
  curves <- predict_m2(m2, d)

  # each wave's boxplot sits at that wave's median age
  wave_ages <- d[, .(wave_age = median(age)), by = wave][order(wave_age)]
  d[wave_ages, on = "wave", wave_age := i.wave_age]
  box_width <- 0.35 * min(diff(wave_ages$wave_age), 4)

  y_label <- paste0(label_outcome(outcome_name), " (z-score)")
  panel_a <- ggplot() +
    geom_point(data = d, aes(age, outcome), colour = "#9c9b96", alpha = 0.25, size = 0.6, stroke = 0) +
    geom_boxplot(data = d, aes(wave_age, outcome, group = wave), width = box_width,
                 outlier.shape = NA, fill = "white", alpha = 0.75, colour = ink, linewidth = 0.35) +
    geom_ribbon(data = curves, aes(age, ymin = lower, ymax = upper, fill = level), alpha = 0.22) +
    geom_line(data = curves, aes(age, fit, colour = level), linewidth = 1) +
    scale_colour_manual(values = gps_colours, name = label_score(score_name), drop = FALSE) +
    scale_fill_manual(values = gps_colours, name = label_score(score_name), drop = FALSE) +
    scale_x_continuous(sec.axis = dup_axis(breaks = wave_ages$wave_age, labels = wave_ages$wave, name = NULL)) +
    labs(x = "Age (years)", y = y_label) +
    theme_figure() +
    theme(panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3),
          axis.line.x.top = element_blank(), axis.ticks.x.top = element_blank(),
          axis.text.x.top = element_text(face = "bold", colour = ink_muted),
          legend.position = "bottom")

  panel_a +
    labs(title = sprintf("%s and %s", label_score(score_name), label_outcome(outcome_name)),
         subtitle = p_text,
         caption = paste0(
           "Grey points: every observation at its age. Boxes: outcome per wave (median, IQR, whiskers 1.5 IQR), ",
           "at the wave's median age.\n",
           "Curves: M2 fixed-effect predictions (step 5) for the GPS at \u22121 SD, the mean and +1 SD, ",
           "covariates at sample means; bands: 95% CI.")) +
    theme(plot.subtitle = element_text(colour = ink_muted, size = 9, face = "plain"))
}
