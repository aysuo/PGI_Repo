library(data.table)
library(ggplot2)

# ---------------------------------------------------------------------------
# Usage: Rscript plot_loglik_margin.R <classified_output.tsv>
# Visualises output from ancestry_classifier.R (dual-filter version).
#
# Produces four plots:
#   1. margin_overall.png        — delta margin distribution, coloured by
#                                  rejection reason
#   2. margin_by_ancestry.png    — delta margin faceted by assigned ancestry,
#                                  with per-ancestry delta threshold line
#   3. score_floor.png           — best log-likelihood score per sample,
#                                  faceted by ancestry, with percentile floor
#                                  line (Filter 1 diagnostic)
#   4. rejection_summary.png     — stacked bar: assigned vs rejection reason
#                                  per ancestry
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript plot_loglik_margin.R <classified_output.tsv>")

data <- fread(args[1])
output <- args[2]

# Must match values used in ancestry_classifier.R
percentile_floor <- 0.9999999
delta_threshold  <- as.numeric(args[4])


# Keep only rows that were classified
classified <- data[!is.na(loglik_margin)]
if (nrow(classified) == 0) stop("No classified samples found in the file.")

ancestries <- c("EUR", "AFR", "EAS", "AMR", "SAS")

anc_colours <- c(
  EUR       = "#378ADD",
  AFR       = "#1D9E75",
  EAS       = "#D85A30",
  AMR       = "#7F77DD",
  SAS       = "#BA7517",
  ambiguous = "#B4B2A9",
  unassigned = "#CCCCCC"
)

reason_colours <- c(
  assigned               = "#534AB7",
  below_percentile_floor = "#E24B4A",
  low_delta_margin       = "#F5A623",
  no_valid_ancestry      = "#AAAAAA"
)

# Derived convenience columns
classified[, status := fcase(
  closest_ancestry %in% ancestries,        "assigned",
  rejection_reason == "below_percentile_floor", "below_percentile_floor",
  rejection_reason == "low_delta_margin",        "low_delta_margin",
  rejection_reason == "no_valid_ancestry",       "no_valid_ancestry",
  default = "assigned"
)]

assigned <- classified[closest_ancestry %in% ancestries]

# ---------------------------------------------------------------------------
# Plot 1 — Overall delta margin distribution, coloured by rejection reason
# ---------------------------------------------------------------------------

p1 <- ggplot(classified, aes(x = loglik_margin, fill = status)) +
  geom_histogram(binwidth = 1, colour = "white", linewidth = 0.2) +
  geom_vline(
    xintercept = delta_threshold,
    linetype = "dashed", colour = "#333333", linewidth = 0.7
  ) +
  annotate(
    "text", x = delta_threshold + 0.3, y = Inf,
    label = paste0("delta threshold = ", delta_threshold),
    hjust = 0, vjust = 1.5, size = 3, colour = "#333333"
  ) +
  scale_fill_manual(
    values = reason_colours,
    labels = c(
      assigned               = "Assigned",
      below_percentile_floor = "Rejected: below score floor",
      low_delta_margin       = "Rejected: low delta margin",
      no_valid_ancestry      = "Rejected: no valid ancestry"
    ),
    name = NULL
  ) +
  labs(
    title    = "Delta margin distribution (all samples)",
    subtitle = "Colour shows which filter caused rejection",
    x        = "loglik_margin  (best − second-best log-likelihood)",
    y        = "Count"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave(paste0(output, "_margin_overall.png"), plot = p1, width = 9, height = 4.5, dpi = 150)
cat("Saved: margin_overall.png\n")

# ---------------------------------------------------------------------------
# Plot 2 — Delta margin faceted by assigned ancestry
# One dashed line per panel at delta_threshold (same for all, but easy to
# extend to per-ancestry values if desired)
# ---------------------------------------------------------------------------

if (nrow(assigned) > 0) {
  p2 <- ggplot(assigned, aes(x = loglik_margin, fill = closest_ancestry)) +
    geom_histogram(binwidth = 1, colour = "white", linewidth = 0.2) +
    geom_vline(
      xintercept = delta_threshold,
      linetype = "dashed", colour = "#333333", linewidth = 0.7
    ) +
    facet_wrap(~closest_ancestry, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = anc_colours, guide = "none") +
    labs(
      title    = "Delta margin by assigned ancestry (assigned samples only)",
      subtitle = paste0("Dashed line = delta threshold (", delta_threshold, ")"),
      x        = "loglik_margin",
      y        = "Count"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(paste0(output,"_margin_by_ancestry.png"), plot = p2, width = 9, height = 5, dpi = 150)
  cat("Saved: margin_by_ancestry.png\n")
}

# ---------------------------------------------------------------------------
# Plot 3 — Best log-likelihood score per sample, faceted by ancestry
# This is the Filter 1 diagnostic: shows where the percentile floor sits
# relative to the score distribution of assigned samples.
# Requires loglik_best_score column produced by the updated classifier.
# ---------------------------------------------------------------------------

if ("loglik_best_score" %in% names(classified) && nrow(assigned) > 0) {

  # Compute per-ancestry score floors from assigned samples as a proxy
  # (in practice these come from the reference panel inside the classifier;
  #  here we approximate by showing the percentile of the assigned scores)
  floor_lines <- assigned[, .(
    floor_approx = quantile(loglik_best_score, probs = percentile_floor, na.rm = TRUE)
  ), by = closest_ancestry]

  p3 <- ggplot(assigned, aes(x = loglik_best_score, fill = closest_ancestry)) +
    geom_histogram(binwidth = 0.5, colour = "white", linewidth = 0.2) +
    geom_vline(
      data = floor_lines,
      aes(xintercept = floor_approx),
      linetype = "dashed", colour = "#E24B4A", linewidth = 0.7
    ) +
    facet_wrap(~closest_ancestry, scales = "free", ncol = 3) +
    scale_fill_manual(values = anc_colours, guide = "none") +
    labs(
      title    = "Best log-likelihood score by ancestry (Filter 1 diagnostic)",
      subtitle = paste0(
        "Red dashed line = approx. p", percentile_floor * 100,
        " floor  |  Samples left of line would be rejected by Filter 1"
      ),
      x = "loglik_best_score  (log-likelihood under best-matching ancestry)",
      y = "Count"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(paste0(output, "_score_floor.png"), plot = p3, width = 9, height = 5, dpi = 150)
  cat("Saved: score_floor.png\n")

} else {
  cat("Skipping score_floor.png: loglik_best_score column not found.\n")
}

# ---------------------------------------------------------------------------
# Plot 4 — Rejection summary stacked bar per ancestry
# Shows the proportion assigned vs rejected (and by which filter) per group.
# ---------------------------------------------------------------------------

# Include ambiguous samples, attributing them back to their "best" ancestry
# via loglik_best_score isn't available, so we use rejection_reason directly
all_with_reason <- data[!is.na(rejection_reason) | closest_ancestry %in% ancestries]

if (nrow(all_with_reason) > 0) {

  # For assigned samples, group by closest_ancestry; for rejected, we don't
  # know the intended ancestry — group them as "ambiguous/unassigned"
  all_with_reason[, anc_group := fcase(
    closest_ancestry %in% ancestries, closest_ancestry,
    default = "ambiguous"
  )]

  all_with_reason[, status2 := fcase(
    closest_ancestry %in% ancestries,                    "assigned",
    rejection_reason == "below_percentile_floor",        "below_percentile_floor",
    rejection_reason == "low_delta_margin",              "low_delta_margin",
    rejection_reason == "no_valid_ancestry",             "no_valid_ancestry",
    default = "assigned"
  )]

  summary_counts <- all_with_reason[, .N, by = .(anc_group, status2)]

  p4 <- ggplot(summary_counts, aes(x = anc_group, y = N, fill = status2)) +
    geom_col(position = "stack", colour = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = reason_colours,
      labels = c(
        assigned               = "Assigned",
        below_percentile_floor = "Rejected: below score floor (F1)",
        low_delta_margin       = "Rejected: low delta margin (F2)",
        no_valid_ancestry      = "Rejected: no valid ancestry"
      ),
      name = NULL
    ) +
    labs(
      title    = "Assignment outcome by ancestry group",
      subtitle = "F1 = percentile floor filter  |  F2 = delta margin filter",
      x        = NULL,
      y        = "Sample count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      legend.position = "bottom",
      legend.text     = element_text(size = 9)
    )

  ggsave(paste0(output, "_rejection_summary.png"), plot = p4, width = 7, height = 4.5, dpi = 150)
  cat("Saved: rejection_summary.png\n")
}

# ---------------------------------------------------------------------------
# Console summary table
# ---------------------------------------------------------------------------

summary_tbl <- classified[, .(
  n                  = .N,
  n_assigned         = sum(closest_ancestry %in% ancestries),
  n_floor_rejected   = sum(rejection_reason == "below_percentile_floor", na.rm = TRUE),
  n_delta_rejected   = sum(rejection_reason == "low_delta_margin",        na.rm = TRUE),
  median_margin      = round(median(loglik_margin, na.rm = TRUE), 2)
)]

cat("\n--- Classification summary ---\n")
print(summary_tbl)

per_anc <- classified[closest_ancestry %in% ancestries, .(
  n             = .N,
  median_margin = round(median(loglik_margin, na.rm = TRUE), 2),
  median_score  = if ("loglik_best_score" %in% names(classified))
                    round(median(loglik_best_score, na.rm = TRUE), 2)
                  else NA_real_
), by = closest_ancestry][order(-n)]

cat("\n--- Per-ancestry breakdown (assigned only) ---\n")
print(per_anc)
