library(data.table)

# ---------------------------------------------------------------------------
# Ancestry classification via Gaussian log-likelihood (Mahalanobis distance)
# Usage: Rscript ancestry_classifier.R <input.tsv> <output.tsv> \
#          <EUR_delta> <AFR_delta> <EAS_delta> <AMR_delta> <SAS_delta>
#
# Example: Rscript ancestry_classifier.R input.tsv output.tsv 2 2 2 5 4
#
# One delta threshold per ancestry. A sample is labelled "ambiguous" if the
# gap between its best and second-best log-likelihood falls below the
# threshold of its best-matching ancestry.
# Set a threshold to 0 to effectively disable rejection for that ancestry.
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop(
  "Usage: script.R <input> <output> <EUR_delta> <AFR_delta> <EAS_delta> <AMR_delta> <SAS_delta>"
)

data   <- fread(args[1])
output <- args[2]

ancestries <- c("EUR", "AFR", "EAS", "AMR", "SAS")
pc_cols    <- paste0("PC", 1:4)
k          <- length(pc_cols)

# ---------------------------------------------------------------------------
# Per-ancestry delta thresholds
# ---------------------------------------------------------------------------

delta_thresholds <- setNames(as.numeric(args[3:7]), ancestries)

message("Delta thresholds:")
for (anc in ancestries) message(sprintf("  %s: %.2f", anc, delta_thresholds[anc]))

# ---------------------------------------------------------------------------
# Early exit if nothing to classify
# ---------------------------------------------------------------------------

unknown_idx <- which(!data$SUPERPOP %in% ancestries)

if (length(unknown_idx) == 0) {
  message("No unknown samples found — writing input unchanged.")
  fwrite(data, output, sep = "\t", quote = FALSE)
  quit(status = 0)
}

unknown   <- data[unknown_idx, ..pc_cols]
n_unknown <- nrow(unknown)

# Storage — pre-fill with -Inf so skipped ancestries never win
loglik_mat <- matrix(-Inf, nrow = n_unknown, ncol = length(ancestries))
colnames(loglik_mat) <- ancestries

# ---------------------------------------------------------------------------
# Fit one Gaussian per reference ancestry; score unknown samples
# ---------------------------------------------------------------------------

for (anc in ancestries) {

  ref <- data[SUPERPOP == anc, ..pc_cols]

  if (nrow(ref) < 50) {
    warning(sprintf("Skipping '%s': only %d reference samples (need >= 50).", anc, nrow(ref)))
    next
  }

  mu    <- colMeans(ref)
  ridge <- mean(diag(cov(ref))) * 1e-4
  Sigma <- cov(ref) + diag(ridge, k)

  Sigma_inv <- solve(Sigma)
  logdet    <- as.numeric(determinant(Sigma, logarithm = TRUE)$modulus)

  message(sprintf("  %s: n_ref=%d", anc, nrow(ref)))

  X  <- as.matrix(unknown)
  Xc <- sweep(X, 2, mu, "-")
  md <- rowSums((Xc %*% Sigma_inv) * Xc)

  loglik_mat[, anc] <- -0.5 * (md + logdet + k * log(2 * pi))
}

# ---------------------------------------------------------------------------
# Assign best ancestry
# ---------------------------------------------------------------------------

all_missing <- rowSums(is.finite(loglik_mat)) == 0

best_anc <- ifelse(
  all_missing,
  "unassigned",
  colnames(loglik_mat)[max.col(loglik_mat, ties.method = "first")]
)

# Best log-likelihood score per sample
best_score <- apply(loglik_mat, 1, function(r) {
  finite_r <- r[is.finite(r)]
  if (length(finite_r) == 0) NA_real_ else max(finite_r)
})

# ---------------------------------------------------------------------------
# Filter — per-ancestry delta margin
# Reject if best and second-best log-likelihoods are too close, using the
# threshold of the best-matching ancestry
# ---------------------------------------------------------------------------

n_valid_anc <- sum(colSums(is.finite(loglik_mat)) > 0)

if (n_valid_anc < 2) {
  warning("Fewer than 2 ancestries have sufficient reference samples — skipping delta filter.")
  delta <- rep(NA_real_, n_unknown)
} else {
  sorted <- t(apply(loglik_mat, 1, function(r) {
    s <- sort(r[is.finite(r)], decreasing = TRUE)
    if (length(s) < 2) c(s[1], -Inf) else s[1:2]
  }))
  delta <- sorted[, 1] - sorted[, 2]

  # Apply per-ancestry threshold
  delta_fail <- mapply(function(anc, d) {
    if (anc %in% c("unassigned", "ambiguous")) return(FALSE)
    thr <- delta_thresholds[anc]
    if (is.na(thr)) return(FALSE)
    d < thr
  }, best_anc, delta)

  best_anc[!all_missing & delta_fail] <- "ambiguous"
}

# ---------------------------------------------------------------------------
# Rejection reason
# ---------------------------------------------------------------------------

rejection_reason <- rep(NA_character_, n_unknown)
rejection_reason[!all_missing & best_anc == "ambiguous"] <- "low_delta_margin"
rejection_reason[all_missing]                             <- "no_valid_ancestry"

# ---------------------------------------------------------------------------
# Write results
# ---------------------------------------------------------------------------

data[unknown_idx, closest_ancestry  := best_anc]
data[unknown_idx, loglik_margin     := delta]
data[unknown_idx, loglik_best_score := best_score]
data[unknown_idx, rejection_reason  := rejection_reason]

fwrite(data, output, sep = "\t", quote = FALSE)

# Summary to console
n_assigned   <- sum(best_anc %in% ancestries)
n_ambiguous  <- sum(best_anc == "ambiguous")
n_unassigned <- sum(best_anc == "unassigned")

message(sprintf(
  "\nDone. %d samples:  %d assigned  |  %d ambiguous  |  %d unassigned",
  n_unknown, n_assigned, n_ambiguous, n_unassigned
))
message(sprintf(
  "Ambiguous breakdown:  %d failed delta margin  |  %d no valid ancestry",
  sum(rejection_reason == "low_delta_margin",  na.rm = TRUE),
  sum(rejection_reason == "no_valid_ancestry", na.rm = TRUE)
))
message(sprintf("Written to: %s", output))
