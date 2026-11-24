#!/usr/bin/env Rscript
# Plot posterior rate distributions per chain, and check ML / median rates for extreme values.
#
# Usage:
#   Rscript plot_rate_posteriors.R <results folder> [key regex] [nchains]
#   e.g. Rscript scripts/plot_rate_posteriors.R results/dwelltime-model-enhancer-12h '3302'
#
# Reads, for every key in the folder (whose name matches the regex, if given):
#   rates_<key>.txt             header + rows: 1 = ML, 2 = mean, 3 = median, 4 = last sample
#   param-stats_<key>.txt       first line = labels of the FITTED rates
#   sampled_rates_<key>.txt     posterior samples (one row per sample, fitted rates only)
#   ll_sampled_rates_<key>.txt  log-likelihood of each saved sample
# The two sampled_* files are written when the fit ran with writesamples = true.
#
# Chains: the saved samples are an evenly spaced subset of all chains concatenated in order,
# so sample i is assigned to chain ceiling(i * nchains / n). This is exact when the chains stored
# the same number of samples and approximate otherwise; the log-likelihood trace shows the
# chain blocks either way.
#
# Writes to <results folder>/rate_posteriors/:
#   <key>.pdf          page 1: log10 posterior density of each fitted rate, one line per chain,
#                              with the ML (black dashed) and median (red) rates
#                      page 2: log-likelihood of the saved samples by chain, ML level dashed
#   rate_checks.csv    per key and rate: ML, median, posterior 2.5-97.5%, flags
#   chain_summary.csv  per key and chain: best and median log-likelihood, median of each rate
#
# Flags in rate_checks.csv:
#   extreme          ML or median rate below SG_RATE_MIN (1e-4 /min) or above SG_RATE_MAX (100 /min)
#   ml_outside_95    ML rate outside the central 95% of the posterior samples
#   ml_far_median    ML and median differ by more than 10-fold
# A good fit has no flags. Many ml_outside_95 / ml_far_median flags mean the best point is not
# where the sampler spends its time (a narrow spike, or chains in different modes).

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript plot_rate_posteriors.R <results folder> [key regex] [nchains]")
folder <- args[1]
pattern <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "."
nchains <- if (length(args) >= 3) as.integer(args[3]) else 16L
rate_min <- as.numeric(Sys.getenv("SG_RATE_MIN", "1e-4"))
rate_max <- as.numeric(Sys.getenv("SG_RATE_MAX", "100"))

outdir <- file.path(folder, "rate_posteriors")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

rate_files <- list.files(folder, pattern = "^rates_.*\\.txt$")
keys <- sub("^rates_(.*)\\.txt$", "\\1", rate_files)
keys <- keys[grepl(pattern, keys)]
if (length(keys) == 0) stop("no rates_<key>.txt files matching '", pattern, "' in ", folder)

read_header <- function(f) strsplit(readLines(f, n = 1), ",")[[1]]

checks <- list(); chains_out <- list(); skipped <- character()
for (key in keys) {
  f_rates <- file.path(folder, paste0("rates_", key, ".txt"))
  f_stats <- file.path(folder, paste0("param-stats_", key, ".txt"))
  f_samp <- file.path(folder, paste0("sampled_rates_", key, ".txt"))
  f_ll <- file.path(folder, paste0("ll_sampled_rates_", key, ".txt"))
  if (!all(file.exists(c(f_rates, f_stats, f_samp, f_ll)))) { skipped <- c(skipped, key); next }

  rates <- read.csv(f_rates, check.names = FALSE)
  fitted_labels <- read_header(f_stats)
  ml <- unlist(rates[1, fitted_labels]); med <- unlist(rates[3, fitted_labels])
  samp <- as.matrix(read.table(f_samp, header = FALSE))
  if (ncol(samp) != length(fitted_labels)) {
    warning(key, ": ", ncol(samp), " sampled columns but ", length(fitted_labels), " fitted labels; skipped")
    skipped <- c(skipped, key); next
  }
  colnames(samp) <- fitted_labels
  ll <- scan(f_ll, quiet = TRUE)
  n <- nrow(samp)
  chain <- pmin(nchains, ceiling(seq_len(n) * nchains / n))
  ml_ll <- max(ll)

  # ---- checks ----
  q <- apply(samp, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  for (lab in fitted_labels) {
    flags <- c(
      if (min(ml[lab], med[lab]) < rate_min || max(ml[lab], med[lab]) > rate_max) "extreme",
      if (ml[lab] < q[1, lab] || ml[lab] > q[2, lab]) "ml_outside_95",
      if (max(ml[lab], med[lab]) / max(min(ml[lab], med[lab]), 1e-300) > 10) "ml_far_median")
    checks[[length(checks) + 1]] <- data.frame(
      key = key, rate = lab, ML = ml[lab], median = med[lab],
      q2.5 = q[1, lab], q97.5 = q[2, lab], flags = paste(flags, collapse = ";"),
      row.names = NULL)
  }
  per_chain <- do.call(rbind, lapply(sort(unique(chain)), function(c) {
    i <- chain == c
    data.frame(key = key, chain = c, n_samples = sum(i), best_ll = max(ll[i]), median_ll = median(ll[i]),
               t(apply(samp[i, , drop = FALSE], 2, median)), check.names = FALSE, row.names = NULL)
  }))
  chains_out[[key]] <- per_chain

  # ---- plots ----
  long <- data.frame(rate = rep(fitted_labels, each = n), value = as.vector(samp),
                     chain = factor(rep(chain, times = length(fitted_labels))))
  long <- long[is.finite(long$value) & long$value > 0, ]
  marks <- data.frame(rate = rep(fitted_labels, 2), value = c(ml, med),
                      what = rep(c("ML", "median"), each = length(fitted_labels)))
  p1 <- ggplot(long, aes(x = log10(value), colour = chain)) +
    geom_density(linewidth = 0.4, adjust = 1.2) +
    geom_vline(data = marks, aes(xintercept = log10(value), linetype = what),
               colour = c(rep("black", length(fitted_labels)), rep("red", length(fitted_labels)))) +
    scale_linetype_manual(values = c(ML = "dashed", median = "solid")) +
    facet_wrap(~rate, scales = "free") +
    labs(title = key, subtitle = "posterior of log10(rate /min) per chain; black dashed = ML, red = median",
         x = "log10 rate", y = "density", linetype = NULL) +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  trace <- data.frame(i = seq_len(n), ll = ll, chain = factor(chain))
  spread <- max(per_chain$best_ll) - min(per_chain$best_ll)
  p2 <- ggplot(trace, aes(i, ll, colour = chain)) + geom_point(size = 0.3) +
    geom_hline(yintercept = ml_ll, linetype = "dashed") +
    labs(title = key, subtitle = sprintf(
      "log-likelihood of saved samples by chain (dashed = ML %.1f); spread of chain best logL = %.1f",
      ml_ll, spread), x = "saved sample (chains concatenated)", y = "log-likelihood") +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  pdf(file.path(outdir, paste0(key, ".pdf")), width = 11, height = 8)
  print(p1); print(p2)
  invisible(dev.off())
}

if (length(checks) == 0) stop("no keys with sampled_rates_/ll_sampled_rates_ files (run with writesamples=true)")
checks <- do.call(rbind, checks)
write.csv(checks, file.path(outdir, "rate_checks.csv"), row.names = FALSE)
chains_all <- do.call(rbind, lapply(chains_out, function(d) d[, c("key", "chain", "n_samples", "best_ll", "median_ll")]))
write.csv(chains_all, file.path(outdir, "chain_summary.csv"), row.names = FALSE)
for (k in names(chains_out))
  write.csv(chains_out[[k]], file.path(outdir, paste0("chains_", k, ".csv")), row.names = FALSE)

cat("keys plotted:", length(chains_out), " skipped (no samples):", length(skipped), "\n")
cat("output:", outdir, "\n\n")
flagged <- checks[nzchar(checks$flags), ]
if (nrow(flagged) == 0) {
  cat("no flagged rates\n")
} else {
  cat("flagged rates:\n")
  print(flagged[, c("key", "rate", "ML", "median", "q2.5", "q97.5", "flags")], row.names = FALSE, digits = 3)
}
cat("\nspread of best logL across chains (>> 2 means some chains are stuck in worse optima):\n")
spr <- aggregate(best_ll ~ key, chains_all, function(x) max(x) - min(x))
names(spr)[2] <- "best_ll_spread"
print(spr[order(-spr$best_ll_spread), ], row.names = FALSE, digits = 4)
