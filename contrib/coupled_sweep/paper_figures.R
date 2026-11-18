#!/usr/bin/env Rscript
# Paper figures from a simulation run (after 3_run_simulations.jl), in five steps:
#   1. Shape map: every simulation's cross-correlation centre against its flanks, one panel per model type,
#      coloured by NRMSE; the shaded area = a Mexican-hat shape (selection 1).
#   2. Folds: for every connection (e.g. CCAT1 state 2 -> MYC R23), the fold in each experimental track (cell)
#      2a: against tracks paired with another cell (permuted: the relationship expected by chance);
#      2b: against the folds the simulations in selection 1 use (points).
#   3. Heatmap of every simulation in selection 1, by model type.
#   4. Heatmap of selection 2: selection 1, keeping only simulations whose every CCAT1-MYC fold lies within the
#      middle DATA_RANGE of the experimental per-track folds of that connection (hidden-unit models have no
#      observable fold and are kept).
#   5. Ranking of the architectures in selection 2 by NRMSE (best parameter set of each), and the candidate list.
#   6. Every architecture of selection 2, in ranking order: the experimental cross-correlation and the simulated one
#      (mean +/- SE over replicates) of their best parameter set, analysed with insert steps 11, 13 and 33.
#
#   Rscript paper_figures.R RUN_DIR INSERT EXPERIMENTAL_CSV CCAT1_PREDICTED.csv MYC_PREDICTED.csv OUTPUT_FOLDER
#   e.g. Rscript paper_figures.R runs/semi_extensive_insert_steps 11 data/experimental.csv \
#            data/predictedtraces_..._enhancer.csv data/predictedtraces_..._gene.csv figures_insert11
#
# INSERT: 11, 13 or 33 reads RUN_DIR/results/insert<INSERT>; none (or "" or -) reads RUN_DIR/results: runs made
# without insert_steps, e.g. earlier runs. RUN_DIR: the run_dir of a regular config (simulations.csv), or
# <run_root>/<profile> of an expanded config (expanded_models.csv and expanded_simulations.csv, binary .f32 curves).
# Needs R with ggplot2 (data.table, if installed, makes reading large runs faster).
# Writes PDF and PNG figures and CSV tables into OUTPUT_FOLDER. Large runs: every simulation's centre and flanks are
# computed while reading, but curves are kept only for the hat-shaped ones; the shape map draws at most
# MAX_POINTS points per model type (its counts use all).

suppressPackageStartupMessages(library(ggplot2))

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
INTERVAL_MIN  <- 5 / 3          # minutes per frame
CENTER_MIN    <- 5              # hat centre: mean cross-correlation at lags within +/- this many minutes
FLANKS_MIN    <- c(20, 60)      # hat flanks: mean cross-correlation at lags between these, both sides
CENTER_AT_LEAST <- 0.01         # selection 1: centre >= this ...
FLANKS_AT_MOST  <- -0.01        # ... and flanks <= this (the shaded area)
DATA_RANGE    <- c(0.25, 0.75)  # selection 2: quantiles of the per-track folds a simulated fold must lie between
MIN_FRAMES    <- 5              # per track: frames needed in and out of the condition to compute a fold
N_NULL        <- 1000           # permutations: each track paired with another cell's track, this many times
MIN_TRACKS    <- 5              # a connection's experimental distribution is drawn (and used in selection 2) only
                                # with at least this many tracks that have a fold
RANK_BY       <- "xcor"         # "xcor" = NRMSE of the cross-correlation; "mean3" = mean NRMSE of it and both autocorrelations
TOP_FRACTION  <- 0.25           # candidates: the best this fraction of architectures in selection 2 ...
NRMSE_CUTOFF  <- NA             # ... or, if a number, every architecture with NRMSE below it
HEATMAP_LAG   <- 120            # heatmap lag range, minutes
MAX_HEATMAP_ROWS <- 20000       # heatmaps: above this many passing simulations, the best one per architecture is drawn
TOP_CURVES    <- Inf            # figure 6: this many best architectures of selection 2 (Inf = all)
PANELS_PER_PAGE <- 24           # figure 6: panels per page (6 rows of 4); one multi-page PDF, one PNG per page
INSERT_COLOURS <- c(`11` = "#2a78d6", `13` = "#eb6834", `33` = "#1baf7a")   # figure 6: insert-step variants
MAX_POINTS    <- 100000         # figure 1: at most this many points drawn per model type
SEED          <- 20261118

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) stop("usage: Rscript paper_figures.R RUN_DIR INSERT EXPERIMENTAL_CSV CCAT1_PREDICTED MYC_PREDICTED OUTPUT_FOLDER")
if (!dir.exists(args[1]) && dir.exists(args[2])) args[1:2] <- args[2:1]   # accept INSERT RUN_DIR as well
if (!dir.exists(args[1])) stop("run folder not found: ", args[1])
run_dir <- args[1]; insert <- if (args[2] %in% c("", "none", "-")) "" else args[2]; out <- args[6]
results <- file.path(run_dir, "results", if (nzchar(insert)) paste0("insert", insert) else "")
has_results <- function(folder) length(list.files(folder, pattern = "^metrics_.*\\.csv$")) > 0
if (nzchar(insert) && !has_results(results)) {          # a run made without insert steps
  message("no results in ", results, ": this run has no insert-step variants; reading ", file.path(run_dir, "results"))
  insert <- ""
  results <- file.path(run_dir, "results")
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
GROUP_ORDER <- c("E_to_G", "G_to_E", "reciprocal", "H2", "H3_symmetric", "H3_asymmetric")
GROUP_COLOURS <- c(E_to_G = "#2a78d6", G_to_E = "#eb6834", reciprocal = "#1baf7a", H2 = "#8e5bd0",
                   H3_symmetric = "#d6448f", H3_asymmetric = "#7a7a7a")
RATES <- c("R12", "R21", "R23", "R32", "Initiate")
SOURCES <- c("1", "2", "3", "R", "3+R")
GENE <- c(E = "CCAT1", G = "MYC")
save_plot <- function(plot, name, width, height) {
  ggsave(file.path(out, paste0(name, ".pdf")), plot, width = width, height = height, limitsize = FALSE)
  ggsave(file.path(out, paste0(name, ".png")), plot, width = width, height = height, dpi = 300, limitsize = FALSE)
}

# ---------------------------------------------------------------------------
# Reading (regular and expanded runs; the newest results file wins for a simulation in several)
# ---------------------------------------------------------------------------
# Columns `keep` (a function of the header) of a CSV file, as text or as read.
read_cols <- function(path, keep = function(header) rep(TRUE, length(header)), text = FALSE) {
  header <- gsub('"', "", strsplit(readLines(path, n = 1), ",", fixed = TRUE)[[1]])
  wanted <- header[keep(header)]
  if (requireNamespace("data.table", quietly = TRUE)) {
    d <- data.table::fread(path, select = wanted, colClasses = if (text) "character" else NULL, showProgress = FALSE,
                           na.strings = c("", "NA"))
    return(as.data.frame(d))
  }
  classes <- ifelse(header %in% wanted, if (text) "character" else NA, "NULL")
  read.csv(path, colClasses = classes, check.names = FALSE)
}
# The simulation list with the columns of a regular simulations.csv: sim_id, model_id, group, <edge>_state,
# <edge>_rates (edge EG, GE, HE, HG) and the folds and hidden rates.
read_simulations <- function(run_dir) {
  if (file.exists(file.path(run_dir, "simulations.csv")))
    return(read_cols(file.path(run_dir, "simulations.csv"), text = TRUE))
  models <- read_cols(file.path(run_dir, "expanded_models.csv"), text = TRUE)
  points <- read_cols(file.path(run_dir, "expanded_simulations.csv"), text = TRUE,
                      keep = function(h) h %in% c("sim_id", "architecture_id") | grepl("^fold_|^h[0-9]{2}$", h))
  at <- match(points$architecture_id, models$architecture_id)
  sims <- data.frame(sim_id = points$sim_id, model_id = points$architecture_id, group = models$group[at])
  for (edge in c("EG", "GE", "HE", "HG")) {
    kind <- models[[paste0(edge, "_source_kind")]][at]
    sims[[paste0(edge, "_state")]] <- ifelse(kind == "reporter_count", "R", models[[paste0(edge, "_states")]][at])
    sims[[paste0(edge, "_rates")]] <- models[[paste0(edge, "_rates")]][at]
  }
  cbind(sims, points[, grepl("^fold_|^h[0-9]{2}$", names(points)), drop = FALSE])
}
# The reporter curves of one results file: sim_id, pair, lag_min, mean, se (CSV, or binary .f32 of expanded runs).
read_curve_file <- function(path) {
  if (grepl("\\.f32$", path)) {
    layout <- read.csv(file.path(dirname(path), "curve_layout.csv"))
    bounds <- as.integer(strsplit(sub("^curves_([0-9]+)-([0-9]+)\\.f32$", "\\1 \\2", basename(path)), " ")[[1]])
    n <- bounds[2] - bounds[1] + 1
    values <- readBin(path, "numeric", n = 2 * n * nrow(layout), size = 4)
    rows <- which(layout$observable == "reporter")
    means <- matrix(values[c(TRUE, FALSE)], nrow = nrow(layout))[rows, , drop = FALSE]
    ses <- matrix(values[c(FALSE, TRUE)], nrow = nrow(layout))[rows, , drop = FALSE]
    return(data.frame(sim_id = rep(bounds[1]:bounds[2], each = length(rows)), pair = layout$pair[rows],
                      lag_min = layout$lag_min[rows], mean = as.vector(means), se = as.vector(ses)))
  }
  d <- read_cols(path, keep = function(h) !grepl("^rep", h))
  d[d$observable == "reporter", c("sim_id", "pair", "lag_min", "mean", "se")]
}
curve_files <- function(folder) {                       # newest first
  files <- list.files(folder, pattern = "^curves_[0-9]+-[0-9]+\\.(csv|f32)$", full.names = TRUE)
  files[order(file.mtime(files), decreasing = TRUE)]
}
hat <- function(lag, value) c(center = mean(value[abs(lag) <= CENTER_MIN + 1e-9], na.rm = TRUE),
                              flanks = mean(value[abs(lag) >= FLANKS_MIN[1] & abs(lag) <= FLANKS_MIN[2]], na.rm = TRUE))
# Centre and flanks of every simulation, per results file; keeps the cross-correlation only of hat-shaped ones.
read_shapes <- function(folder) {
  files <- curve_files(folder)
  if (!length(files)) stop("no curves files in ", folder)
  seen <- integer(0); shapes <- list(); kept <- list()
  for (f in files) {
    d <- read_curve_file(f)
    d <- d[d$pair == "EG" & !(d$sim_id %in% seen), ]
    if (!nrow(d)) next
    d$lag_min <- round(d$lag_min, 2)
    ok <- is.finite(d$mean)
    in_center <- ok & abs(d$lag_min) <= CENTER_MIN + 1e-9
    in_flanks <- ok & abs(d$lag_min) >= FLANKS_MIN[1] & abs(d$lag_min) <= FLANKS_MIN[2]
    by <- factor(d$sim_id)
    shape <- data.frame(sim_id = as.integer(levels(by)),
                        center = tapply(ifelse(in_center, d$mean, 0), by, sum) / tapply(in_center, by, sum),
                        flanks = tapply(ifelse(in_flanks, d$mean, 0), by, sum) / tapply(in_flanks, by, sum))
    shape$hat_shape <- shape$center >= CENTER_AT_LEAST & shape$flanks <= FLANKS_AT_MOST
    shapes[[length(shapes) + 1]] <- shape
    kept[[length(kept) + 1]] <- d[d$sim_id %in% shape$sim_id[which(shape$hat_shape)], c("sim_id", "lag_min", "mean")]
    seen <- c(seen, shape$sim_id)
  }
  list(shape = do.call(rbind, shapes), eg = do.call(rbind, kept))
}

sims <- read_simulations(run_dir)
sims$sim_id <- as.integer(sims$sim_id)
metric_files <- list.files(results, pattern = "^metrics_.*\\.csv$", full.names = TRUE)
if (!length(metric_files)) stop("no metrics files in ", results)
metric_files <- metric_files[order(file.mtime(metric_files), decreasing = TRUE)]
wanted <- c("sim_id", "nrmse_reporter_EG", "nrmse_reporter_EE", "nrmse_reporter_GG", "on_E", "on_G")
metrics <- do.call(rbind, lapply(metric_files, read_cols, keep = function(h) h %in% wanted))
metrics <- metrics[!duplicated(metrics$sim_id), ]
read <- read_shapes(results)
eg <- read$eg
experiment <- read.csv(args[3])
experiment <- experiment[experiment$observable == "reporter", ]
experiment$lag_min <- round(experiment$lag_min, 2)

table <- merge(merge(sims, metrics, by = "sim_id"), read$shape, by = "sim_id")
table$nrmse_xcor <- table$nrmse_reporter_EG
table$nrmse_mean3 <- rowMeans(table[, c("nrmse_reporter_EG", "nrmse_reporter_EE", "nrmse_reporter_GG")])
table$rank_value <- if (RANK_BY == "mean3") table$nrmse_mean3 else table$nrmse_xcor
table$hat_shape[is.na(table$hat_shape)] <- FALSE
exp_eg <- experiment[experiment$pair == "EG", ]
exp_hat <- hat(exp_eg$lag_min, exp_eg$mean)
table$group <- factor(table$group, levels = c(GROUP_ORDER, "uncoupled"))
cat(nrow(table), "simulations with results in", results, "\n")

# ---------------------------------------------------------------------------
# 1. Shape map
# ---------------------------------------------------------------------------
coupled <- table[table$group != "uncoupled", ]
coupled$group <- droplevels(coupled$group)
counts <- do.call(rbind, lapply(split(coupled, coupled$group), function(d) data.frame(
  group = d$group[1], label = sprintf("%s\n%d of %d architectures, %s of %s simulations", d$group[1],
    length(unique(d$model_id[d$hat_shape])), length(unique(d$model_id)),
    format(sum(d$hat_shape), big.mark = ","), format(nrow(d), big.mark = ",")))))
labels <- setNames(counts$label, counts$group)
uncoupled <- table[table$group == "uncoupled", ]
coupled <- do.call(rbind, lapply(split(coupled, coupled$group), function(d)
  if (nrow(d) > MAX_POINTS) d[sample.int(nrow(d), MAX_POINTS), ] else d))
coupled <- coupled[order(-coupled$nrmse_xcor), ]           # best drawn last, on top
limits <- function(v, extra) range(c(quantile(v, c(0.01, 0.99), na.rm = TRUE), extra))
p1 <- ggplot(coupled, aes(center, flanks)) +
  annotate("rect", xmin = CENTER_AT_LEAST, xmax = Inf, ymin = -Inf, ymax = FLANKS_AT_MOST, fill = "grey88") +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = nrmse_xcor), size = 0.6, alpha = 0.7, stroke = 0) +
  { if (nrow(uncoupled)) geom_point(data = uncoupled, shape = 23, size = 2.2, fill = "white", colour = "black") } +
  annotate("point", x = exp_hat["center"], y = exp_hat["flanks"], shape = 4, size = 3, stroke = 1.4) +
  facet_wrap(~group, ncol = 3, labeller = labeller(group = labels)) +
  scale_colour_viridis_c("NRMSE of the\ncross-correlation", trans = "log10", direction = -1) +
  coord_cartesian(xlim = limits(coupled$center, c(exp_hat["center"], 0)), ylim = limits(coupled$flanks, c(exp_hat["flanks"], 0))) +
  labs(x = sprintf("centre: mean cross-correlation within +/-%g min", CENTER_MIN),
       y = sprintf("flanks: mean from %g to %g min, both sides", FLANKS_MIN[1], FLANKS_MIN[2]),
       title = sprintf("Cross-correlation shape (insert %s). Shaded: centre >= %g and flanks <= %g. X = experiment, diamond = uncoupled",
                       insert, CENTER_AT_LEAST, FLANKS_AT_MOST)) +
  theme_classic(base_size = 9) + theme(strip.background = element_blank(), plot.title = element_text(size = 9, face = "bold"))
save_plot(p1, "1_shape_map", 11, 7)

# ---------------------------------------------------------------------------
# 2. Folds in the experimental tracks, in shifted / permuted tracks, and in the hat-shaped simulations
# ---------------------------------------------------------------------------
read_tracks <- function(path) {
  d <- read.csv(path, check.names = FALSE, na.strings = c("", "NA", "NaN", "missing"))
  ids <- sort(as.integer(sub("Reporters_", "", grep("^Reporters_[0-9]+$", names(d), value = TRUE))))
  lapply(ids, function(i) {
    reporters <- d[[paste0("Reporters_", i)]]
    n <- sum(!is.na(reporters))
    first <- as.integer(sub("^\\[\\s*([0-9]+).*", "\\1", d[[paste0("Rstate_", i)]][seq_len(n)]))
    list(state = as.integer(d[[paste0("Gstate_", i)]][seq_len(n)]), first = first, reporters = reporters[seq_len(n)])
  })
}
tracks <- list(E = read_tracks(args[4]), G = read_tracks(args[5]))
stopifnot(length(tracks$E) == length(tracks$G))
# Per track: the number of coupling effects of each source condition at frames 1..n-1, as in the model:
# states 1, 2, 3: 1 while in that state; R: the number of reporters; 3+R: 1 in state 3 plus the number of reporters.
conditions <- function(t, n) {
  s <- t$state[seq_len(n - 1)]; r <- t$reporters[seq_len(n - 1)]
  cbind(`1` = s == 1, `2` = s == 2, `3` = s == 3, R = r, `3+R` = (s == 3) + r) * 1
}
rates <- function(t, n) {
  s <- t$state; f <- t$first; now <- seq_len(n - 1); nxt <- now + 1
  risk <- cbind(R12 = s[now] == 1, R21 = s[now] == 2, R23 = s[now] == 2, R32 = s[now] == 3, Initiate = s[now] == 3 & f[now] == 0)
  event <- cbind(R12 = s[nxt] == 2, R21 = s[nxt] == 1, R23 = s[nxt] == 3, R32 = s[nxt] == 2, Initiate = f[nxt] == 1) & risk
  list(risk = risk * 1, event = event * 1)
}
# Fold of one track, as in the model: the rate is base x (1 + (fold - 1) x n), n = number of coupling effects.
# The base comes from the frames with n = 0; the fold makes the expected events in the frames with n > 0 equal the
# observed ones:  fold = 1 + (k1 x T0 / k0 - T1) / N1   (k = events, T = frames at risk, N1 = sum of n over the
# frames with n > 0; 0.5 added to the event counts). For states 1, 2 and 3 (n = 1) this is simply
# (rate while in the state) / (rate otherwise). Folds below 0.001 are held at 0.001 (the model clips too).
# NA without MIN_FRAMES frames on each side or without any event. Returns a 5 (condition) x 5 (rate) matrix.
track_fold <- function(C, R) {
  on <- (C > 0) * 1
  T1 <- crossprod(on, R$risk); k1 <- crossprod(on, R$event); N1 <- crossprod(C, R$risk)
  T0 <- matrix(colSums(R$risk), nrow(T1), ncol(T1), byrow = TRUE) - T1
  k0 <- matrix(colSums(R$event), nrow(T1), ncol(T1), byrow = TRUE) - k1
  fold <- pmax(1 + ((k1 + 0.5) * T0 / (k0 + 0.5) - T1) / N1, 1e-3)
  fold[T1 < MIN_FRAMES | T0 < MIN_FRAMES | (k1 + k0) == 0] <- NA
  fold
}
as_rows <- function(folds, direction, kind) {         # list of 5 x 5 matrices -> long data frame
  values <- do.call(rbind, lapply(folds, as.vector))
  src <- if (direction == "EG") "E" else "G"; tgt <- if (direction == "EG") "G" else "E"
  names <- as.vector(outer(SOURCES, RATES, function(s, r) sprintf("%s %s -> %s %s", GENE[src], s, GENE[tgt], r)))
  data.frame(component = rep(names, each = nrow(values)), fold = as.vector(values), kind = kind)
}
ntracks <- length(tracks$E)
fold_rows <- list()
for (direction in c("EG", "GE")) {
  source <- tracks[[substr(direction, 1, 1)]]; target <- tracks[[substr(direction, 2, 2)]]
  n <- vapply(seq_len(ntracks), function(i) min(length(source[[i]]$state), length(target[[i]]$state)), 1L)
  C <- lapply(seq_len(ntracks), function(i) conditions(source[[i]], n[i]))
  R <- lapply(seq_len(ntracks), function(i) rates(target[[i]], n[i]))
  fold_rows[[length(fold_rows) + 1]] <- as_rows(lapply(seq_len(ntracks), function(i) track_fold(C[[i]], R[[i]])), direction, "experiment")
  permuted <- list()
  for (b in seq_len(N_NULL)) {
    for (i in seq_len(ntracks)) {                     # the source track of another cell
      j <- sample(setdiff(seq_len(ntracks), i), 1)
      len <- min(n[i], n[j])
      Rj <- lapply(R[[i]], function(x) x[seq_len(len - 1), , drop = FALSE])
      permuted[[length(permuted) + 1]] <- track_fold(conditions(source[[j]], len), Rj)
    }
  }
  fold_rows[[length(fold_rows) + 1]] <- as_rows(permuted, direction, "permuted cells")
}
folds_data <- do.call(rbind, fold_rows)
folds_data <- folds_data[is.finite(folds_data$fold) & folds_data$fold > 0, ]
tracks_with_fold <- table(folds_data$component[folds_data$kind == "experiment"])
measured <- names(tracks_with_fold)[tracks_with_fold >= MIN_TRACKS]   # connections with an experimental distribution
data_summary <- do.call(rbind, lapply(split(folds_data, list(folds_data$component, folds_data$kind), drop = TRUE), function(d) {
  q <- quantile(d$fold, c(0.025, DATA_RANGE, 0.5, 0.975))
  data.frame(component = d$component[1], kind = d$kind[1], values = nrow(d), median = q[4], low = q[2], high = q[3],
             q025 = q[1], q975 = q[5])
}))
write.csv(data_summary, file.path(out, "2_data_folds.csv"), row.names = FALSE)
range_of <- with(data_summary[data_summary$kind == "experiment" & data_summary$component %in% measured, ],
                 setNames(Map(c, low, high), component))

# Folds of the hat-shaped simulations, one row per simulation and connection component
components_of <- function(row) {
  rows <- list()
  for (edge in c("EG", "GE")) {
    set <- row[[paste0(edge, "_rates")]]
    if (is.na(set) || !nzchar(set)) next
    src <- substr(edge, 1, 1); tgt <- substr(edge, 2, 2)
    for (rate in strsplit(set, "+", fixed = TRUE)[[1]]) {
      value <- row[[paste0("fold_", edge, "_", rate)]]
      if (is.null(value) || is.na(value) || !nzchar(value)) value <- row[[paste0("fold_", edge)]]
      rows[[length(rows) + 1]] <- data.frame(sim_id = row$sim_id, model_id = row$model_id, group = as.character(row$group),
        component = sprintf("%s %s -> %s %s", GENE[src], row[[paste0(edge, "_state")]], GENE[tgt], rate), fold = as.numeric(value))
    }
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}
# How many simulations (of the whole run) use each connection component
count_components <- function(tab) {
  counts <- setNames(numeric(0), character(0))
  for (edge in c("EG", "GE")) {
    state <- tab[[paste0(edge, "_state")]]; set <- tab[[paste0(edge, "_rates")]]
    ok <- !is.na(set) & nzchar(set) & tab$group != "uncoupled" & !grepl("gamma0", tab$model_id)
    combos <- table(paste(state[ok], set[ok], sep = "|"))
    for (k in names(combos)) {
      parts <- strsplit(k, "|", fixed = TRUE)[[1]]
      for (rate in strsplit(parts[2], "+", fixed = TRUE)[[1]]) {
        name <- sprintf("%s %s -> %s %s", GENE[substr(edge, 1, 1)], parts[1], GENE[substr(edge, 2, 2)], rate)
        counts[name] <- sum(counts[name], combos[[k]], na.rm = TRUE)
      }
    }
  }
  counts
}
selected1 <- table[table$hat_shape, ]
sim_folds <- do.call(rbind, lapply(seq_len(nrow(selected1)), function(i) components_of(selected1[i, ])))
if (is.null(sim_folds)) sim_folds <- data.frame(sim_id = character(), model_id = character(), group = character(), component = character(), fold = numeric())
write.csv(sim_folds, file.path(out, "2_simulation_folds.csv"), row.names = FALSE)

order_components <- as.vector(t(outer(c("CCAT1", "MYC"), as.vector(outer(SOURCES, RATES, paste)), function(g, sr) {
  parts <- strsplit(sr, " ", fixed = TRUE)
  vapply(seq_along(g), function(k) sprintf("%s %s -> %s %s", g[k], parts[[k]][1], ifelse(g[k] == "CCAT1", "MYC", "CCAT1"), parts[[k]][2]), "")
})))
used <- count_components(table)
hat_used <- table(sim_folds$component)
label <- function(component) {
  tracks <- if (component %in% names(tracks_with_fold)) tracks_with_fold[[component]] else 0
  sprintf("%s   (tracks %d; hat-shaped %s of %s)", component, tracks,
          format(if (component %in% names(hat_used)) hat_used[[component]] else 0, big.mark = ","),
          format(if (component %in% names(used)) used[[component]] else 0, big.mark = ","))
}
violins <- function(d, rows) {                         # at most 20,000 values per violin, rows as labelled factor
  d <- d[d$component %in% rows, ]
  d <- do.call(rbind, lapply(split(d, list(d$component, d$kind), drop = TRUE),
                             function(x) if (nrow(x) > 20000) x[sample.int(nrow(x), 20000), ] else x))
  d$component <- factor(d$component, levels = rev(rows), labels = vapply(rev(rows), label, ""))
  d
}
fold_axis <- list(geom_hline(yintercept = 1, colour = "grey50", linewidth = 0.3), scale_y_log10(), coord_flip(),
                  theme_classic(base_size = 8), theme(plot.title = element_text(size = 8.5, face = "bold")))

# 2a. Experiment against permuted cells, for the connections with an experimental distribution
rows_a <- intersect(order_components, measured)
data_a <- violins(folds_data, rows_a)
data_a$kind <- factor(data_a$kind, levels = c("permuted cells", "experiment"))
medians_a <- aggregate(fold ~ component + kind, data_a, median)
p2a <- ggplot(data_a, aes(component, fold)) + fold_axis +
  geom_violin(aes(colour = kind), fill = NA, position = position_dodge(width = 0.8), scale = "width", width = 0.75, linewidth = 0.35) +
  geom_point(data = medians_a, aes(colour = kind), position = position_dodge(width = 0.8), shape = 124, size = 3) +
  scale_colour_manual(NULL, values = c(experiment = "#2a78d6", `permuted cells` = "grey55"), breaks = c("experiment", "permuted cells")) +
  labs(x = NULL, y = "fold: rate while the source condition holds / rate otherwise (log scale)",
       title = sprintf("Folds per experimental track (blue) and with each track paired with another cell (grey, %d permutations); | = median",
                       N_NULL))
save_plot(p2a, "2a_folds_data_vs_permuted", 10, max(3, 0.3 * length(rows_a) + 1.3))

# 2b. Experiment against the folds of the hat-shaped simulations
rows_b <- intersect(order_components, unique(c(measured, as.character(sim_folds$component))))
data_b <- violins(folds_data[folds_data$kind == "experiment", ], intersect(rows_b, measured))
data_b$component <- factor(as.character(data_b$component), levels = vapply(rev(rows_b), label, ""))
points_b <- sim_folds[sim_folds$component %in% rows_b, ]
points_b$component <- factor(vapply(as.character(points_b$component), label, ""), levels = vapply(rev(rows_b), label, ""))
medians_b <- if (nrow(data_b)) aggregate(fold ~ component, data_b, median) else data_b
p2b <- ggplot(data_b, aes(component, fold)) + fold_axis +
  geom_violin(colour = "#2a78d6", fill = NA, scale = "width", width = 0.8, linewidth = 0.4) +
  geom_point(data = medians_b, colour = "#2a78d6", shape = 124, size = 3.5) +
  geom_jitter(data = points_b, aes(fill = group), shape = 21, colour = "white", stroke = 0.15, size = 1.8, alpha = 0.8,
              width = 0.25, height = 0) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual("hat-shaped\nsimulations", values = GROUP_COLOURS, drop = TRUE) +
  labs(x = NULL, y = "fold: rate while the source condition holds / rate otherwise (log scale)",
       title = "Folds per experimental track (blue violin, | = median) and the folds the hat-shaped simulations use (points)")
save_plot(p2b, "2b_folds_data_vs_simulations", 10, max(3, 0.3 * length(rows_b) + 1.3))

# ---------------------------------------------------------------------------
# 3 and 4. Heatmaps of selection 1 and selection 2
# ---------------------------------------------------------------------------
in_range <- function(component, fold) {
  r <- range_of[[as.character(component)]]
  is.null(r) || (fold >= r[1] & fold <= r[2])          # no experimental estimate: not constrained
}
sim_folds$in_range <- mapply(in_range, sim_folds$component, sim_folds$fold)
outside <- unique(sim_folds$sim_id[!sim_folds$in_range])
selected2 <- selected1[!(selected1$sim_id %in% outside), ]
selected1$selection2 <- !(selected1$sim_id %in% outside)
write.csv(selected1, file.path(out, "3_selection1_hat_shape.csv"), row.names = FALSE)
write.csv(selected2, file.path(out, "4_selection2_folds_in_data_range.csv"), row.names = FALSE)

heatmap <- function(selected, title) {
  if (!nrow(selected)) return(ggplot() + annotate("text", x = 0, y = 0, label = "no simulations") + theme_void())
  if (nrow(selected) > MAX_HEATMAP_ROWS) {               # too many rows to draw: the best per architecture
    selected <- selected[order(selected$rank_value), ]
    selected <- selected[!duplicated(selected$model_id), ]
    title <- paste0(title, sprintf(" (more than %s: best parameter set per architecture shown)", format(MAX_HEATMAP_ROWS, big.mark = ",")))
  }
  selected <- selected[order(selected$group, selected$rank_value), ]
  selected$row <- ave(seq_len(nrow(selected)), selected$group, FUN = seq_along)
  d <- merge(eg[eg$sim_id %in% selected$sim_id & abs(eg$lag_min) <= HEATMAP_LAG, ], selected[, c("sim_id", "group", "row")])
  e <- exp_eg[abs(exp_eg$lag_min) <= HEATMAP_LAG, ]
  e <- data.frame(sim_id = "experiment", lag_min = e$lag_min, mean = e$mean, group = "experiment", row = 1)
  d <- rbind(e, d[, names(e)])
  n_rows <- table(selected$group)
  d$group <- factor(d$group, levels = c("experiment", names(n_rows)[n_rows > 0]),
                    labels = c("experiment", sprintf("%s\n%d", names(n_rows)[n_rows > 0], n_rows[n_rows > 0])))
  limit <- 2 * max(abs(e$mean), na.rm = TRUE)
  lags <- sort(unique(d$lag_min))                        # multi-tau lags are uneven: each tile spans to the midpoints
  edges <- c(lags[1] - (lags[2] - lags[1]) / 2, (lags[-1] + lags[-length(lags)]) / 2,
             lags[length(lags)] + (lags[length(lags)] - lags[length(lags) - 1]) / 2)
  d$xmin <- edges[match(d$lag_min, lags)]; d$xmax <- edges[match(d$lag_min, lags) + 1]
  d$ymin <- d$row - 0.5; d$ymax <- d$row + 0.5
  strip <- d$group == "experiment"                       # the experiment strip: 4% of the rows, at least 3
  d$ymax[strip] <- 0.5 + max(3, round(0.04 * nrow(selected)))
  ggplot(d, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = pmax(pmin(mean, limit), -limit))) +
    geom_rect() +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2("cross-\ncorrelation", low = "#104281", mid = "#f0efec", high = "#8f1f1f", limits = c(-limit, limit)) +
    scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
    labs(x = "lag (min); positive = MYC after CCAT1", y = NULL, title = title) +
    theme_minimal(base_size = 8) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(), strip.text.y.left = element_text(angle = 0, hjust = 1),
          panel.spacing.y = unit(2, "pt"), plot.title = element_text(size = 9, face = "bold"))
}
height <- function(n) min(48, max(4, 2 + 0.012 * n))
save_plot(heatmap(selected1, sprintf("Selection 1: %s hat-shaped simulations, best NRMSE on top within each type", format(nrow(selected1), big.mark = ","))),
          "3_heatmap_selection1", 8, height(nrow(selected1)))
save_plot(heatmap(selected2, sprintf("Selection 2: %s of them with every CCAT1-MYC fold within the middle %g%% of the experimental tracks",
                                     format(nrow(selected2), big.mark = ","), 100 * diff(DATA_RANGE))),
          "4_heatmap_selection2", 8, height(nrow(selected2)))

# ---------------------------------------------------------------------------
# 5. Ranking of the architectures in selection 2 and the candidate list
# ---------------------------------------------------------------------------
if (nrow(selected2)) {
  best <- selected2[order(selected2$rank_value), ]
  best <- best[!duplicated(best$model_id), ]
  best$rank <- seq_len(nrow(best))
  best$candidate <- if (is.na(NRMSE_CUTOFF)) best$rank <= ceiling(TOP_FRACTION * nrow(best)) else best$rank_value < NRMSE_CUTOFF
  keep_cols <- c("rank", "candidate", "model_id", "group", "sim_id", "rank_value", "nrmse_reporter_EG", "nrmse_reporter_EE",
                 "nrmse_reporter_GG", "center", "flanks", "on_E", "on_G",
                 grep("^(EG|GE|HE|HG)_(state|rates)$|^fold_|^h[0-9]{2}$", names(best), value = TRUE))
  ranked <- best[, keep_cols]
  ranked <- ranked[, colSums(!is.na(ranked) & ranked != "") > 0]
  write.csv(ranked, file.path(out, "5_ranked_architectures.csv"), row.names = FALSE)
  write.csv(ranked[ranked$candidate, ], file.path(out, "5_candidates_to_fit.csv"), row.names = FALSE)
  best$model_id <- factor(best$model_id, levels = rev(best$model_id))
  p5 <- ggplot(best, aes(rank_value, model_id, colour = group, shape = candidate)) +
    geom_point(size = 1.6) +
    scale_colour_manual(NULL, values = GROUP_COLOURS) + scale_shape_manual("candidate", values = c(`TRUE` = 16, `FALSE` = 1)) +
    labs(x = if (RANK_BY == "mean3") "mean NRMSE (cross-correlation and both autocorrelations)" else "NRMSE of the cross-correlation",
         y = NULL, title = sprintf("Selection 2 architectures (best parameter set each), best first; candidates = %s",
                                   if (is.na(NRMSE_CUTOFF)) sprintf("top %g%%", 100 * TOP_FRACTION) else sprintf("NRMSE < %g", NRMSE_CUTOFF))) +
    theme_classic(base_size = 7) + theme(plot.title = element_text(size = 8.5, face = "bold"))
  save_plot(p5, "5_ranking", 9, max(3, 0.13 * nrow(best) + 1.2))

  # 6. Top architectures: experiment and simulation, with the reporter at insert steps 11, 13 and 33
  top <- head(best[order(best$rank), ], TOP_CURVES)
  ids <- as.integer(top$sim_id)
  read_top <- function(folder) {                       # mean and SE of the cross-correlation of the top simulations
    files <- curve_files(folder)
    if (!length(files)) return(NULL)
    bounds <- do.call(rbind, lapply(strsplit(sub("^curves_([0-9]+)-([0-9]+)\\..*$", "\\1 \\2", basename(files)), " "), as.integer))
    files <- files[vapply(seq_along(files), function(i) any(ids >= bounds[i, 1] & ids <= bounds[i, 2]), TRUE)]
    if (!length(files)) return(NULL)
    d <- do.call(rbind, lapply(files, function(f) {
      part <- read_curve_file(f)
      part[part$pair == "EG" & part$sim_id %in% ids, c("sim_id", "lag_min", "mean", "se")]
    }))
    d$lag_min <- round(d$lag_min, 2)
    d[!duplicated(d[, c("sim_id", "lag_min")]), ]
  }
  variants <- if (nzchar(insert)) c("11", "13", "33") else ""
  variants <- variants[vapply(file.path(run_dir, "results", ifelse(nzchar(variants), paste0("insert", variants), "")), has_results, TRUE)]
  sim_curves <- do.call(rbind, lapply(variants, function(v) {
    d <- read_top(file.path(run_dir, "results", if (nzchar(v)) paste0("insert", v) else ""))
    if (is.null(d)) NULL else cbind(d, insert = if (nzchar(v)) v else "all positions")
  }))
  panel <- setNames(sprintf("%d. %s", top$rank, gsub("__", " ", top$model_id)), top$sim_id)
  sim_curves <- sim_curves[abs(sim_curves$lag_min) <= HEATMAP_LAG, ]
  sim_curves$panel <- factor(panel[as.character(sim_curves$sim_id)], levels = panel)
  e <- exp_eg[abs(exp_eg$lag_min) <= HEATMAP_LAG, ]
  e <- do.call(rbind, lapply(panel, function(p) cbind(e[, c("lag_min", "mean", "se")], panel = p)))
  e$panel <- factor(e$panel, levels = panel)
  colours <- if (nzchar(insert)) INSERT_COLOURS else c(`all positions` = "#2a78d6")
  pages <- split(levels(sim_curves$panel), ceiling(seq_along(levels(sim_curves$panel)) / PANELS_PER_PAGE))
  name <- "6_selection2_architectures_insert_steps"
  pdf(file.path(out, paste0(name, ".pdf")), width = 11, height = 1.9 * ceiling(PANELS_PER_PAGE / 4) + 1)
  for (k in seq_along(pages)) {
    on_page <- pages[[k]]
    p6 <- ggplot(sim_curves[sim_curves$panel %in% on_page, ], aes(lag_min, mean)) +
      geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
      geom_ribbon(data = e[e$panel %in% on_page, ], aes(ymin = mean - se, ymax = mean + se), fill = "grey80") +
      geom_line(data = e[e$panel %in% on_page, ], colour = "black", linewidth = 0.5) +
      geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = insert), alpha = 0.2) +
      geom_line(aes(colour = insert), linewidth = 0.5) +
      facet_wrap(~panel, ncol = 4, labeller = label_wrap_gen(38)) +
      scale_colour_manual("insert step\n(CCAT1, MYC)", values = colours) +
      scale_fill_manual("insert step\n(CCAT1, MYC)", values = colours) +
      labs(x = "lag (min); positive = MYC after CCAT1", y = "cross-correlation (reporter count)",
           title = sprintf("Selection 2: %d architectures in ranking order (insert %s), page %d of %d: experiment (black, +/- SE), simulation (mean +/- SE)",
                           nrow(top), if (nzchar(insert)) insert else "all positions", k, length(pages))) +
      theme_classic(base_size = 7) + theme(strip.background = element_blank(), strip.text = element_text(size = 6),
                                           plot.title = element_text(size = 8.5, face = "bold"))
    print(p6)
    ggsave(file.path(out, sprintf("%s_page%02d.png", name, k)), p6, width = 11,
           height = 1.9 * ceiling(length(on_page) / 4) + 1, dpi = 300)
  }
  invisible(dev.off())
}
cat(sprintf("Selection 1 (hat shape): %d simulations, %d architectures\n", nrow(selected1), length(unique(selected1$model_id))))
cat(sprintf("Selection 2 (folds within the data's middle %g%%): %d simulations, %d architectures\n",
            100 * diff(DATA_RANGE), nrow(selected2), length(unique(selected2$model_id))))
cat("Wrote figures and tables to", out, "\n")
