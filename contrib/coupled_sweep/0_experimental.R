#!/usr/bin/env Rscript
# Step 0. Compute the experimental curves and ON fractions that the simulations are compared with.
#
#   Rscript 0_experimental.R CORRELATION.R CCAT1_PREDICTED.csv MYC_PREDICTED.csv data 200 1e-7 1.6666666666666667
#
# Arguments: your fcs_multitau_xcor_final.R, the CCAT1 and MYC predictedtraces_*.csv files, the output folder,
# the number of bootstrap replicates, epsilon, and the frame interval in minutes.
# Writes into the output folder:
#   experimental.csv           the six experimental correlation curves (mean and bootstrap SE)
#   experimental_metrics.csv   experimental ON fractions and mean reporter counts
#   experimental_lengths.csv   number of frames of every experimental trace
#   experimental_reporters_long.csv   the traces in long format
# The code below is the delivered experimental.R, unchanged.
# Invoke the user's actual R functions in a private environment (no CLI side effects).
args <- commandArgs(TRUE)
stopifnot(length(args) >= 5)
source.file <- args[1]; ccat1.file <- args[2]; myc.file <- args[3]; out <- args[4]
B <- as.integer(args[5]); epsilon <- if (length(args)>5) as.numeric(args[6]) else 1e-7
interval.min <- if (length(args)>6) as.numeric(args[7]) else 5/3
env <- new.env(parent=globalenv()); sys.source(source.file, envir=env)
e <- read.csv(ccat1.file, check.names=FALSE, na.strings=c("", "NA", "NaN", "missing"))
g <- read.csv(myc.file, check.names=FALSE, na.strings=c("", "NA", "NaN", "missing"))
ids <- function(d) sort(as.integer(sub("Reporters_", "", grep("^Reporters_[0-9]+$", names(d), value=TRUE))))
stopifnot(identical(ids(e), ids(g)), length(ids(e)) == 137)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
traces <- lapply(ids(e), function(i) {
  x <- e[[paste0("Reporters_",i)]]; y <- g[[paste0("Reporters_",i)]]
  stopifnot(identical(is.na(x),is.na(y)))
  n <- sum(!is.na(x)); stopifnot(n >= 17, all(!is.na(x[seq_len(n)])), all(is.na(x[-seq_len(n)])))
  a <- cbind(x[seq_len(n)], y[seq_len(n)])
  stopifnot(all(is.finite(a)), all(a>=0), all(a<=3), all(a==round(a)))
  a
})
long <- do.call(rbind,lapply(seq_along(traces), function(i) data.frame(trace=ids(e)[i],
  frame=seq_len(nrow(traces[[i]])), CCAT1=traces[[i]][,1], MYC=traces[[i]][,2])))
write.csv(long,file.path(out,"experimental_reporters_long.csv"),row.names=FALSE)
write.csv(data.frame(trace=ids(e),nframes=vapply(traces,nrow,integer(1))),file.path(out,"experimental_lengths.csv"),row.names=FALSE)
set.seed(20261106)
selections <- replicate(B,sample(seq_along(traces),length(traces),replace=TRUE),simplify=FALSE)
all.rows <- list(); k <- 0
for (obs in c("binary","reporter")) for (pair in c("EE","GG","EG")) {
  ij <- switch(pair, EE=c(1,1), GG=c(2,2), EG=c(1,2))
  curves <- lapply(traces,function(a) {
    if (obs=="binary") a <- 1*(a>0)
    env$multi_tau_xcorr(a[,ij[1]]+epsilon, a[,ij[2]]+epsilon,
                       framerate=1/(60*interval.min),M=4,bg_correct=FALSE)
  })
  avg <- env$coulon_average(curves)
  boots <- vapply(selections,function(pick) env$coulon_average(curves[pick])$Gn_mean, numeric(length(avg$lag_s)))
  k <- k+1
  all.rows[[k]] <- data.frame(observable=obs,pair=pair,lag_frames=round(avg$lag_s/(60*interval.min)),
    lag_min=avg$lag_s/60,mean=avg$Gn_mean,se=apply(boots,1,sd),n_traces=length(traces),
    se_kind="paired_trace_bootstrap",bootstrap_replicates=B)
}
write.csv(do.call(rbind,all.rows),file.path(out,"experimental.csv"),row.names=FALSE,na="NA")
stats <- function(selected) {
  a <- do.call(rbind,traces[selected]); c(on_CCAT1=mean(a[,1]>0),on_MYC=mean(a[,2]>0),
    reporters_CCAT1=mean(a[,1]),reporters_MYC=mean(a[,2]))
}
st <- stats(seq_along(traces)); sb <- vapply(selections,stats,numeric(4))
write.csv(data.frame(metric=names(st),mean=st,se=apply(sb,1,sd),se_kind="paired_trace_bootstrap"),
          file.path(out,"experimental_metrics.csv"),row.names=FALSE)
cat("Computed six experimental curves from",length(traces),"paired traces; lengths",paste(sort(unique(vapply(traces,nrow,integer(1)))),collapse=","),"\n")
