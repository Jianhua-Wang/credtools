#!/usr/bin/env Rscript
# SuSiE-inf wrapper script for credtools.
# Called by credtools/wrappers/susie_inf.py via subprocess.
#
# Input files (in temp_dir):
#   sumstats.csv - SNP, BHAT, SHAT (per-SNP marginal effect and SE)
#   ld.bin       - LD matrix in float64 binary format (column-major)
#   ld_dim.txt   - LD matrix dimension (single integer)
#   snpids.txt   - SNP ID list (one per line)
#
# Output files (in temp_dir):
#   susie_inf_pips.csv   - SNP, PIP
#   susie_inf_cs.csv     - CS_ID, SNP
#   susie_inf_status.txt - tab-separated key/value with `converged` and `niter`

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  params <- list()
  i <- 1
  while (i <= length(args)) {
    if (startsWith(args[i], "--")) {
      key <- sub("^--", "", args[i])
      if (i + 1 <= length(args) && !startsWith(args[i + 1], "--")) {
        params[[key]] <- args[i + 1]
        i <- i + 2
      } else {
        params[[key]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  return(params)
}

params <- parse_args(args)

temp_dir                   <- params$temp_dir
n                          <- as.integer(params$n)
L                          <- as.integer(params$L)
coverage                   <- as.numeric(params$coverage)
max_iter                   <- as.integer(params$max_iter)
estimate_residual_variance <- as.logical(params$estimate_residual_variance)
min_abs_corr               <- as.numeric(params$min_abs_corr)
tol                        <- as.numeric(params$tol)
# Use `optim` for the prior estimator: it is the only estimator supported by
# all `unmappable_effects` modes, so we hardcode it for symmetry with the
# SuSiE-ash wrapper.
estimate_prior_method      <- "optim"

cat("SuSiE-inf wrapper: starting\n")
cat(sprintf("  temp_dir: %s\n", temp_dir))
cat(sprintf("  n: %d\n", n))
cat(sprintf("  L: %d\n", L))
cat(sprintf("  coverage: %f\n", coverage))
cat(sprintf("  max_iter: %d\n", max_iter))

if (!requireNamespace("susieR", quietly = TRUE)) {
  stop("susieR R package is not installed. Please install it with:\n",
       "  install.packages('susieR')\n",
       "Or for the development version with SuSiE-inf:\n",
       "  remotes::install_github('stephenslab/susieR')")
}
suppressMessages(library(susieR))

ss <- read.csv(file.path(temp_dir, "sumstats.csv"), stringsAsFactors = FALSE)
snpids <- readLines(file.path(temp_dir, "snpids.txt"))
snpids <- snpids[nchar(snpids) > 0]
n_snps <- as.integer(trimws(readLines(file.path(temp_dir, "ld_dim.txt"))[1]))
ld_raw <- readBin(file.path(temp_dir, "ld.bin"),
                  what = "double", n = n_snps * n_snps, size = 8)
R <- matrix(ld_raw, nrow = n_snps, ncol = n_snps, byrow = FALSE)

if (!"unmappable_effects" %in% names(formals(susieR::susie_rss))) {
  stop("Installed susieR (", as.character(packageVersion("susieR")),
       ") does not support `unmappable_effects`. ",
       "Please upgrade to susieR >= 0.16.1 (the SuSiE 2.0 release).")
}

cat(sprintf("Running susie_rss(unmappable_effects='inf') with %d SNPs...\n",
            n_snps))

fit <- tryCatch({
  susieR::susie_rss(
    bhat                       = ss$BHAT,
    shat                       = ss$SHAT,
    n                          = n,
    R                          = R,
    L                          = L,
    coverage                   = coverage,
    max_iter                   = max_iter,
    estimate_residual_variance = estimate_residual_variance,
    min_abs_corr               = min_abs_corr,
    tol                        = tol,
    estimate_prior_method      = estimate_prior_method,
    unmappable_effects         = "inf"
  )
}, error = function(e) {
  cat(sprintf("Error in susie_rss: %s\n", e$message))
  stop(e)
})

cat("SuSiE-inf completed.\n")

pip_vec <- as.numeric(fit$pip)
if (length(pip_vec) != length(snpids)) {
  stop(sprintf(
    "PIP length (%d) does not match SNP count (%d)",
    length(pip_vec), length(snpids)
  ))
}
pip_df <- data.frame(SNP = snpids, PIP = pip_vec, stringsAsFactors = FALSE)
write.csv(pip_df, file.path(temp_dir, "susie_inf_pips.csv"), row.names = FALSE)

cs_rows <- data.frame(CS_ID = integer(), SNP = character(), stringsAsFactors = FALSE)
sets <- fit$sets
if (!is.null(sets) && !is.null(sets$cs) && length(sets$cs) > 0) {
  cs_list <- sets$cs
  for (cs_i in seq_along(cs_list)) {
    idx <- as.integer(unlist(cs_list[[cs_i]]))
    idx <- idx[idx >= 1 & idx <= length(snpids)]
    if (length(idx) == 0) next
    cs_rows <- rbind(
      cs_rows,
      data.frame(
        CS_ID = cs_i,
        SNP = snpids[idx],
        stringsAsFactors = FALSE
      )
    )
  }
}
write.csv(cs_rows, file.path(temp_dir, "susie_inf_cs.csv"), row.names = FALSE)

converged <- isTRUE(fit$converged)
niter_val <- if (!is.null(fit$niter)) as.integer(fit$niter) else NA_integer_
status_path <- file.path(temp_dir, "susie_inf_status.txt")
writeLines(
  c(
    sprintf("converged\t%s", if (converged) "TRUE" else "FALSE"),
    sprintf("niter\t%s", if (is.na(niter_val)) "" else niter_val)
  ),
  status_path
)

cat(sprintf(
  "SuSiE-inf wrapper: done. %d credible sets, converged=%s, niter=%s.\n",
  length(unique(cs_rows$CS_ID)),
  if (converged) "TRUE" else "FALSE",
  if (is.na(niter_val)) "NA" else niter_val
))
