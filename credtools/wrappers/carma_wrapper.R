#!/usr/bin/env Rscript
# CARMA wrapper script for credtools
# Called by credtools/wrappers/carma.py via subprocess
#
# Input files (in temp_dir):
#   sumstats.csv - SNP, Z (z-scores)
#   ld.bin       - LD matrix in float64 binary format (column-major)
#   ld_dim.txt   - LD matrix dimension (single integer)
#   snpids.txt   - SNP ID list (one per line)
#
# Output files (in temp_dir):
#   carma_pips.csv      - SNP, PIP
#   carma_cs.csv        - CS_ID, SNP
#   carma_outliers.csv  - SNP

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

temp_dir            <- params$temp_dir
num_causal          <- as.integer(params$num_causal)
rho_index           <- as.numeric(params$rho_index)
bf_index            <- as.numeric(params$bf_index)
outlier_switch      <- as.logical(params$outlier_switch)
outlier_bf_index    <- as.numeric(params$outlier_bf_index)
effect_size_prior   <- as.character(params$effect_size_prior)
em_dist             <- as.character(params$em_dist)
max_model_dim       <- as.integer(params$max_model_dim)
all_iter            <- as.integer(params$all_iter)
all_inner_iter      <- as.integer(params$all_inner_iter)
input_alpha         <- as.numeric(params$input_alpha)
epsilon_threshold   <- as.numeric(params$epsilon_threshold)
tau                 <- as.numeric(params$tau)
y_var               <- as.numeric(params$y_var)

cat("CARMA wrapper: starting\n")
cat(sprintf("  temp_dir: %s\n", temp_dir))
cat(sprintf("  num_causal: %d\n", num_causal))
cat(sprintf("  rho_index: %f\n", rho_index))
cat(sprintf("  outlier_switch: %s\n", outlier_switch))

if (!requireNamespace("CARMA", quietly = TRUE)) {
  stop("CARMA R package is not installed. Please install it with:\n",
       "  devtools::install_github('ZikunY/CARMA')")
}
suppressMessages(library(CARMA))

# Load inputs
ss <- read.csv(file.path(temp_dir, "sumstats.csv"), stringsAsFactors = FALSE)
snpids <- readLines(file.path(temp_dir, "snpids.txt"))
snpids <- snpids[nchar(snpids) > 0]
n_snps <- as.integer(trimws(readLines(file.path(temp_dir, "ld_dim.txt"))[1]))
ld_raw <- readBin(file.path(temp_dir, "ld.bin"),
                  what = "double", n = n_snps * n_snps, size = 8)
ld <- matrix(ld_raw, nrow = n_snps, ncol = n_snps, byrow = FALSE)

# CARMA expects list-form inputs (it supports meta-analysis across studies).
z.list <- list(ss$Z)
ld.list <- list(ld)
lambda.list <- list(1)

cat(sprintf("Running CARMA with %d SNPs...\n", n_snps))

fit <- tryCatch({
  CARMA(
    z.list = z.list,
    ld.list = ld.list,
    lambda.list = lambda.list,
    effect.size.prior = effect_size_prior,
    rho.index = rho_index,
    BF.index = bf_index,
    EM.dist = em_dist,
    Max.Model.Dim = max_model_dim,
    all.iter = all_iter,
    all.inner.iter = all_inner_iter,
    input.alpha = input_alpha,
    epsilon.threshold = epsilon_threshold,
    num.causal = num_causal,
    y.var = y_var,
    tau = tau,
    outlier.switch = outlier_switch,
    outlier.BF.index = outlier_bf_index,
    output.labels = NULL,
    printing.log = FALSE
  )
}, error = function(e) {
  cat(sprintf("Error in CARMA: %s\n", e$message))
  stop(e)
})

cat("CARMA completed.\n")

# CARMA returns a list of length n_studies; we use the first (and only) study.
res <- fit[[1]]

pip_vec <- as.numeric(res$PIP)
if (length(pip_vec) != length(snpids)) {
  stop(sprintf(
    "PIP length (%d) does not match SNP count (%d)",
    length(pip_vec), length(snpids)
  ))
}
pip_df <- data.frame(SNP = snpids, PIP = pip_vec, stringsAsFactors = FALSE)
write.csv(pip_df, file.path(temp_dir, "carma_pips.csv"), row.names = FALSE)

cs_rows <- data.frame(CS_ID = integer(), SNP = character(), stringsAsFactors = FALSE)
cs_obj <- res[["Credible set"]]
if (!is.null(cs_obj) && length(cs_obj) >= 2) {
  # cs_obj[[2]] is a list of credible sets; each element is a vector of SNP indices.
  cs_list <- cs_obj[[2]]
  if (length(cs_list) > 0) {
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
}
write.csv(cs_rows, file.path(temp_dir, "carma_cs.csv"), row.names = FALSE)

outlier_snps <- character(0)
outlier_obj <- res$Outliers
if (!is.null(outlier_obj)) {
  if (is.list(outlier_obj) && length(outlier_obj) >= 1) {
    cand <- outlier_obj[[1]]
  } else {
    cand <- outlier_obj
  }
  cand_idx <- suppressWarnings(as.integer(unlist(cand)))
  cand_idx <- cand_idx[!is.na(cand_idx) & cand_idx >= 1 & cand_idx <= length(snpids)]
  if (length(cand_idx) > 0) {
    outlier_snps <- snpids[unique(cand_idx)]
  }
}
write.csv(
  data.frame(SNP = outlier_snps, stringsAsFactors = FALSE),
  file.path(temp_dir, "carma_outliers.csv"),
  row.names = FALSE
)

cat(sprintf(
  "CARMA wrapper: done. %d credible sets, %d outliers.\n",
  length(unique(cs_rows$CS_ID)), length(outlier_snps)
))
