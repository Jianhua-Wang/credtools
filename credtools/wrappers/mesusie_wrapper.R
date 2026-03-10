#!/usr/bin/env Rscript
# MESuSiE wrapper script for credtools
# Called by credtools/wrappers/mesusie.py via subprocess
#
# Input files (in temp_dir):
#   pop_{i}_sumstats.csv  - Summary stats per population (SNP, Beta, Se, Z, N)
#   pop_{i}_ld.bin        - LD matrix in float64 binary format
#   pop_{i}_ld_dim.txt    - LD matrix dimension (single integer)
#   pop_{i}_snpids.txt    - SNP ID list (one per line)
#   pop_names.txt         - Population names (one per line)
#
# Output files (in temp_dir):
#   mesusie_pips.csv       - (SNP, PIP)
#   mesusie_cs.csv         - (CS_ID, SNP)
#   mesusie_purity.csv     - (CS_ID, PURITY, CS_TYPE)
#   mesusie_converged.txt  - "TRUE" or "FALSE"

# Parse command line arguments
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

# Required parameters
temp_dir <- params$temp_dir
n_pop <- as.integer(params$n_pop)
L <- as.integer(params$L)
coverage <- as.numeric(params$coverage)
max_iter <- as.integer(params$max_iter)
purity <- as.numeric(params$purity)
estimate_residual_variance <- as.logical(params$estimate_residual_variance)

cat("MESuSiE wrapper: starting\n")
cat(sprintf("  temp_dir: %s\n", temp_dir))
cat(sprintf("  n_pop: %d\n", n_pop))
cat(sprintf("  L: %d\n", L))
cat(sprintf("  coverage: %f\n", coverage))
cat(sprintf("  max_iter: %d\n", max_iter))
cat(sprintf("  purity: %f\n", purity))
cat(sprintf("  estimate_residual_variance: %s\n", estimate_residual_variance))

# Check MESuSiE package
if (!requireNamespace("MESuSiE", quietly = TRUE)) {
  stop("MESuSiE R package is not installed. Please install it with:\n",
       "  devtools::install_github('borangao/MESuSiE')")
}
library(MESuSiE)

# Read population names
pop_names_file <- file.path(temp_dir, "pop_names.txt")
pop_names <- readLines(pop_names_file)

# Read input files for each population
summary_stat_list <- list()
R_mat_list <- list()
snpid_list <- list()

for (i in seq_len(n_pop)) {
  idx <- i - 1  # Python uses 0-based indexing
  pop_name <- pop_names[i]

  # Read summary stats
  ss_file <- file.path(temp_dir, sprintf("pop_%d_sumstats.csv", idx))
  ss <- read.csv(ss_file, stringsAsFactors = FALSE)

  # Read SNP IDs
  snpid_file <- file.path(temp_dir, sprintf("pop_%d_snpids.txt", idx))
  snpids <- readLines(snpid_file)
  # Remove empty trailing lines
  snpids <- snpids[nchar(snpids) > 0]

  # Read LD matrix dimension
  dim_file <- file.path(temp_dir, sprintf("pop_%d_ld_dim.txt", idx))
  n_snps <- as.integer(trimws(readLines(dim_file)[1]))

  # Read LD matrix (float64 binary)
  ld_file <- file.path(temp_dir, sprintf("pop_%d_ld.bin", idx))
  ld_raw <- readBin(ld_file, what = "double", n = n_snps * n_snps, size = 8)
  R_mat <- matrix(ld_raw, nrow = n_snps, ncol = n_snps, byrow = FALSE)

  # Build summary stat dataframe for MESuSiE
  # MESuSiE expects dataframe with columns: SNP, Beta, Se, Z, N
  # (used by meSuSieData class: XtX_diag uses Z/Se/N, Xty_pro uses Beta)
  summary_stat_list[[pop_name]] <- data.frame(
    SNP = ss$SNP,
    Beta = ss$Beta,
    Se = ss$Se,
    Z = ss$Z,
    N = ss$N,
    stringsAsFactors = FALSE
  )
  R_mat_list[[pop_name]] <- R_mat
  snpid_list[[i]] <- snpids

  cat(sprintf("  Population %s: %d SNPs, N=%d\n", pop_name, n_snps, ss$N[1]))
}

# Use the SNP IDs from the first population (all should be aligned by Python)
snpids <- snpid_list[[1]]

cat("Running meSuSie_core...\n")

# Call meSuSie_core
# Signature: meSuSie_core(R_mat_list, summary_stat_list, L,
#   residual_variance, prior_weights, ancestry_weight,
#   optim_method, estimate_residual_variance, max_iter,
#   cor_method, cor_threshold)
fit <- tryCatch({
  meSuSie_core(
    R_mat_list = R_mat_list,
    summary_stat_list = summary_stat_list,
    L = L,
    estimate_residual_variance = estimate_residual_variance,
    max_iter = max_iter,
    cor_threshold = purity
  )
}, error = function(e) {
  cat(sprintf("Error in meSuSie_core: %s\n", e$message))
  stop(e)
})

cat("meSuSie_core completed.\n")

# Extract PIPs from the fit object
pip <- fit$pip
pip_df <- data.frame(SNP = snpids, PIP = pip, stringsAsFactors = FALSE)
write.csv(pip_df, file.path(temp_dir, "mesusie_pips.csv"), row.names = FALSE)

# Extract credible sets
# meSuSie_core stores the full meSuSie_get_cs output in fit$cs
# fit$cs is a list with: $cs, $cs_category, $purity, $cs_index, $coverage, $requested_coverage
cs_result <- fit$cs

# Write credible set results
cs_rows <- data.frame(CS_ID = integer(), SNP = character(), stringsAsFactors = FALSE)
purity_rows <- data.frame(CS_ID = integer(), PURITY = numeric(), CS_TYPE = character(), stringsAsFactors = FALSE)

if (!is.null(cs_result$cs) && length(cs_result$cs) > 0) {
  cs_list <- cs_result$cs
  cs_names <- names(cs_list)
  for (cs_i in seq_along(cs_list)) {
    cs_snp_indices <- unlist(cs_list[[cs_i]])
    cs_snp_names <- snpids[cs_snp_indices]
    for (snp in cs_snp_names) {
      cs_rows <- rbind(cs_rows, data.frame(CS_ID = cs_i, SNP = snp, stringsAsFactors = FALSE))
    }

    # Extract purity (min.abs.corr column from purity dataframe)
    cs_purity <- NA
    if (!is.null(cs_result$purity) && is.data.frame(cs_result$purity)) {
      cs_name <- cs_names[cs_i]
      if (cs_name %in% rownames(cs_result$purity)) {
        cs_purity <- cs_result$purity[cs_name, "min.abs.corr"]
      }
    }

    # Extract CS type (shared/ancestry-specific) from cs_category
    cs_type <- "unknown"
    if (!is.null(cs_result$cs_category)) {
      cs_name <- cs_names[cs_i]
      if (cs_name %in% names(cs_result$cs_category)) {
        cs_type <- as.character(cs_result$cs_category[cs_name])
      }
    }

    purity_rows <- rbind(purity_rows, data.frame(
      CS_ID = cs_i, PURITY = cs_purity, CS_TYPE = cs_type, stringsAsFactors = FALSE
    ))
  }
}

write.csv(cs_rows, file.path(temp_dir, "mesusie_cs.csv"), row.names = FALSE)
write.csv(purity_rows, file.path(temp_dir, "mesusie_purity.csv"), row.names = FALSE)

# Write convergence status
converged <- ifelse(!is.null(fit$converged), fit$converged, TRUE)
writeLines(as.character(toupper(converged)), file.path(temp_dir, "mesusie_converged.txt"))

cat(sprintf("MESuSiE wrapper: done. Found %d credible sets.\n", nrow(purity_rows)))
