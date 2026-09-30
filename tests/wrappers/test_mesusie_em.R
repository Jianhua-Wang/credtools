# Small real-R integration/permutation checks, without modifying MESuSiE.
args <- commandArgs(trailingOnly = TRUE)
library(MESuSiE)
source(file.path(args[1], "mesusie_em.R"))
original_ser <- body(MESuSiE:::single_effect_regression)
original_core <- body(MESuSiE::meSuSie_core)
cache <- file.path(args[2], "cache")
set.seed(664)
for (n in c(2L, 3L, 5L)) {
  p <- 16L
  ss <- lapply(seq_len(n), function(k) {
    z <- c(8 + k/10, rep(.1, p - 1))
    data.frame(SNP = paste0("s", seq_len(p)), Beta = z / 100,
               Se = rep(.01, p), Z = z, N = rep(10000, p))
  })
  ld <- replicate(n, diag(p), simplify = FALSE)
  names(ss) <- names(ld) <- paste0("P", seq_len(n))
  results <- list()
  for (order in list(seq_len(n), rev(seq_len(n)))) {
    a <- build_credtools_mesusie(args[1], cache, tol = .001, coverage = .9)
    posterior_calls <- 0L
    adapter_env <- environment(a$core)
    adapter_env$mvlmm_reg <- local({
      stable <- adapter_env$mvlmm_reg
      function(...) { posterior_calls <<- posterior_calls + 1L; stable(...) }
    })
    f <- a$core(ld[order], ss[order], L = 2L, max_iter = 50L,
                estimate_residual_variance = FALSE, cor_threshold = 0)
    stopifnot(all(is.finite(f$pip)), all(f$pip >= 0 & f$pip <= 1),
              isTRUE(a$telemetry$converged), f$cs$requested_coverage == .9,
              a$telemetry$inner_calls > 0L, posterior_calls > 0L)
    results[[length(results) + 1L]] <- list(pip = f$pip, cs = f$cs$cs)
  }
  stopifnot(max(abs(results[[1]]$pip - results[[2]]$pip)) < 1e-6,
            identical(results[[1]]$cs, results[[2]]$cs))
}
# Exercise the original optimizer and nondefault controls on a cheap two-pop fit.
a <- build_credtools_mesusie(args[1], cache, optimizer = "native", coverage = .8)
f <- a$core(ld[1:2], ss[1:2], L = 1L, max_iter = 3L,
            estimate_residual_variance = FALSE, cor_threshold = 0)
stopifnot(f$cs$requested_coverage == .8, a$telemetry$inner_calls == 0L,
          identical(original_ser, body(MESuSiE:::single_effect_regression)),
          identical(original_core, body(MESuSiE::meSuSie_core)))
cat("MESUSIE_EM_INTEGRATION_PASS\n")
