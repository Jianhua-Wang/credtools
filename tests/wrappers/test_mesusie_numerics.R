# Compare stable Gaussian evaluations to native and independent analytic cases.
args <- commandArgs(trailingOnly = TRUE)
library(MESuSiE)
Rcpp::sourceCpp(file.path(args[1], "mesusie_numerics.cpp"),
                cacheDir = file.path(args[2], "numeric_cache"))
vec_cov <- function(v) {
  r <- v / sqrt(outer(diag(v), diag(v)))
  diag(r) <- log(diag(v))
  r[upper.tri(r, diag = TRUE)]
}
set.seed(991)
cases <- 0L
for (n in 1:7) {
  p <- 20L
  b <- matrix(rnorm(p*n), p, n)
  s <- matrix(runif(p*n, .2, 2), p, n)
  z <- matrix(rnorm(n*n), n, n)
  V <- z %*% t(z) / n + diag(.5, n)
  stable <- mes_stable_mvlmm(b, s, V)
  native <- if (n == 1L) MESuSiE:::uni_reg(b, s, V[1,1]) else
    MESuSiE:::mvlmm_reg(b, s, V)
  for (key in c("lbf", "post_mean", "post_mean2"))
    stopifnot(max(abs(stable[[key]] - native[[key]])) < 1e-9)
  configs <- unlist(lapply(seq_len(n), function(k) combn(n, k, simplify=FALSE)),
                    recursive=FALSE)
  prior <- matrix(runif(p*length(configs)), p)
  prior <- prior / sum(prior)
  prior[1,] <- 0
  old <- MESuSiE:::loglik_cpp(vec_cov(V), b, s, prior, n,
                             matrix(cumsum(seq_len(n))), configs)
  value <- mes_stable_loglik(b, s, V, prior, configs)
  stopifnot(abs(value-old) < 1e-9)
  for (j in seq_len(p)) {
    P <- solve(solve(V) + diag(1/s[j,], n))
    m <- P %*% (b[j,] / s[j,])
    stopifnot(max(abs(stable$post_mean[j,]-m)) < 1e-9,
              max(abs(stable$post_mean2[j,]-m^2-diag(P))) < 1e-9)
  }
  cases <- cases + 1L
}
# Zero and rank-deficient prior covariances are valid without added variance.
b <- matrix(c(1, 2, -1, 3), 2, 2)
s <- matrix(.01, 2, 2)
zero <- mes_stable_mvlmm(b, s, matrix(0, 2, 2))
stopifnot(max(abs(zero$lbf)) < 1e-10, all(zero$post_mean == 0),
          all(zero$post_mean2 == 0))
V <- matrix(1, 2, 2)
rankone <- mes_stable_mvlmm(b, s, V)
m <- rowSums(b/s)/(1+rowSums(1/s))
variance <- 1/(1+rowSums(1/s))
stopifnot(max(abs(rankone$post_mean-m)) < 1e-9,
          max(abs(rankone$post_mean2-m^2-variance)) < 1e-9)
# Ill-conditioned full-rank prior: determinant-lemma/Sherman-Morrison reference.
# Independent closed forms avoid both the native determinant and a Cholesky solve.
v <- c(100, 14, 75, 110, 85)
epsilon <- 2e-8
V <- tcrossprod(v) + diag(epsilon, 5)
b <- rbind(c(.01, -.2, .003, -.005, .01), c(10, -25, 2, 4, -3))
s <- rbind(c(2e-5, .3, 8, 1e-4, .01), c(.001, .02, .5, .008, .01))
r <- mes_stable_mvlmm(b, s, V)
for (j in 1:2) {
  d <- s[j,] + epsilon
  denominator <- 1 + sum(v^2/d)
  m <- epsilon/d*b[j,] + s[j,]/d*v*sum(v*b[j,]/d)/denominator
  variance <- epsilon*s[j,]/d + (s[j,]/d*v)^2/denominator
  score <- .5*(sum(b[j,]/s[j,]*m) - sum(log1p(epsilon/s[j,])) -
                 log1p(sum(v^2/d)))
  stopifnot(max(abs(r$post_mean[j,]-m)) < 2e-6,
            max(abs(r$post_mean2[j,]-r$post_mean[j,]^2-variance)) < 1e-8,
            abs(r$lbf[j]-score) < 1e-7*max(1,abs(score)))
}
fails <- function(expr) inherits(tryCatch(force(expr), error=identity), "error")
stopifnot(fails(mes_stable_mvlmm(b[,1:2], s[,1:2], diag(c(1,-.01)))),
          fails(mes_stable_mvlmm(b, s*0, V)),
          fails(mes_stable_loglik(b, s, V, matrix(0,2,1), list(1:5))),
          fails(mes_stable_loglik(b, s, V, matrix(1,2,1), list(c(1L,1L)))))
cat("MESUSIE_NUMERICS_PASS", cases, "regular dimensions\n")
