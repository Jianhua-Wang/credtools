# Process-local MESuSiE adapter: never modify the installed R namespace.
# Preserve the Gaussian model; use stable likelihood/posterior arithmetic for EM.
build_credtools_mesusie <- function(script_dir, cache_dir, optimizer = "em",
                                    tol = 1e-3, coverage = .95,
                                    em_max_iter = 100L, em_tol = 1e-9) {
  stopifnot(optimizer %in% c("em", "native"), is.finite(tol), tol > 0,
            is.finite(coverage), coverage > 0, coverage < 1,
            is.finite(em_tol), em_tol > 0, em_max_iter >= 1)
  env <- new.env(parent = asNamespace("MESuSiE"))
  telemetry <- new.env(parent = emptyenv())
  telemetry$converged <- FALSE
  telemetry$n_iter <- 0L
  telemetry$inner_calls <- 0L
  telemetry$inner_maxiter <- 0L
  telemetry$fallbacks <- 0L
  env$telemetry <- telemetry
  env$outer_tol <- tol
  env$check_elbo <- function(new, old, iter) {
    if (!is.finite(new)) stop("MESuSiE produced a nonfinite ELBO")
    delta <- new - old
    if (is.finite(delta) && delta < -tol)
      warning(sprintf("ELBO decreased by %.3g at iteration %d", -delta, iter))
    is.finite(delta) && delta >= 0 && delta < tol
  }
  ser <- MESuSiE:::single_effect_regression
  if (optimizer == "em") {
    if (!requireNamespace("Rcpp", quietly = TRUE) ||
        !requireNamespace("RcppArmadillo", quietly = TRUE))
      stop("MESuSiE EM requires Rcpp, RcppArmadillo and a working R C++ toolchain")
    # The cache belongs to this invocation; parallel jobs never share its index.
    Rcpp::sourceCpp(file.path(script_dir, "mesusie_covariance_em.cpp"),
                    env = env, cacheDir = cache_dir, rebuild = FALSE)
    Rcpp::sourceCpp(file.path(script_dir, "mesusie_numerics.cpp"),
                    env = env, cacheDir = cache_dir, rebuild = FALSE)
    env$mvlmm_reg <- function(betahat, shat2, V_mat)
      env$mes_stable_mvlmm(betahat, shat2, V_mat)
    objective <- function(v, beta, se2, prior, configs) {
      env$mes_stable_loglik(beta, se2, v, prior, configs)
    }
    env$update_cov <- function(beta, se2, obj, effect) {
      old <- obj$V[[effect]]
      have <- all(diag(old) > 0)
      initial <- if (have) old else
        diag(pmax(1e-10, apply(beta^2 - se2, 2, max)), ncol(beta))
      fit <- env$mes_covariance_em(beta, se2, initial, obj$pi,
                                    obj$column_config, em_max_iter, em_tol)
      value <- objective(fit$V, beta, se2, obj$pi, obj$column_config)
      if (!is.finite(value) || abs(value - fit$objective) >= 1e-6)
        stop("MESuSiE EM objective disagrees with the stable Gaussian likelihood")
      oldvalue <- if (have) objective(old, beta, se2, obj$pi, obj$column_config) else Inf
      telemetry$inner_calls <- telemetry$inner_calls + 1L
      telemetry$inner_maxiter <- telemetry$inner_maxiter + as.integer(fit$status == 5L)
      if (is.finite(oldvalue) && oldvalue < value) {
        telemetry$fallbacks <- telemetry$fallbacks + 1L
        return(old)
      }
      fit$V
    }
    statements <- as.list(body(ser))
    at <- which(vapply(statements, function(x) is.call(x) &&
      identical(x[[1]], as.name("if")) &&
      grepl("estimate_prior_method", paste(deparse(x[[2]]), collapse = ""), fixed = TRUE), logical(1)))
    if (length(at) != 1L)
      stop("Unsupported MESuSiE single_effect_regression; cannot safely install EM adapter")
    statements[[at]] <- quote({V_mat <- update_cov(betahat, shat2, meSuSieObject_obj, l_index)})
    body(ser) <- as.call(statements)
  }
  environment(ser) <- env
  env$single_effect_regression <- ser
  core <- MESuSiE::meSuSie_core
  counts <- c(stop = 0L, output = 0L, cs = 0L)
  rewrite <- function(node) {
    if (!is.call(node)) return(node)
    if (identical(node[[1]], as.name("if")) &&
        grepl("ELBO", paste(deparse(node[[2]]), collapse = ""), fixed = TRUE)) {
      node[[2]] <- quote(check_elbo(meSuSieObject_obj$ELBO[iter + 1],
                                    meSuSieObject_obj$ELBO[iter], iter))
      node[[3]] <- quote({telemetry$converged <- TRUE; break})
      counts["stop"] <<- counts["stop"] + 1L
    }
    if (identical(node[[1]], as.name("meSuSie_get_cs"))) {
      node$coverage <- coverage
      counts["cs"] <<- counts["cs"] + 1L
    }
    if (identical(node[[1]], as.name("return")) &&
        identical(node[[2]], as.name("meSuSieObject_obj"))) {
      counts["output"] <<- counts["output"] + 1L
      return(quote({telemetry$n_iter <- iter
        telemetry$elbo <- meSuSieObject_obj$ELBO[seq_len(iter + 1L)]
        return(meSuSieObject_obj)}))
    }
    if (length(node) > 1L) for (i in 2:length(node))
      if (!identical(node[[i]], quote(expr = ))) node[i] <- list(rewrite(node[[i]]))
    node
  }
  body(core) <- rewrite(body(core))
  if (any(counts != 1L))
    stop("Unsupported MESuSiE meSuSie_core; convergence/CS adapter pattern changed")
  environment(core) <- env
  list(core = core, telemetry = telemetry)
}
