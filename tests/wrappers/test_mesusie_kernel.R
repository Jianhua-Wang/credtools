# Numerical equivalence to the frozen v0.9.7 EM, including the generic path.
args <- commandArgs(trailingOnly=TRUE)
reference <- new.env(); candidate <- new.env()
Rcpp::sourceCpp(args[2], env=reference, cacheDir=file.path(args[3],"reference"))
Rcpp::sourceCpp(file.path(args[1],"mesusie_covariance_em.cpp"),
                env=candidate, cacheDir=file.path(args[3],"candidate"))
configs_for <- function(n) unlist(lapply(seq_len(n),function(k)
  combn(seq_len(n),k,simplify=FALSE)),recursive=FALSE)
set.seed(20260929)
for(n in 1:7) for(scenario in 1:3) {
  p <- 8L; configs <- lapply(configs_for(n),matrix)
  A <- matrix(rnorm(n*n),n,n)
  V <- crossprod(A)*.001+diag(.0001,n)
  if(scenario==2) V <- matrix(.003,n,n)+diag(1e-10,n)
  prior <- matrix(runif(p*length(configs),.2,2),p,length(configs))
  if(scenario==3 && n>1) prior[,1] <- 0
  prior <- prior/sum(prior)
  beta <- matrix(rnorm(p*n,sd=.04),p,n)
  se2 <- matrix(runif(p*n,.0001,.006),p,n)
  a <- reference$mes_covariance_em(beta,se2,V,prior,configs,100L,1e-9)
  b <- candidate$mes_covariance_em(beta,se2,V,prior,configs,100L,1e-9)
  stopifnot(max(abs(a$V-b$V))<1e-8,abs(a$objective-b$objective)<1e-6,
            a$iterations==b$iterations,a$status==b$status)
}
# Input failures must not turn into plausible-looking results.
bad <- se2;bad[1,1] <- 0
stopifnot(inherits(try(candidate$mes_covariance_em(beta,bad,V,prior,configs),silent=TRUE),"try-error"))
bad <- beta;bad[1,1] <- NaN
stopifnot(inherits(try(candidate$mes_covariance_em(bad,se2,V,prior,configs),silent=TRUE),"try-error"))
stopifnot(inherits(try(candidate$mes_covariance_em(beta,se2,V,prior*0,configs),silent=TRUE),"try-error"))
cat("MESUSIE_KERNEL_EQUIVALENCE_PASS\n")
