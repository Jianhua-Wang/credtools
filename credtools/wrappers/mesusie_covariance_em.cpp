// Covariance EM for the MESuSiE SER mixture, validated in CRED-633/T2DF-362.
// Inactive ancestry effects are latent, not zero-valued observations.
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List mes_covariance_em(const arma::mat& beta, const arma::mat& se2,
                            arma::mat V, const arma::mat& prior,
                            Rcpp::List configs, int maxiter=100,
                            double tol=1e-9) {
  const int p=beta.n_rows, n=beta.n_cols, nc=configs.size();
  if (p < 1 || n < 1 || nc < 1 || maxiter < 1 || !std::isfinite(tol) || tol <= 0 ||
      se2.n_rows != beta.n_rows || se2.n_cols != beta.n_cols ||
      V.n_rows != beta.n_cols || V.n_cols != beta.n_cols ||
      prior.n_rows != beta.n_rows || prior.n_cols != (unsigned)nc ||
      !beta.is_finite() || !se2.is_finite() || !V.is_finite() ||
      !prior.is_finite() || arma::any(arma::vectorise(se2) <= 0) ||
      arma::any(arma::vectorise(prior) < 0) || arma::accu(prior) <= 0) {
    Rcpp::stop("Invalid covariance-EM dimensions, data, weights or controls");
  }
  double previous=-arma::datum::inf, objective=previous;
  int iteration=0, status=5;
  arma::mat lbf(p,nc);
  arma::cube second(n,n,p*nc);
  for(iteration=0;iteration<maxiter;iteration++) {
    for(int k=0;k<nc;k++) {
      arma::uvec active=Rcpp::as<arma::uvec>(configs[k]);
      if (active.is_empty() || active.min() < 1 || active.max() > (unsigned)n)
        Rcpp::stop("Invalid ancestry configuration");
      active-=1;
      arma::mat VA=V.cols(active), VV=V.submat(active,active);
      for(int j=0;j<p;j++) {
        arma::vec y=beta.row(j).t(), s=se2.row(j).t();
        y=y.elem(active);s=s.elem(active);
        arma::mat T=VV+arma::diagmat(s);
        arma::mat precision=arma::inv_sympd(T);
        arma::vec m=VA*precision*y;
        arma::mat C=V-VA*precision*VA.t();
        second.slice(j+p*k)=arma::symmatu(C+m*m.t());
        lbf(j,k)=0.5*(arma::accu(arma::square(y)/s)-
          arma::as_scalar(y.t()*precision*y)+arma::accu(arma::log(s))-
          arma::log_det_sympd(T));
      }
    }
    double shift=lbf.max();
    arma::mat w=arma::exp(lbf-shift)%prior;
    double mass=arma::accu(w);
    if (!std::isfinite(shift) || !std::isfinite(mass) || mass <= 0)
      Rcpp::stop("Nonfinite covariance-EM objective");
    w/=mass;
    objective=shift+std::log(mass);
    if(iteration>0 && std::abs(objective-previous)<tol) {status=3;break;}
    if(iteration==maxiter-1)break;
    previous=objective;
    arma::mat next(n,n,arma::fill::zeros);
    for(int k=0;k<nc;k++)for(int j=0;j<p;j++)next+=w(j,k)*second.slice(j+p*k);
    // This floor is on effect covariance, not an LD-matrix repair.
    arma::vec eig;arma::mat U;arma::eig_sym(eig,U,arma::symmatu(next));
    eig.transform([](double x){return std::max(x,1e-12);});
    V=arma::symmatu(U*arma::diagmat(eig)*U.t());
  }
  return Rcpp::List::create(Rcpp::Named("V")=V,
      Rcpp::Named("objective")=-objective,Rcpp::Named("iterations")=iteration+1,
      Rcpp::Named("status")=status);
}
