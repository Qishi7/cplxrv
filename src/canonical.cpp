#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]
// Uses the R toolchain's default C++ standard.
using namespace Rcpp;

namespace {
void positive(double x, const char *what) {
  if (!std::isfinite(x) || x <= 0)
    stop("%s must be finite and positive", what);
}
double gamma_rate(double shape, double rate) {
  positive(shape, "gamma shape");
  positive(rate, "gamma rate");
  double x = R::rgamma(shape, 1.0 / rate);
  positive(x, "gamma draw");
  return x;
}
double inv_gamma(double shape, double rate) {
  double x = 1.0 / gamma_rate(shape, rate);
  positive(x, "inverse-gamma draw");
  return x;
}
// w ~ GIG(1/2, chi, psi). Reciprocal is IGau(sqrt(psi/chi),psi).
// Stable Michael-Schucany-Haas roots expressed directly on the w scale.
double local_w(double chi, double psi) {
  if (!std::isfinite(chi) || chi < 0)
    stop("GIG chi must be finite and nonnegative");
  positive(psi, "GIG psi");
  if (chi == 0)
    return gamma_rate(0.5, psi / 2.0);
  double b = std::sqrt(chi) / std::sqrt(psi);
  positive(b, "GIG root scale");
  double z = R::rnorm(0, 1), t = (0.5 * (z * z)) / psi;
  double hi = b + t + std::sqrt(t) * std::sqrt(2.0 * b + t);
  positive(hi, "GIG upper root");
  double w = R::runif(0, 1) < 1.0 / (1.0 + b / hi) ? hi : b * (b / hi);
  positive(w, "GIG draw");
  return w;
}
arma::mat covariance(int n, int d, double re, double im) {
  if (n < 1 || (d != 1 && d != 2))
    stop("invalid dimensions");
  if (!std::isfinite(re) || !std::isfinite(im) || re * re + im * im >= 1)
    stop("rho must satisfy |rho| < 1");
  if (d == 1) {
    if (re != 0 || im != 0)
      stop("real samplers require rho=0");
    return arma::eye(n, n);
  }
  arma::mat C(2 * n, 2 * n, arma::fill::zeros);
  C.submat(0, 0, n - 1, n - 1) = (1 + re) * arma::eye(n, n);
  C.submat(n, n, 2 * n - 1, 2 * n - 1) = (1 - re) * arma::eye(n, n);
  C.submat(0, n, n - 1, 2 * n - 1) = im * arma::eye(n, n);
  C.submat(n, 0, 2 * n - 1, n - 1) = im * arma::eye(n, n);
  return C;
}
arma::mat precision(const arma::mat &C) {
  arma::mat K;
  if (!arma::inv_sympd(K, C))
    stop("covariance Cholesky inversion failed; no ridge added");
  return K;
}
arma::vec B_values(const arma::vec &b, int p, int d) {
  arma::vec B = arma::square(b.head(p));
  if (d == 2)
    B += arma::square(b.tail(p));
  return B;
}
arma::vec coefficient(const arma::mat &T, const arma::vec &h, const arma::vec &w,
                      double s, int d) {
  // beta | w,s,y ~ N(A^{-1} h, s A^{-1}), A = H^T K H + D_w^{-1}.
  arma::vec iw = 1.0 / w;
  if (d == 2)
    iw = arma::join_cols(iw, iw);
  arma::mat A = T + arma::diagmat(iw), L;
  if (!A.is_finite() || !h.is_finite() || !arma::chol(L, A, "lower"))
    stop("coefficient precision is nonfinite or not SPD; no ridge added");
  arma::vec z(h.n_elem);
  for (arma::uword j = 0; j < z.n_elem; ++j)
    z(j) = R::rnorm(0, 1);
  arma::vec m = arma::solve(arma::trimatu(L.t()), arma::solve(arma::trimatl(L), h));
  arma::vec out = m + std::sqrt(s) * arma::solve(arma::trimatu(L.t()), z);
  if (!out.is_finite())
    stop("nonfinite coefficient draw");
  return out;
}
double noise_scale(const arma::vec &y, const arma::mat &H, const arma::mat &K,
                   const arma::vec &b, const arma::vec &B, const arma::vec &w, double a_s,
                   double b_s, int n, int p, int d) {
  // s | beta,w,y is inverse gamma; b_s is the inverse-gamma scale.
  arma::vec e = y - H * b;
  double q = arma::as_scalar(e.t() * K * e) + arma::sum(B / w);
  return inv_gamma(a_s + 0.5 * d * (n + p), b_s + 0.5 * q);
}
void gdp_restore(const arma::vec &B, double s, double alpha, double eta, int d,
                 arma::vec &lam, arma::vec &w) {
  // Joint block: lambda | b,s (w marginalized), then w | lambda,b,s.
  // No beta/scale update or storage is permitted between these two loops.
  for (arma::uword j = 0; j < w.n_elem; ++j)
    lam(j) = gamma_rate(alpha + d, eta + std::sqrt(B(j) / s));
  for (arma::uword j = 0; j < w.n_elem; ++j)
    w(j) = local_w(B(j) / s, lam(j) * lam(j));
}
void gdp_sweep(const arma::vec &y, const arma::mat &H, const arma::mat &K,
               const arma::mat &T, const arma::vec &h, arma::vec &b, arma::vec &w,
               double &s, arma::vec &lam, double alpha, double eta, double a_s,
               double b_s, int n, int p, int d) {
  // Partially collapsed GDP sweep: beta, s, then the joint (lambda,w) block.
  b = coefficient(T, h, w, s, d);
  arma::vec B = B_values(b, p, d);
  s = noise_scale(y, H, K, b, B, w, a_s, b_s, n, p, d);
  gdp_restore(B, s, alpha, eta, d, lam, w);
}
} // namespace

// Reference-test helpers use the same kernels as the fitting engine.
// [[Rcpp::export]]
Rcpp::List covariance_cpp(int n, int d, double re, double im) {
  arma::mat C = covariance(n, d, re, im);
  return List::create(_["C"] = C, _["K"] = precision(C));
}
// [[Rcpp::export]]
arma::vec local_w_cpp(int count, double chi, double psi) {
  if (count < 1)
    stop("count must be positive");
  arma::vec x(count);
  for (int i = 0; i < count; ++i)
    x(i) = local_w(chi, psi);
  return x;
}
// [[Rcpp::export]]
arma::mat gdp_restore_cpp(int count, double B, double s, double alpha, double eta,
                          int d) {
  if (count < 1 || (d != 1 && d != 2) || !std::isfinite(B) || B < 0)
    stop("invalid restore arguments");
  positive(s, "s");
  positive(alpha, "alpha");
  positive(eta, "eta");
  arma::mat out(count, 2);
  arma::vec lam(1), w(1), b(1);
  b(0) = B;
  for (int i = 0; i < count; ++i) {
    gdp_restore(b, s, alpha, eta, d, lam, w);
    out(i, 0) = lam(0);
    out(i, 1) = w(0);
  }
  return out;
}
// Independent-state one-sweep endpoint for invariant-distribution certification (p=1).
// [[Rcpp::export]]
arma::mat gdp_invariant_cpp(arma::vec y, arma::mat H, int d, double re, double im,
                            arma::mat bs, double alpha, double eta, double a_s,
                            double b_s) {
  if ((d != 1 && d != 2) || H.n_cols != static_cast<arma::uword>(d) ||
      H.n_rows != y.n_elem || H.n_rows % d ||
      bs.n_cols != static_cast<arma::uword>(d + 1) || !bs.is_finite() ||
      arma::any(bs.col(d) <= 0))
    stop("invalid invariant-test dimensions/state");
  positive(alpha, "alpha");
  positive(eta, "eta");
  positive(a_s, "a_s");
  positive(b_s, "b_s");
  arma::mat K = precision(covariance(H.n_rows / d, d, re, im)), T = H.t() * K * H;
  arma::vec h = H.t() * K * y;
  arma::mat out(bs.n_rows, d + 3);
  for (arma::uword i = 0; i < bs.n_rows; ++i) {
    arma::vec b = bs.row(i).head(d).t(), w(1), lam(1);
    double s = bs(i, d);
    gdp_restore(B_values(b, 1, d), s, alpha, eta, d, lam, w);
    gdp_sweep(y, H, K, T, h, b, w, s, lam, alpha, eta, a_s, b_s, H.n_rows / d, 1, d);
    out.submat(i, 0, i, d - 1) = b.t();
    out(i, d) = s;
    out(i, d + 1) = lam(0);
    out(i, d + 2) = w(0);
  }
  return out;
}
// [[Rcpp::export]]
Rcpp::List canonical_cpp(arma::vec y, arma::mat H, int d, std::string family, double re,
                         double im, int iterations, int warm, int thin, arma::vec beta,
                         arma::vec w, double s, double global_L, arma::vec lam,
                         double a_tau, double b_tau, double r, double delta, double alpha,
                         double eta, double a_s, double b_s) {
  if ((d != 1 && d != 2) || H.n_rows != y.n_elem || H.n_rows % d || H.n_cols % d ||
      !H.n_cols)
    stop("incompatible real-block dimensions");
  int n = H.n_rows / d, p = H.n_cols / d;
  if (iterations < 1 || warm < 0 || warm >= iterations || thin < 1)
    stop("invalid storage schedule");
  if (!y.is_finite() || !H.is_finite() || !beta.is_finite() || beta.n_elem != H.n_cols ||
      w.n_elem != static_cast<arma::uword>(p) || !w.is_finite() || arma::any(w <= 0))
    stop("invalid data or starting state");
  if (family != "student" && family != "lasso" && family != "gdp")
    stop("unknown family");
  positive(s, "s");
  positive(a_s, "a_s");
  positive(b_s, "b_s");
  positive(a_tau, "a_tau");
  positive(b_tau, "b_tau");
  positive(r, "r");
  positive(delta, "delta");
  positive(alpha, "alpha");
  positive(eta, "eta");
  positive(global_L, "L");
  if (lam.n_elem != static_cast<arma::uword>(p) || !lam.is_finite() ||
      arma::any(lam <= 0))
    stop("invalid lambda start");
  arma::mat C = covariance(n, d, re, im), K = precision(C);
  // These caches depend only on fixed H, y, rho. No w/s-dependent cache survives a sweep.
  arma::mat T = H.t() * K * H;
  arma::vec h = H.t() * K * y;
  if (!T.is_finite() || !h.is_finite())
    stop("fixed likelihood products overflowed");
  int kept = (iterations - warm) / thin;
  if (kept < 1)
    stop("schedule retains no draws");
  arma::mat bd(kept, d * p), wd(kept, p),
      ld(kept, family == "gdp" ? p : (family == "lasso" ? 1 : 0));
  arma::vec sd(kept);
  IntegerVector idx(kept);
  int row = 0;
  // Gibbs sweeps. Preserve RNG order and restore GDP local scales before storage.
  for (int iter = 1; iter <= iterations; ++iter) {
    if (iter % 256 == 0)
      checkUserInterrupt();
    if (family == "gdp") {
      gdp_sweep(y, H, K, T, h, beta, w, s, lam, alpha, eta, a_s, b_s, n, p, d);
    } else {
      // Reconstruct A using current w; coefficient covariance uses current s.
      beta = coefficient(T, h, w, s, d);
      arma::vec B = B_values(beta, p, d);
      // Student: inverse gamma local scales; lasso: GIG local scales.
      for (int j = 0; j < p; ++j) {
        w(j) = family == "student" ? inv_gamma(a_tau + 0.5 * d, b_tau + 0.5 * (B(j) / s))
                                   : local_w(B(j) / s, global_L);
      }
      if (family == "lasso")
        global_L = gamma_rate(0.5 * (d + 1) * p + r, delta + 0.5 * arma::sum(w));
      s = noise_scale(y, H, K, beta, B, w, a_s, b_s, n, p, d);
    }
    // Keep every thin-th post-warmup state, starting at warm + thin.
    if (iter > warm && (iter - warm) % thin == 0) {
      bd.row(row) = beta.t();
      wd.row(row) = w.t();
      sd(row) = s;
      idx[row] = iter;
      if (family == "gdp")
        ld.row(row) = lam.t();
      if (family == "lasso")
        ld(row, 0) = global_L;
      ++row;
    }
  }
  if (row != kept)
    stop("storage invariant failed");
  return List::create(_["beta"] = bd, _["w"] = wd, _["s"] = sd, _["hyper"] = ld,
                      _["sampleidx"] = idx,
                      _["final"] = List::create(_["beta"] = beta, _["w"] = w, _["s"] = s,
                                                _["L"] = global_L, _["lambda"] = lam));
}
