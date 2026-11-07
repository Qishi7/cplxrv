#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double get_tune_cpp(double tune, double keep, int iter, double step, 
                    double exp_term, double target) {
  double a;
  if (step < pow(iter, -exp_term)) {
    a = step;
  } else {
    a = pow(iter, -exp_term);
  }
  
  double d;
  if (keep < target) {
    d = log(tune) - a;
  } else {
    d = log(tune) + a;
  }
  
  return exp(d);
}
