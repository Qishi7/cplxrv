#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::uvec keep_elem(const arma::vec& x, const arma::ivec& to_remove) {
  std::unordered_set<double> remove_set(to_remove.begin(), to_remove.end());
  std::vector<unsigned int> indices;
  for(unsigned int i = 0; i < x.n_elem; ++i) {
    if(remove_set.find(x[i]) == remove_set.end()) {
      indices.push_back(i);
    }
  }
  return arma::uvec(indices);
}
