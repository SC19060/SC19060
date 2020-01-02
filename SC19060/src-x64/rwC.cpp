#include <Rcpp.h>
using namespace Rcpp;

//' @title A Metropolis sampler using Rcpp
//' @description A Metropolis sampler using Rcpp
//' @param sigma variance of the proposal distribution
//' @param x0 the initial value
//' @param N the number of samples
//' @return a randpm sample of Cauchy distribution and the number of candidate points rejected
//' @examples
//' \dontrun{
//' dir_cpp <- 'D:/Rstudio/Rcpp/'
//' source(paste0(dir_cpp,'rwR.R'))
//' sourceCpp(paste0(dir_cpp,'rwC.cpp'))
//' N <- 2000
//' sigma <- .5
//' x0 <- 25
//' rwC <- rwC(sigma, x0, N)
//' print(rwC[N+1]/N)
//' }
//' @export
// [[Rcpp::export]]
NumericVector rwC(double sigma, double x0, int N) {
  NumericVector x(N+1);
  x[0] = x0;
  NumericVector u = as<NumericVector>(runif(N));
  int k = 0;
  for (int i = 1; i < N; i++) {
    double y = as<double>(rnorm(1, x[i-1], sigma));
    if (u[i] <= (0.5*exp(-abs(y)) / (0.5*exp(-abs(x[i-1])))))
      x[i] = y ;
    else {
      x[i] = x[i-1];
      k = k + 1;
    } 
  }
  x[N] = k;
  return(x);
}