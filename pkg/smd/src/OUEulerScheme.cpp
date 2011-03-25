#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

RcppExport SEXP cOUEuler ( SEXP A, SEXP B, SEXP D, SEXP x0, SEXP sqdelta, SEXP n );


SEXP cOUEuler ( SEXP A, SEXP B, SEXP D, SEXP x0, SEXP sqdelta, SEXP n ) {

  NumericVector AA(A);
  NumericMatrix BB(B), DD(D);
  int p = AA.size();
////  if(p != BB.nrow() || p != BB.ncol()) {
////    error("B has wrong dimensions.");
////  }
////  if(p != D.nrow() || p != D.ncol()) {
////    error("D has wrong dimensions.");
////  } 
  double sqd = REAL(sqdelta)[0], d = sqd*sqd, ddrift, dnoise;
  int N = REAL(n)[0];
  NumericMatrix x(N, p);
  NumericVector w(p);
  int k, i, j;
  GetRNGstate();

  for(i = 0; i < p; i++)
    x(0, i) = REAL(x0)[i];

  for(k = 1; k < N; k++) {
    for(i = 0; i < p; i++) {
      w = rnorm(p);
      ddrift = 0;
      dnoise = 0;
      for(j = 0; j < p; j++) {
        ddrift += BB(i, j) * (x(k-1, j) - AA(j));
        dnoise += DD(i, j) * w(j);
      }
      x(k, i) = x(k-1, i) + d*ddrift + sqd*dnoise;
    }
  }

  PutRNGstate();
  return(x);       
}

