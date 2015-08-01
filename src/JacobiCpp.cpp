#include <Rcpp.h>
using namespace Rcpp;

SEXP machine_double_eps(std::string value = "double.eps") // not exported.
{
    return (as<List>(Environment::base_env()[".Machine"]))[value];
}

NumericMatrix Ident(int n) // not exported.
{
    NumericMatrix I(n, n);
    for(int i = 0; i < n; i++) I(i, i) = 1.0;
    return I;
}

//' The Classical Jacobi Algorithm in Compiled Code
//'
//' Eigenvalues and optionally, eigenvectore, of a real symmetric matrix using the
//' classical Jacobi algorithm, (Jacobi, 1854)
//' @import Rcpp
//' @useDynLib JacobiEigen
//' @title The Jacobi Algorithm using Rcpp
//' @param x A real symmetric matrix
//' @param only_values A logical value: do you want eigenvalues only?
//' @param eps an error tolerance. 0.0 implies \code{.Machine$double.eps} and
//'   \code{sqrt(.Machine$double.eps)} if \code{only_values = TRUE
//'   }
//' @export JacobiCpp
//' @examples
//' V <- crossprod(matrix(1:25, 5))
//' JacobiCpp(V)
//' identical(JacobiCpp(V), JacobiR(V))
//' all.equal(JacobiCpp(V)$values, base::eigen(V)$values)
//' @return a list of two components as for \code{base::eigen}
// [[Rcpp::export]]
List JacobiCpp(NumericMatrix x, bool only_values = false, double eps = 0.0)
{
    NumericMatrix S(clone(x));
    int nr = S.nrow();
    bool vectors = !only_values;
    NumericMatrix H;

    if(vectors) {
      H = Ident(nr);
    }

    bool def = only_values & (eps == 0.0);
    double eps0 = as<double>(machine_double_eps());
    eps = eps > eps0 ? eps : eps0;  // i.e. tol. no lower than .Machine$double.eps
    if(def) eps = sqrt(eps); // only a lower accuracy is needed for eigenvalues only.

    while(true) {
	    double maxS = 0.0;
	    int i=0, j=0;
	    for(int row = 1; row < nr; row++) {  // find value & position of maximum |off-diagonal|
	        for(int col = 0; col < row; col++) {
	        	double val = fabs(S(row, col));
		        if(maxS < val) {
		          maxS = val;
		          i = row;
		          j = col;
	        	}
	       }
	    }
	    if(maxS < eps) break;

	    NumericVector Si = S(_, i), Sj = S(_, j);

	    double theta = 0.5*atan2(2.0*Si(j), Sj(j) - Si(i));
	    double s = sin(theta), c = cos(theta);

	    S(i, _) = S(_, i) = c*Si - s*Sj;
	    S(j, _) = S(_, j) = s*Si + c*Sj;
	    S(i, j) = S(j, i) = 0.0;
	    S(i, i) = c*c*Si(i) - 2.0*s*c*Si(j) + s*s*Sj(j);
	    S(j, j) = s*s*Si(i) + 2.0*s*c*Si(j) + c*c*Sj(j);

      if(vectors) {
	        NumericVector Hi = H(_, i);
	        H(_, i) = c*Hi - s*H(_, j);
	        H(_, j) = s*Hi + c*H(_, j);
      }
    }
    if(vectors) {
      return List::create(_["values"] = diag(S),
                          _["vectors"] = H);
    } else {
      return List::create(_["values"] = diag(S),
                          _["vectors"] = R_NilValue);
    }
}
