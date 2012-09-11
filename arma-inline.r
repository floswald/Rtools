
# inline ARMA tests
rm(list=ls(all=T))
library(inline)
library(RcppArmadillo)


# stupid way to maximize over each row
# stupid: cannot use member functions max(mat,1) which would find max by row, because also want to have the index
# so: loop over rows of matrix and do operations one by one.
cppcode <- '
#include <iostream>
using namespace arma;
mat A = Rcpp::as<arma::mat>(X);	// mat R matrix
Rcpp::Rcout << "matrix A" << endl;
Rcpp::Rcout << A << endl << endl;

uword j;
uvec iy(A.n_rows);
vec onevec(A.n_rows);
vec y(A.n_rows);
rowvec tmp(A.n_cols);

for (int i=0; i<A.n_rows; i++) {
	tmp = A.row(i);
	y(i) = tmp.max( j );
	iy(i) = j + 1; // go to 1 based indices
}

Rcpp::List list = Rcpp::List::create( _["values"] = y, _["indices"] = iy );
return list;
'
armafun <- cxxfunction(signature(X="matrix"),body=cppcode,plugin='RcppArmadillo')


# add a vector of limiting indices
# convection: b in {0,1,...,ncols(mat)-1} stands for how many elements are not admissible (excluded) from the left
# i.e. if b=0, none are excluded; max over entire vector
#      if b=1, exclude first element
#      if b=4, exclude elements 1,2,3 and 4.
cppcode <- '
#include <iostream>
using namespace arma;
mat A = Rcpp::as<arma::mat>(X);	// map R matrix
uvec b = Rcpp::as<arma::uvec>(B);	// map R index vector
Rcpp::Rcout << "matrix A" << endl;
Rcpp::Rcout << A << endl << endl;
Rcpp::Rcout << "vector b" << endl;
Rcpp::Rcout << b << endl << endl;

uword j;
uvec iy(A.n_rows);
vec y(A.n_rows);
rowvec tmp(A.n_cols);
rowvec tmp2;
int m = A.n_cols - 1;

for (int i=0; i<A.n_rows; i++) {
	if ( b( i ) == 0 ) {
		// no borrowing limit. max over entire row.
		tmp = A.row(i);
		y(i) = tmp.max( j );
		iy(i) = j + 1; // go to 1 based indices
	} else {
		tmp2 = A( i , span( b(i), m ) );
		y(i) = tmp2.max( j );
		iy(i) = j + b(i) + 1; // go to 1 based indices
	}
}


Rcpp::List list = Rcpp::List::create( _["values"] = y, _["indices"] = iy );
return list;
'

armafun <- cxxfunction(signature(X="matrix",B="integer"),body=cppcode,plugin='RcppArmadillo')

armafun(matrix(rnorm(40),nrow=4,ncol=10),sample(0:9,size=4,replace=TRUE))



