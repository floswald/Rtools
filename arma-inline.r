
# inline ARMA tests
rm(list=ls(all=T))
library(inline)
library(RcppArmadillo)
library(rbenchmark)


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

if ( b.is_empty() ){
	throw std::logic_error( "ArmaKronProd: not all matrices are square!" );
	return 0;
}
//Rcpp::Rcout << "matrix A" << endl;
//Rcpp::Rcout << A << endl << endl;
//Rcpp::Rcout << "vector b" << endl;
//Rcpp::Rcout << b << endl << endl;

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

# benchmark against cppmax

library(RcppEigen)
library(rbenchmark)

cppcode <- '
using namespace Eigen;
using namespace std;
//define a type for mapping R to Eigen objects
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;

//initiate objects as constant types so that changing them is impossible.
const MapMatd a0(Rcpp::as<MapMatd>(aa0));

int n = a0.rows();
Rcpp::NumericVector retvals( n );
Rcpp::IntegerVector indVec( n );
Rcpp::IntegerVector ones( n );
ones.fill( 1 );

std::ptrdiff_t j;

for (int i=0; i<n; i++) {
	retvals( i ) = a0.row(i).maxCoeff(&j);
	indVec( i ) = j;
}

indVec += ones;	// go from zero based index to one based index.

Rcpp::List list = Rcpp::List::create( _["values"] = retvals, _["indices"] = indVec );
return list;
'
cppmax <- cxxfunction( signature( aa0 = "matrix" ), body = cppcode, plugin = "RcppEigen" )

mat <- matrix(rnorm(20000),nrow=4000,ncol=50)
blim <- rep(0L,nrow(mat))

benchmark(RcppEigen = cppmax(mat), RcppArma = armafun(mat,blim),replications = 1000)


# arma shot at kron.prod
# ======================

source("~/Dropbox/git/Rtools/tools.r")	# adds kron.prod()


# what is needed
# a "list" of matrices. pass with Rcpp list?

code <- '
using namespace arma;
// we get a vector of values and a list of matrices from R
arma::vec y = Rcpp::as<arma::vec>(Ry)	;
Rcpp::List mats(Rmats);
int nmat =  mats.length();	// how many matrices
int nall = y.size();

// need to map to arma::mat type
std::vector<arma::mat> list;
for (size_t i=0; i<nmat; i++){
	list.push_back( Rcpp::as<arma::mat>( mats[ i ] ) );
	// check matrices are square
	if ( !list.at(i).is_square() ){
        throw std::logic_error( "ArmaKronProd: not all matrices are square!" );
	}
//	Rcpp::Rcout << list.at( i ) << std::endl;
}

// product for first matrix
vec y0 = y;
vec y1(nall);
y1.zeros();
int n = list.at(0).n_rows;
int m = nall/n;
uvec lhs_out = linspace<uvec>(0,n-1,n);

uvec lhs_in = lhs_out * m;
uvec rhs(n);
uvec lhs(n);

for (int i=0; i<m; i++){
	lhs = lhs_in + i;
	rhs    = lhs_out + i*n;
//	Rcpp::Rcout << "lhs = " << lhs.t() << std::endl;
//	Rcpp::Rcout << "rhs = " << rhs.t() << std::endl;
//	Rcpp::Rcout << "y1.elem( lhs ) = " << y1.elem( lhs ).t() << std::endl;
//	Rcpp::Rcout << "y0.elem( rhs ) = " << y0.elem( rhs ).t() << std::endl;
//	Rcpp::Rcout << "list.at( 0 ) = " << list.at( 0 ) << std::endl << std::endl;
	y1.elem( lhs ) = list.at( 0 ) * y0.elem( rhs ) ;
}

if (nmat > 1){
	for(int imat=1; imat < nmat; imat++){
		y0 = y1;
		n  = list.at(imat).n_rows;
		m  = nall/n;
		lhs_out.resize(n);
		lhs_out = linspace<uvec>(0,n-1,n);
		uvec lhs_in = lhs_out * m;
		rhs.resize(n);
		lhs.resize(n);
		for (int i=0; i<m; i++){
			lhs = lhs_in + i;
			rhs    = lhs_out + i*n;
			y1.elem( lhs ) = list.at( imat ) * y0.elem( rhs ) ;
		}
	}
}
return Rcpp::wrap(y1);
'




fun <- cxxfunction(signature(Ry="numeric",Rmats = "list"),code, plugin = "RcppArmadillo")
input <- list(matrix(rnorm(25),nrow=5,ncol=5),matrix(rnorm(10),nrow=5,ncol=2),matrix(rnorm(4),nrow=2,ncol=2))
# input <- list(x = seq(1, 10, by = 0.5),y=matrix(rnorm(10),nrow=5,ncol=2))

# data
x1 = matrix(rnorm(25),nrow=5,ncol=5)
x2 = matrix(rnorm(9),nrow=3,ncol=3)
x3 = matrix(rnorm(4),nrow=2,ncol=2)

y = rnorm(5*3*2)

all.equal(as.numeric(fun(y,list(x1,x2,x3))),kron.prod(y=y,matrices=list(x1,x2,x3)))

# benchmark

x1 = matrix(rnorm(2500),nrow=50,ncol=50)
x2 = matrix(rnorm(100),nrow=10,ncol=10)
x3 = matrix(rnorm(36),nrow=6,ncol=6)
x4 = matrix(rnorm(25),nrow=5,ncol=5)
x5 = matrix(rnorm(16),nrow=4,ncol=4)

y = rnorm(50 * 10 * 6 * 5 * 4)

all.equal(as.numeric(fun(y,list(x1,x2,x3,x4,x5))),kron.prod(y=y,matrices=list(x1,x2,x3,x4,x5)))

benchmark(kron.prod = kron.prod(y=y,matrices=list(x1,x2,x3,x4,x5)), arma.kron = fun(y,list(x1,x2,x3,x4,x5)), replications=10)
