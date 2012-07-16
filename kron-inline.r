# deBoor algorithm in C
# Rcpp inline

# aim: compute kron(mat1,mat2) * y = x
# mat1 is (n1,m1)
# mat2 is (n2,m2)
# y    is (m,1)
# result:
# x    is (m,1)
# m = m1 * m2

# 2 matrices and 1 vector

# looped implementation.
# 
rm(list=ls(all=T))


library(RcppEigen)
library(inline)
source("~/Dropbox/git/Rtools/tools.r")	# uses: kron.prod()
library(rbenchmark)

incl.code <- '
#include <Eigen/Dense>

Eigen::VectorXd kronproddense(
		Eigen::MatrixXd a0,
		Eigen::MatrixXd a1,
		Eigen::VectorXd y )
{
	Eigen::VectorXd retvec = Eigen::VectorXd::Zero( a0.rows() * a1.rows() );

	//iterate over rows
	for (int row_idx0=0; row_idx0<a0.rows(); ++row_idx0) {

		int row_offset1 = row_idx0;		// how much to offset index for second matrix a1?
		row_offset1    *= a1.rows();	// tells you size of jump you have to make

		for (int row_idx1=0; row_idx1<a1.rows(); ++row_idx1) {

			//start looping over columns now
			//columns of a0:
			for (int col_idx0=0; col_idx0<a0.cols(); ++col_idx0 ) {

				int col_offset1 = col_idx0;		// same offsetting story for columns
				col_offset1    *= a1.cols();
				double factor1  = a0(row_idx0,col_idx0);	// precompute parts of the multiplication. if many dimensions, that saves a lot of access operations!

				// columns of a1:
				for (int col_idx1=0; col_idx1<a1.cols(); ++col_idx1 ) {

					// compute product at the corresponding index
					retvec( row_offset1 + row_idx1 ) += factor1 * a1(row_idx1,col_idx1) * y( col_offset1 + col_idx1 );
				}
			}
		}
	}
	return retvec;
}'

cppcode <- '
//define a type for mapping R to Eigen objects
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

//initiate objects as constant types so that changing them is impossible.
const MapMatd a0(Rcpp::as<MapMatd>(aa0));
const MapMatd a1(Rcpp::as<MapMatd>(aa1));
const MapMatd y(Rcpp::as<MapMatd>(yy));

Eigen::VectorXd result( kronproddense(a0,a1,y) );
return wrap(result);
'

# create R objects
n1 <- 20
m1 <- 20
n2 <- 30
m2 <- 30
m <- m1*m2
set.seed <- 12344

aa0 <- matrix(rnorm(n1*m1),n1,m1)
aa1 <- matrix(rnorm(n2*m2),n2,m2)
yy <- array(1,dim=c(m,1))

truth <- kronecker(aa0,aa1) %*% yy

kroncpp <- cxxfunction(signature(aa0="matrix",aa1="matrix",yy="numeric"),body=cppcode, plugin="RcppEigen",includes=incl.code)

result <- kroncpp(aa0,aa1,yy)

all.equal(as.numeric(truth),result)

library(rbenchmark)

benchmark(cpp=kroncpp(aa0,aa1,yy),R=kron.prod(y=yy,matrices=list(aa0,aa1)))


# 3 dimensions
# ============

# compute kronecker(a0,kronecker(a1,a2)) * y
# a0,a1,a2 are matrices, y is a vector of length a0.cols() * a1.cols() * a2.cols()
# return retvec is length a0.rows() * a1.rows() * a2.rows()

incl.code <- '
#include <Eigen/Dense>

Eigen::VectorXd kronproddense3(
		Eigen::MatrixXd a0,
		Eigen::MatrixXd a1,
		Eigen::MatrixXd a2,
		Eigen::VectorXd y )
{
	int m = a0.cols() * a1.cols() * a2.cols();
	int n = a0.rows() * a1.rows() * a2.rows();
	Eigen::VectorXd retvec = Eigen::VectorXd::Zero( n );

	if ( y.size() != m ) {
        throw std::runtime_error( "y and kronecker product dimensions are not comformable" );
	}

	//iterate over rows
	for (int row_idx0=0; row_idx0<a0.rows(); ++row_idx0) {

		int row_offset1 = row_idx0;		// how much to offset index for second matrix a1?
		row_offset1    *= a1.rows();	// tells you size of jump you have to make

		for (int row_idx1=0; row_idx1<a1.rows(); ++row_idx1) {
			
			int row_offset2 = row_offset1 + row_idx1;		// how much to offset index for third matrix a2?
			row_offset2    *= a2.rows();	// tells you size of jump you have to make

			for (int row_idx2=0; row_idx2<a2.rows(); ++row_idx2) {

				//start looping over columns now
				//columns of a0:
				for (int col_idx0=0; col_idx0<a0.cols(); ++col_idx0 ) {

					int col_offset1 = col_idx0;		// same offsetting story for columns
					col_offset1    *= a1.cols();
					double factor1  = a0(row_idx0,col_idx0);	// precompute parts of the multiplication. if many dimensions, that saves a lot of access operations!

					// columns of a1:
					for (int col_idx1=0; col_idx1<a1.cols(); ++col_idx1 ) {
						
						int col_offset2 = col_offset1 + col_idx1;		// same offsetting story for columns
						col_offset2    *= a2.cols();
						double factor2  = factor1 * a1(row_idx1,col_idx1);	// precompute parts of the multiplication. if many dimensions, that saves a lot of access operations!

						// columns of a2:
						for (int col_idx2=0; col_idx2<a2.cols(); ++col_idx2 ) {

							// compute product at the corresponding index
							retvec( row_offset2 + row_idx2 ) += factor2 * a2(row_idx2,col_idx2) * y( col_offset2 + col_idx2 );
						}
					}
				}
			}
		}
	}
	return retvec;
}'

cppcode <- '
//define a type for mapping R to Eigen objects
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

//initiate objects as constant types so that changing them is impossible.
const MapMatd a0(Rcpp::as<MapMatd>(aa0));
const MapMatd a1(Rcpp::as<MapMatd>(aa1));
const MapMatd a2(Rcpp::as<MapMatd>(aa2));
const MapMatd y(Rcpp::as<MapMatd>(yy));

Eigen::VectorXd result( kronproddense3(a0,a1,a2,y) );
return wrap(result);
'

# create R objects
n1 <- 10
m1 <- 10
n2 <- 30
m2 <- 30
n3 <- 9
m3 <- 9
m <- m1*m2*m3
set.seed <- 12344

aa0 <- matrix(rnorm(n1*m1),n1,m1)
aa1 <- matrix(rnorm(n2*m2),n2,m2)
aa2 <- matrix(rnorm(n3*m3),n3,m3)
yy <- array(1,dim=c(m,1))

truth <- kronecker(aa0,kronecker(aa1,aa2)) %*% yy

kroncpp <- cxxfunction(signature(aa0="matrix",aa1="matrix",aa2="matrix",yy="numeric"),body=cppcode, plugin="RcppEigen",includes=incl.code)

result <- kroncpp(aa0,aa1,aa2,yy)

result.r <- kron.prod(y=yy,matrices=list(aa2,aa1,aa0))

all.equal(as.numeric(truth),result)
all.equal(as.numeric(truth),result.r)


benchmark(cpp=kroncpp(aa0,aa1,aa2,yy),R=kron.prod(y=yy,matrices=list(aa2,aa1,aa0)),replications=10)

# 4 dimensions
# ============

# compute kronecker(a0,kronecker(a1,kronecker(a2,a3))) * y
# a0,a1,a2,a3 are matrices, y is a vector of length a0.cols() * a1.cols() * a2.cols()* a3.cols()
# return retvec is length a0.rows() * a1.rows() * a2.rows()* a3.rows()

incl.code <- '
#include <Eigen/Dense>

Eigen::VectorXd kronproddense4(
		Eigen::MatrixXd a0,
		Eigen::MatrixXd a1,
		Eigen::MatrixXd a2,
		Eigen::MatrixXd a3,
		Eigen::VectorXd y )
{
	int m = a0.cols() * a1.cols() * a2.cols() * a3.cols();
	int n = a0.rows() * a1.rows() * a2.rows() * a3.rows();
	Eigen::VectorXd retvec = Eigen::VectorXd::Zero( n );

	if ( y.size() != m ) {
        throw std::runtime_error( "y and kronecker product dimensions are not comformable" );
	}

	//iterate over rows
	for (int row_idx0=0; row_idx0<a0.rows(); ++row_idx0) {

		int row_offset1 = row_idx0;		// how much to offset index for second matrix a1?
		row_offset1    *= a1.rows();	// tells you size of jump you have to make

		for (int row_idx1=0; row_idx1<a1.rows(); ++row_idx1) {
			
			int row_offset2 = row_offset1 + row_idx1;		// how much to offset index for third matrix a2?
			row_offset2    *= a2.rows();	// tells you size of jump you have to make

			for (int row_idx2=0; row_idx2<a2.rows(); ++row_idx2) {
				
				int row_offset3 = row_offset2 + row_idx2;		// how much to offset index for third matrix a2?
				row_offset3    *= a3.rows();	// tells you size of jump you have to make

				for (int row_idx3=0; row_idx3<a3.rows(); ++row_idx3) {

					//start looping over columns now
					//columns of a0:
					for (int col_idx0=0; col_idx0<a0.cols(); ++col_idx0 ) {

						int col_offset1 = col_idx0;		// same offsetting story for columns
						col_offset1    *= a1.cols();
						double factor1  = a0(row_idx0,col_idx0);	// precompute parts of the multiplication. if many dimensions, that saves a lot of access operations!

						// columns of a1:
						for (int col_idx1=0; col_idx1<a1.cols(); ++col_idx1 ) {
							
							int col_offset2 = col_offset1 + col_idx1;		// same offsetting story for columns
							col_offset2    *= a2.cols();
							double factor2  = factor1 * a1(row_idx1,col_idx1);	// precompute parts of the multiplication. if many dimensions, that saves a lot of access operations!

							// columns of a2:
							for (int col_idx2=0; col_idx2<a2.cols(); ++col_idx2 ) {
								
								int col_offset3 = col_offset2 + col_idx2;		// same offsetting story for columns
								col_offset3    *= a3.cols();
								double factor3  = factor2 * a2(row_idx2,col_idx2);	// precompute parts of the multiplication. if many dimensions, that saves a lot of access operations!

								// columns of a3:
								for (int col_idx3=0; col_idx3<a3.cols(); ++col_idx3 ) {

									// compute product at the corresponding index
									retvec( row_offset3 + row_idx3 ) += factor3 * a3(row_idx3,col_idx3) * y( col_offset3 + col_idx3 );
								}
							}
						}
					}
				}
			}
		}
	}
	return retvec;
}'

cppcode <- '
//define a type for mapping R to Eigen objects
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

//initiate objects as constant types so that changing them is impossible.
const MapMatd a0(Rcpp::as<MapMatd>(aa0));
const MapMatd a1(Rcpp::as<MapMatd>(aa1));
const MapMatd a2(Rcpp::as<MapMatd>(aa2));
const MapMatd a3(Rcpp::as<MapMatd>(aa3));
const MapMatd y(Rcpp::as<MapMatd>(yy));

Eigen::VectorXd result( kronproddense4(a0,a1,a2,a3,y) );
return wrap(result);
'

# create R objects
n1 <- 10
m1 <- 10
n2 <- 30
m2 <- 30
n3 <- 9
m3 <- 9
n4 <- 9
m4 <- 9
m <- m1*m2*m3*m4
set.seed <- 12344

aa0 <- matrix(rnorm(n1*m1),n1,m1)
aa1 <- matrix(rnorm(n2*m2),n2,m2)
aa2 <- matrix(rnorm(n3*m3),n3,m3)
aa3 <- matrix(rnorm(n4*m4),n4,m4)
yy <- array(1,dim=c(m,1))

# truth <- kronecker(aa0,kronecker(aa1,kronecker(aa2,aa3))) %*% yy

kroncpp <- cxxfunction(signature(aa0="matrix",aa1="matrix",aa2="matrix",aa3="matrix",yy="numeric"),body=cppcode, plugin="RcppEigen",includes=incl.code)

result <- kroncpp(aa0,aa1,aa2,aa3,yy)

result.r <- kron.prod(y=yy,matrices=list(aa3,aa2,aa1,aa0))

all.equal(result.r,result)


benchmark(cpp=kroncpp(aa0,aa1,aa2,aa3,yy),R=kron.prod(y=yy,matrices=list(aa3,aa2,aa1,aa0)),replications=10)

# the R code is 24 times as fast as the c implementation in that case.



# TODO: sparse.
# now with sparse implementation. fact is that working with bases, that's going to make things fast
# ===============================


// KronProdMat is a switch. number of passed matrices decides about which function to call.
Eigen::VectorXd KronProdMat(
        std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor> > a )	// a is a collection of sparse matrices
{
	int nmats  = a.size();	// number of matrices in list a

	if ( nmats == 2 ) {
		return KronProdMat2( a );
	} else if ( nmats == 3) {
		return KronProdMat3( a );
	} else if (nmats == 4 ) {
		return KronProdMat4( a );
	} else if (nmats == 5 ) {
		return KronProdMat5( a );
	} else if (nmats > 5 ) {
        throw std::runtime_error( "KronProdMat: too many (> 5) dimensions specified." );
	}
}


