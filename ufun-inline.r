

# test file for cpp functions
library(RcppEigen)
library(inline)


# utility function with labor supply
# ================

incl.code <- '
#include <iostream>
#include <Eigen/Core>
#include <vector>
using namespace std;
using namespace Eigen;

// implementation if house size is a vector
Eigen::ArrayXXd utilfun(
		Eigen::MatrixXd a0,
		Eigen::ArrayXd size,
		int age,
		Eigen::MatrixXd vecparams,
		NumericVector dparams,
		int owner)
{
	//declarations
	int m = a0.cols();
	int n = a0.rows();
	int nd = dparams.size();
	Eigen::ArrayXd tmpvec(m);
	double diff;
	double g;
	double u;
	double grad, hess;
	Eigen::ArrayXXd ret;
	Eigen::ArrayXXd a0_arr = a0.array();
	ret.setZero(n,m);
	
	const double cutoff = dparams( nd-1 );
	Eigen::ArrayXd sizefac = size.pow(dparams(1));
	Eigen::ArrayXd sizefac2 =  dparams(4) + dparams(5) / vecparams.col(2).coeff(age-1) * size ;
	if (size.size() != n){
        throw std::runtime_error( "length of house sizes is not equal length of cash" );
	}

	for (int i=0; i<n; i++){
		if((a0_arr.row(i) < cutoff).any()){	// if any entry in row i is less than cutoff
			for (int j=0; j<m; j++){
				if(a0_arr(i,j) < cutoff){ // compute approximated utility at all entries of row i below cutoff
					g        = pow(cutoff,dparams(0)) * sizefac(i);
					diff     = a0_arr(i,j) - cutoff;
					grad     = vecparams.col(0).coeff(age-1) * dparams(3) * g / cutoff;
					hess     = dparams(2) * grad / cutoff;
					ret(i,j) = vecparams.col(1).coeff(age-1) * (g - 1) + (grad * diff) + 0.5 * hess * pow(diff,2);
					if (owner==1) {
						ret(i,j) += sizefac2(i);
					}
				} else {
					g = pow(static_cast<double>(a0_arr(i,j)),dparams(0)) * sizefac(i);	//pow(double x, double expo)
					ret(i,j) = vecparams.col(1).coeff(age-1) * (g - 1);
					if (owner==1) {
						ret(i,j) += sizefac2(i);
					}
				}
			}
		} else {
			tmpvec = a0_arr.row(i).pow(dparams(0)) * sizefac(i);
			ret.row(i) = vecparams.col(1).coeff(age-1) * (tmpvec - 1.0); // can substract a scalar from an array (not a matrix)
			if (owner==1) {
				ret.row(i) += sizefac2(i);
			}
		}
	}
	return ret;
}
'
cppcode <- '
Eigen::MatrixXd Res = Rcpp::as<Eigen::Map<Eigen::MatrixXd> > >(Res_);	// map R matrix to Eigen type
int s               = Rcpp::as<int>(s_);	// house size
int age             = Rcpp::as<int>(age_);	// age
int own             = Rcpp::as<int>(age_);	// owner yes/no
// now get parameters out of list par
cparams <- list(xi=up$xi,alpha=alpha,theta=theta,cutoff=cutoff,fsize=up$fsize)

Rcpp::List par( par_ ) ;
Rcpp::List xi( par["xi"] );
	
dparams_R   <-  with(up,c(xi[[3]],xi[[4]],xi[[5]],alpha,theta[[1]],theta[[2]],cutoff))
 		             dparams( 0      1       2      3      4          5         6     )
vecparams_R <-  with(up,cbind(xi[[1]],xi[[2]],fsize))
 				 vecparams(   0        1        2   )


//define a type for mapping R to Eigen objects
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;
typedef Eigen::Map<Eigen::ArrayXd> MapArrd;

//initiate objects as constant types so that changing them is impossible.
const MapMatd a0(Rcpp::as<MapMatd>(aa0));
const MapMatd vecparams(Rcpp::as<MapMatd>(vecparams_R));
const MapArrd size(Rcpp::as<MapArrd>(size_R));
const int age(Rcpp::as<int>(age_R));
const int owner(Rcpp::as<int>(owner_R));
const Rcpp::NumericVector dparams(dparams_R);

Eigen::ArrayXXd result( utilfun(a0,size,age,vecparams,dparams,owner) );

return wrap(result);
'

cutil <- cxxfunction(signature(aa0="matrix",size_R="numeric",vecparams_R = "matrix", dparams_R = "numeric", age_R="integer", owner_R="integer"),body=cppcode, plugin="RcppEigen",includes=incl.code)








incl.coded <- '
// implementation if house size is a double
Eigen::ArrayXXd utilfun(
		Eigen::MatrixXd a0,
		double size,
		int age,
		Eigen::MatrixXd vecparams,
		NumericVector dparams,
		int owner)
{
	//declarations
	int m = a0.cols();
	int n = a0.rows();
	int nd = dparams.size();
	Eigen::ArrayXd tmpvec(m);
	double diff;
	double g;
	double u;
	double grad, hess;
	Eigen::ArrayXXd ret;
	Eigen::ArrayXXd a0_arr = a0.array();
	ret.setZero(n,m);
	
	const double cutoff = dparams( nd-1 );
	double sizefac = pow(size,dparams(1));
	double sizefac2 = dparams(4) + dparams(5) * size / vecparams.col(2).coeff(age-1);

	for (int i=0; i<n; i++){
		if((a0_arr.row(i) < cutoff).any()){	// if any entry in row i is less than cutoff
			for (int j=0; j<m; j++){
				if(a0_arr(i,j) < cutoff){ // compute approximated utility at all entries of row i below cutoff
					g        = pow(cutoff,dparams(0)) * sizefac;
					diff     = a0_arr(i,j) - cutoff;
					grad     = vecparams.col(0).coeff(age-1) * dparams(3) * g / cutoff;
					hess     = dparams(2) * grad / cutoff;
					ret(i,j) = vecparams.col(1).coeff(age-1) * (g - 1) + (grad * diff) + 0.5 * hess * pow(diff,2);
					if (owner==1) {
						ret(i,j) += sizefac2 ;
					}
				} else {
					g = pow(static_cast<double>(a0_arr(i,j)),dparams(0)) * sizefac;	//pow(double x, double expo)
					ret(i,j) = vecparams.col(1).coeff(age-1) * (g - 1);
					if (owner==1) {
						ret(i,j) += sizefac2 ;
					}
				}
			}
		} else {
			tmpvec = a0_arr.row(i).pow(dparams(0)) * sizefac;
			ret.row(i) = vecparams.col(1).coeff(age-1) * (tmpvec - 1.0); // can substract a scalar from an array (not a matrix)
			if (owner==1) {
				ret.row(i) += sizefac2;
			}
		}
	}
	return ret;
}
'


cppcode_d <- '
//define a type for mapping R to Eigen objects
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

//initiate objects as constant types so that changing them is impossible.
const MapMatd a0(Rcpp::as<MapMatd>(aa0));
const MapMatd vecparams(Rcpp::as<MapMatd>(vecparams_R));
const double sized(Rcpp::as<double>(sized_R));
const int age(Rcpp::as<int>(age_R));
const int owner(Rcpp::as<int>(owner_R));
const Rcpp::NumericVector dparams(dparams_R);

Eigen::ArrayXXd result( utilfun(a0,sized,age,vecparams,dparams,owner) );

return wrap(result);
'

# create C++ functions:

# house size is a vector

# if house size is a scalar:
cutild <- cxxfunction(signature(aa0="matrix",sized_R="numeric",vecparams_R = "matrix", dparams_R = "numeric", age_R="integer", owner_R="integer"),body=cppcode_d, plugin="RcppEigen",includes=incl.coded)

# KronProdSPMat4: c++ kronecker product for 4 dimensions. geared towards spline basis matrices.
# =============================================================================================

include <- '
using namespace Eigen;
using namespace std;
Eigen::VectorXd KronProdSPMat4(
		Eigen::SparseMatrix<double, Eigen::RowMajor> a0,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a1,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a2,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a3,
		Eigen::VectorXd y) {

	Eigen::VectorXd retvec;
	retvec.setZero( a0.rows() * a1.rows() * a2.rows() * a3.rows()  );
	if ( y.rows() != a0.cols() * a1.cols() * a2.cols()  * a3.cols()  ) {
		cout << "KronProdMat5 error: y and matrices not conformable" << endl;
	}

	//loop rows a0
	for (int row_idx0=0; row_idx0<a0.outerSize(); ++row_idx0) {
		int row_offset1 = row_idx0;
		row_offset1    *= a1.rows();

		// loop rows a1
		for (int row_idx1=0; row_idx1<a1.outerSize(); ++row_idx1) {
			int row_offset2 = row_offset1 + row_idx1;
			row_offset2    *= a2.rows();

			// loop rows a2
			for (int row_idx2=0; row_idx2<a2.outerSize(); ++row_idx2) {
				int row_offset3 = row_offset2 + row_idx2;
				row_offset3    *= a3.rows();

				// loop rows a3
				for (int row_idx3=0; row_idx3<a3.outerSize(); ++row_idx3) {

					// loop cols a0 (non-zero elements only)
					for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it0(a0,row_idx0); it0; ++it0) {
						int col_offset1 = it0.index();
						col_offset1    *= a1.innerSize();
						double factor1 = it0.value();

						// loop cols a1
						for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it1(a1,row_idx1); it1; ++it1) {
							int col_offset2 = col_offset1 + it1.index();
							col_offset2    *= a2.innerSize();
							double factor2  = factor1 * it1.value();

							//loop cols a2
							for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it2(a2,row_idx2); it2; ++it2) {
								int col_offset3 = col_offset2 + it2.index();
								col_offset3    *= a3.innerSize();
								double factor3  = factor2 * it2.value();

								for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it3(a3,row_idx3); it3; ++it3){
										retvec( row_offset3 + row_idx3 ) += factor3 * it3.value() * y( col_offset3 + it3.index() );
								}
							}
						}
					}
				}
			}
		}
	}
	return retvec;
}
'

cpp.code <- '
using Eigen::MappedSparseMatrix;
using Eigen::SparseMatrix;
using Eigen::Map;
typedef Map<Eigen::VectorXd> MapVecd;
const MappedSparseMatrix<double> a0(as<MappedSparseMatrix<double> >(aa0));
const MappedSparseMatrix<double> a1(as<MappedSparseMatrix<double> >(aa1));
const MappedSparseMatrix<double> a2(as<MappedSparseMatrix<double> >(aa2));
const MappedSparseMatrix<double> a3(as<MappedSparseMatrix<double> >(aa3));
const MapVecd y(as<MapVecd>(yy));

VectorXd result( KronProdSPMat4( a0,a1,a2,a3,y ) );
return wrap( result );
'

# create c function

kroncpp4 <- cxxfunction( signature(aa0="dgCMatrix",aa1="dgCMatrix",aa2="dgCMatrix",aa3="dgCMatrix",yy="numeric"), body=cpp.code, includes = include, plugin="RcppEigen")



# KronProdSPMat5: c++ kronecker product for 5 dimensions. geared towards spline basis matrices.
# =============================================================================================

include <- '
using namespace Eigen;
using namespace std;
Eigen::VectorXd KronProdSPMat5(
		Eigen::SparseMatrix<double, Eigen::RowMajor> a0,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a1,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a2,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a3,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a4,
		Eigen::VectorXd y) {

	Eigen::VectorXd retvec;
	retvec.setZero( a0.rows() * a1.rows() * a2.rows() * a3.rows() * a4.rows() );
	if ( y.rows() != a0.cols() * a1.cols() * a2.cols()  * a3.cols() * a4.cols() ) {
		cout << "KronProdMat5 error: y and matrices not conformable" << endl;
	}

	//loop rows a0
	for (int row_idx0=0; row_idx0<a0.outerSize(); ++row_idx0) {
		int row_offset1 = row_idx0;
		row_offset1    *= a1.rows();

		// loop rows a1
		for (int row_idx1=0; row_idx1<a1.outerSize(); ++row_idx1) {
			int row_offset2 = row_offset1 + row_idx1;
			row_offset2    *= a2.rows();

			// loop rows a2
			for (int row_idx2=0; row_idx2<a2.outerSize(); ++row_idx2) {
				int row_offset3 = row_offset2 + row_idx2;
				row_offset3    *= a3.rows();

				// loop rows a3
				for (int row_idx3=0; row_idx3<a3.outerSize(); ++row_idx3) {
					int row_offset4 = row_offset3 + row_idx3;
					row_offset4    *= a4.rows();

					//loop rows a4
					for (int row_idx4=0; row_idx4<a4.outerSize(); ++row_idx4) {

						// loop cols a0 (non-zero elements only)
						for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it0(a0,row_idx0); it0; ++it0) {
							int col_offset1 = it0.index();
							col_offset1    *= a1.innerSize();
							double factor1 = it0.value();

							// loop cols a1
							for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it1(a1,row_idx1); it1; ++it1) {
								int col_offset2 = col_offset1 + it1.index();
								col_offset2    *= a2.innerSize();
								double factor2  = factor1 * it1.value();

								//loop cols a2
								for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it2(a2,row_idx2); it2; ++it2) {
									int col_offset3 = col_offset2 + it2.index();
									col_offset3    *= a3.innerSize();
									double factor3  = factor2 * it2.value();

									for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it3(a3,row_idx3); it3; ++it3){
										int col_offset4 = col_offset3 + it3.index();
										col_offset4    *= a4.innerSize();
										double factor4  = factor3 * it3.value();

										for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it4(a4,row_idx4); it4; ++it4){
											retvec( row_offset4 + row_idx4 ) += factor4 * it4.value() * y( col_offset4 + it4.index() );
										}
									}
								}
							}

						}
					}
				}
			}
		}
	}
	return retvec;
}
'

cpp.code <- '
using Eigen::MappedSparseMatrix;
using Eigen::SparseMatrix;
using Eigen::Map;
typedef Map<Eigen::VectorXd> MapVecd;
const MappedSparseMatrix<double> a0(as<MappedSparseMatrix<double> >(aa0));
const MappedSparseMatrix<double> a1(as<MappedSparseMatrix<double> >(aa1));
const MappedSparseMatrix<double> a2(as<MappedSparseMatrix<double> >(aa2));
const MappedSparseMatrix<double> a3(as<MappedSparseMatrix<double> >(aa3));
const MappedSparseMatrix<double> a4(as<MappedSparseMatrix<double> >(aa4));
const MapVecd y(as<MapVecd>(yy));

VectorXd result( KronProdSPMat5( a0,a1,a2,a3,a4,y ) );
return wrap( result );
'

# create c function

kroncpp <- cxxfunction( signature(aa0="dgCMatrix",aa1="dgCMatrix",aa2="dgCMatrix",aa3="dgCMatrix",aa4="dgCMatrix",yy="numeric"), body=cpp.code, includes = include, plugin="RcppEigen")


# KronProdSPMat5MEM: c++ kronecker product for 5 dimensions where 4 dimensions are stored in memory (because only dim 5 changes)
# =============================================================================================

include <- '
using namespace Eigen;
using namespace std;
Eigen::VectorXd KronProdSPMat5_makeMEM(
		Eigen::SparseMatrix<double, Eigen::RowMajor> a0,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a1,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a2,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a3 {

		Eigen::MatrixXd retmat( a0.rows() * a1.rows() * a2.rows() * a3.rows(), a0.cols() * a1.cols() * a2.cols()  * a3.cols() );


		//loop rows a0
		for (int row_idx0=0; row_idx0<a0.outerSize(); ++row_idx0) {
			int row_offset1 = row_idx0;
			row_offset1    *= a1.rows();

			// loop rows a1
			for (int row_idx1=0; row_idx1<a1.outerSize(); ++row_idx1) {
				int row_offset2 = row_offset1 + row_idx1;
				row_offset2    *= a2.rows();

				// loop rows a2
				for (int row_idx2=0; row_idx2<a2.outerSize(); ++row_idx2) {
					int row_offset3 = row_offset2 + row_idx2;
					row_offset3    *= a3.rows();

					// loop rows a3
					for (int row_idx3=0; row_idx3<a3.outerSize(); ++row_idx3) {
						int row_offset4 = row_offset3 + row_idx3;
						row_offset4    *= a4.rows();

						// loop cols a0 (non-zero elements only)
						for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it0(a0,row_idx0); it0; ++it0) {
							int col_offset1 = it0.index();
							col_offset1    *= a1.innerSize();
							double factor1 = it0.value();

							// loop cols a1
							for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it1(a1,row_idx1); it1; ++it1) {
								int col_offset2 = col_offset1 + it1.index();
								col_offset2    *= a2.innerSize();
								double factor2  = factor1 * it1.value();

								//loop cols a2
								for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it2(a2,row_idx2); it2; ++it2) {
									int col_offset3 = col_offset2 + it2.index();
									col_offset3    *= a3.innerSize();
									double factor3  = factor2 * it2.value();

									for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it3(a3,row_idx3); it3; ++it3){
										int col_offset4 = col_offset3 + it3.index();
										double factor4  = factor3 * it3.value();

										retmat( row_offset4, col_offset4 ) = factor4;

										}
									}
								}
							}
						}
					}
				}
			}
		return retmat;
		}
'

cpp.code <- '
using Eigen::MappedSparseMatrix;
using Eigen::SparseMatrix;
using Eigen::Map;
typedef Map<Eigen::VectorXd> MapVecd;
const MappedSparseMatrix<double> a0(as<MappedSparseMatrix<double> >(aa0));
const MappedSparseMatrix<double> a1(as<MappedSparseMatrix<double> >(aa1));
const MappedSparseMatrix<double> a2(as<MappedSparseMatrix<double> >(aa2));
const MappedSparseMatrix<double> a3(as<MappedSparseMatrix<double> >(aa3));

MatrixXd result( KronProdSPMat5_makeMEM( a0,a1,a2,a3 ) );
return wrap( result );
'

# create c function

# kroncpp_makeMEM <- cxxfunction( signature(aa0="dgCMatrix",aa1="dgCMatrix",aa2="dgCMatrix",aa3="dgCMatrix"), body=cpp.code, includes = include, plugin="RcppEigen")


# and then use this precomputed matrix:

include <- '
using namespace Eigen;
using namespace std;
Eigen::VectorXd KronProdSPMat5_useMEM(
		Eigen::MatrixXdMatrix a_MEM,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a4,
		Eigen::VectorXd y) {

	int n = a_MEM.rows();
	int m = a_MEM.cols();
	Eigen::VectorXd retvec;
	retvec.setZero( n * a4.rows() );

	// loop over rows of memory matrix
	for (int rows_a=0; rows_a<n; ++rows_a) {

		// loop over indices of a4

		for (int row_idx4=0; row_idx4<a4.outerSize(); ++row_idx4) {

			// loop over cols of memory matrix

			for (int cols_a=0; cols_a<m; ++cols_a){

				// loop over non-zero cols of a4

				for (Eigen::SparseMatrix<double,RowMajor>::InnerIterator it4(a4,row_idx4); it4; ++it4){
					int col_offset4 =
					col_offset3    *= a3.innerSize();
					retvec( rows_a + row_idx4 ) += a_MEM(rows_a, cols_a) * it4.value() * y( cols_a + it4.index() );
				}
			}
		}
	}
	return retvec;
}
'

cpp.code <- '
using Eigen::MappedSparseMatrix;
using Eigen::SparseMatrix;
using Eigen::Map;
typedef Map<Eigen::VectorXd> MapVecd;
const MappedSparseMatrix<double> a4(as<MappedSparseMatrix<double> >(aa4));
const MapMatd a_MEM(Rcpp::as<MapMatd>(aa_MEM));
const MapVecd y(as<MapVecd>(yy));

VectorXd result( KronProdSPMat5_useMEM( a_MEM, a4 ) );
return wrap( result );
'

# create c function

# kroncpp_useMEM <- cxxfunction( signature(aa_MEM="matrix",aa4="dgCMatrix",yy="numeric"), body=cpp.code, includes = include, plugin="RcppEigen")



# A c++ alternative to 
#     owner[,maxval.stay := apply(Wo,1,max)]
#     owner[,save.stay   := apply(Wo,1,which.max)]
# is a c++ function that takes a matrix Wo and returns a list with maximum value and column index for each row of Wo 

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


# here's the armadillo implementation that takes a vector of borrowing limits
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
# armamax <- cxxfunction(signature(X="matrix",B="integer"),body=cppcode,plugin='RcppArmadillo')

