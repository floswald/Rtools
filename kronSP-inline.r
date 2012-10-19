
# test file for kronecker products on sparse matrices
# R inline

library(RcppEigen)
library(inline)
library(rbenchmark)
library(Matrix)
source("~/Dropbox/git/Rtools/tools.r")


# test with 4 dimensions. 
# ========================

# make R data
a = matrix(c(1,2,3,4,5,0,0,6,0,0,0,7,8,0,9,0,0,0,10,0),nrow=5,ncol=4,byrow=T)
b <- matrix(c(1,2,0,0,0,3),nrow=2,ncol=3,byrow=T)
c <- matrix(c(0,1,0,3,4,0),nrow=3,ncol=2,byrow=T)
d <- matrix(c(0,1,2,0,2,0,0,1),nrow=2,ncol=4,byrow=T)
aa = as(a,"dgCMatrix")
bb = as(b,"dgCMatrix")
cc = as(c,"dgCMatrix")
dd = as(d,"dgCMatrix")

y <- as.numeric(1:(ncol(aa)*ncol(bb)*ncol(cc)*ncol(dd)))

# R kronecker product
r.result <- kronecker(aa,kronecker(bb,kronecker(cc,dd))) %*% y

# C++ code

include <- '
// KronProdSPMat4Print
// computes kronecker(a0,kronecker(a1,kronecker(a2,a3))) * y, a0, a1, a2, a3 SPARSE matrices!
// dim(a0) = (n0,m0)
// dim(a1) = (n1,m1)
// dim(a2) = (n2,m2)
// dim(a3) = (n3,m3)
// length(y) = m0 * m1 * m2 * m3
using namespace Eigen;
using namespace std;
Eigen::VectorXd KronProdSPMat4(
		Eigen::SparseMatrix<double, Eigen::RowMajor> a0,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a1,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a2,
		Eigen::SparseMatrix<double, Eigen::RowMajor> a3,
		Eigen::VectorXd y) {

	Eigen::VectorXd retvec;
	retvec.setZero( a0.rows() * a1.rows() * a2.rows() * a3.rows() );
	if ( y.rows() != a0.cols() * a1.cols() * a2.cols()  * a3.cols()) {
		cout << "KronProdMat4 error: y and matrices not conformable" << endl;
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

kroncpp <- cxxfunction( signature(aa0="dgCMatrix",aa1="dgCMatrix",aa2="dgCMatrix",aa3="dgCMatrix",yy="numeric"), body=cpp.code, includes = include, plugin="RcppEigen")

# call c functino

c.result <- kroncpp( aa,bb,cc,dd,y )

# check output
all.equal(as.numeric(r.result),c.result)

# five dimensions
# =================

e <- matrix(c(0,1,0,0,2,0,1,0),nrow=4,ncol=2,byrow=TRUE)
ee <- as(e,"dgCMatrix")

y <- as.numeric(1:(ncol(aa)*ncol(bb)*ncol(cc)*ncol(dd) * ncol(ee)))

# R kronecker product
r.result <- kronecker(aa,kronecker(bb,kronecker(cc,kronecker(dd,ee)))) %*% y

# cpp codes
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

# call c functino

c.result <- kroncpp( aa,bb,cc,dd,ee,y )

# check output
all.equal(as.numeric(r.result),c.result)






# benchmark against kron.prod. new data with all square mats needed
# benchmark 4 dimensions
# -----------------------------------------------------------------

# mapping kron.prod to kroncpp: 
# kron.prod: fastest index first
# kroncpp: fastest index last

a = matrix(c(1,2,3,4,5,0,0,6,0,0,0,7,8,0,9,0),nrow=4,ncol=4,byrow=T)
b <- matrix(c(1,2,0,0,0,3,0,0,2),nrow=3,ncol=3,byrow=T)
c <- matrix(c(0,1,0,3),nrow=2,ncol=2,byrow=T)
d <- matrix(rep(c(0,1,2,0,2,0,0,1),2),nrow=4,ncol=4,byrow=T)
aa = as(a,"dgCMatrix")
bb = as(b,"dgCMatrix")
cc = as(c,"dgCMatrix")
dd = as(d,"dgCMatrix")

y <- as.numeric(1:(ncol(aa)*ncol(bb)*ncol(cc)*ncol(dd)))

# R kronecker product
r.result <- kronecker(aa,kronecker(bb,kronecker(cc,dd))) %*% y
c.result <- kroncpp( aa,bb,cc,dd,y )
all.equal(as.numeric(r.result),c.result)

# bring on kron.prod.r
rr.result <- kron.prod( y = y, matrices = list(d,c,b,a) )
all.equal( as.numeric(r.result), rr.result)

benchmark(R = kron.prod( y = y, matrices = list(d,c,b,a) ), cpp = kroncpp( aa,bb,cc,dd,y ), replications = 100)


# benchmark 5 dimensions
# -----------------------------------------------------------------

a = matrix(c(1,2,3,4,5,0,0,6,0,0,0,7,8,0,9,0),nrow=4,ncol=4,byrow=T)
b <- matrix(c(1,2,0,0,0,3,0,0,2),nrow=3,ncol=3,byrow=T)
c <- matrix(c(0,1,0,3),nrow=2,ncol=2,byrow=T)
d <- matrix(rep(c(0,1,2,0),9),nrow=6,ncol=6,byrow=T)
e <- matrix(c(rep(c(0,1,2,0,2),2),rep(c(3,1,2,0,0),3)),nrow=5,ncol=5,byrow=T)
aa = as(a,"dgCMatrix")
bb = as(b,"dgCMatrix")
cc = as(c,"dgCMatrix")
dd = as(d,"dgCMatrix")
ee = as(e,"dgCMatrix")

y <- as.numeric(1:(ncol(aa)*ncol(bb)*ncol(cc)*ncol(dd)*ncol(ee)))

# R kronecker product
r.result <- kronecker(aa,kronecker(bb,kronecker(cc,kronecker(dd,ee)))) %*% y
c.result <- kroncpp( aa,bb,cc,dd,ee,y )
all.equal(as.numeric(r.result),c.result)

# bring on kron.prod.r
rr.result <- kron.prod( y = y, matrices = list(e,d,c,b,a) )
all.equal( as.numeric(r.result), rr.result)

benchmark(R = kron.prod( y = y, matrices = list(e,d,c,b,a) ), cpp = kroncpp( aa,bb,cc,dd,ee,y ), replications = 500)

