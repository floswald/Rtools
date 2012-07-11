
# compare performance of lm() with Eigen::sparse.solve()
# want to solve Ax = b for x, where A is a sparse matrix, x and b are dense.

# aside: i think there is no benefit from sparse vector classes when doing
# the miranda+Fackler thing. say we have
# Phi*c = y
# and Phi = kron(F,kron(E,kron(D,kron(C,kron(B,A))))), for example for 
# basis matrices. The miranda fackler
# approach is to take inverses and multiply with y, instead of forming inv(Phi):
# c = kron(f,kron(e,kron(d,kron(c,kron(b,a))))) * y, small letters are inv of big ones.
# problem: inv(A)=a is not sparse. 

# here want to solve Phi*c = y with sparse solvers in Eigen.
# Phi will eventually have to be formed outside of R because of memory problems.
# until now the tests go up to a 10000 by 9000 matrix Phi
# Phi is just one huge basis matrices rather than the actual kron of several.

rm(list=ls(all=T))
library('Rcpp')
library('RcppEigen')
library('inline')
library('benchmark')

### define Rcpp functions inline:

# include text. contains: Eigen namespaces, Rcpp::as, Dense2Sparse, FillMat, FillMatRow

incl <- '
#include <Eigen/Sparse>
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::as;
typedef Map<MatrixXd> MapMatd; 
typedef Map<VectorXd> MapVecd;
	
// utility function: Rcpp::NumericVector -> Eigen::SparseVector
// (c) Jelmer Ypma 2012
Eigen::SparseVector<double> Dense2Sparse(
    Rcpp::NumericVector vIn )
{
    Eigen::SparseVector<double> vOut( vIn.size() );
    for (int i=0;i<vIn.size();i++) {
        if (vIn( i ) != 0.0 ) {
            vOut.insert( i ) = vIn( i );
        }
    }
    return vOut;
}

// utility: fill sparse matrix column major
// adapted from Jelmer Ypma 2012
Eigen::SparseMatrix<double> FillMat( Rcpp::NumericMatrix Rmat )
{
	// fills SparseMatrix retmat with non zero elts of Rmat
	// COLUMN MAJOR
    // define the return matrix
    Eigen::SparseMatrix<double> retmat( Rmat.rows(), Rmat.cols() );
	
	// loop over each colum of Rmat filling in those elements != 0
    for (int icol=0;icol< Rmat.cols();icol++) {
	// call dense2sparse!!!
	Eigen::SparseVector<double> colvec = Dense2Sparse ( Rmat.column( icol) );
        
        // loop over non-zero elements of colvec and copy these to retmat
        for (Eigen::SparseVector<double>::InnerIterator row_it(colvec); row_it; ++row_it) {
            retmat.insert( icol, row_it.index() ) = row_it.value();
        }
    }    
	// finish filling: finalize
    retmat.finalize();
    return retmat;
}

// utility: fill sparse matrix row major
// (c) Jelmer Ypma 2012
Eigen::SparseMatrix<double> FillMatRow( Rcpp::NumericMatrix Rmat )
{
	// fills SparseMatrix retmat with non zero elts of Rmat
	// ROW MAJOR
    // define the return matrix
    Eigen::SparseMatrix<double> retmat( Rmat.rows(), Rmat.cols() );
	
	// loop over each row of Rmat filling in those elements != 0
    for (int irow=0;irow< Rmat.rows();irow++) {
	// call dense2sparse!!!
	Eigen::SparseVector<double> rowvec = Dense2Sparse ( Rmat.row( irow ) );
        
        // loop over non-zero elements of rowvec and copy these to retmat
        for (Eigen::SparseVector<double>::InnerIterator col_it(rowvec); col_it; ++col_it) {
            retmat.insert( irow, col_it.index() ) = col_it.value();
        }
    }    
	// finish filling: finalize
    retmat.finalize();
    return retmat;
}
'


# compare performance of lm() and RcppEigen Cholesky decomp
# Eigen setup. Based on RcppEigen vignette.

# Cpp code

sparseLSCpp <- '
typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMatRow;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SimplicialLDLT<SpMat> SpChol;
typedef Eigen::CholmodDecomposition<SpMat> CholMD;

int myswitch = as<int>(nmethod);

const SpMatRow A(FillMatRow(Rmat));
//Rcpp::Rcout << "A :" << A << std::endl;

const VectorXd yy(as<MapVecd>(y));
//Rcpp::Rcout << "yy :" << yy << std::endl;

int n(A.rows()), m(A.cols()), ny(yy.size());
// check if right dims for t(A)*y
eigen_assert( n == ny );


const SpMatRow At( A.adjoint());	//t(a)
//Rcpp::Rcout << "At:" << At << std::endl;

const VectorXd Aty(At * yy);	//t(a)*y
//check output
//Rcpp::Rcout << "Aty:" << Aty << std::endl;

// only Cholmod method
if (myswitch == 1){
	const SpChol Ch(At * At.adjoint());			//t(a)*a
	if (Ch.info() != Eigen::Success) {
		Rcpp::Rcout << "no success in Cholmod" << std::endl;
		return R_NilValue;
	}
	return List::create(Named("betaSpChol") = Ch.solve(Aty));
}

// only L-method
if (myswitch == 2){
	const CholMD L(At);		//Cholesky decomp of t(a)
	if (L.info() != Eigen::Success) {
		Rcpp::Rcout << "no success in L" << std::endl;
		return R_NilValue;
	}
	return List::create(Named("betaCholMD") = L.solve(Aty));
}

//both
if (myswitch == 3){
	const SpChol Ch(At * At.adjoint());			//t(a)*a
	if (Ch.info() != Eigen::Success) {
		Rcpp::Rcout << "no success in Cholmod" << std::endl;
		return R_NilValue;
	}
	const CholMD L(At);		//Cholesky decomp of t(a)
	if (L.info() != Eigen::Success) {
		Rcpp::Rcout << "no success in Cholmod" << std::endl;
		return R_NilValue;
	}
	return List::create(Named("L") = wrap(L),
						Named("betaSpChol") = Ch.solve(Aty),
						Named("betaCholMD") = L.solve(Aty));
}
'

### define cxxfunction

sparseLS <- cxxfunction(signature(Rmat="matrix",y="numeric",nmethod="integer"),sparseLSCpp,"RcppEigen",incl)

# nmethod = 1 : SimplicialLDLT solution
# nmethod = 2 : CholmodDecomposition solution
# nmethod = 3 : both solutions
nmethod <- 3

### test performance

# create data
library(FunctionApproximation)

# test params
dims    <- matrix(NA,ncol=2,nrow=4)
MB      <- array(NA,dim=c(4,1))
lmtimes <- matrix(NA,nrow=4,ncol=5)
sptimes <- matrix(NA,nrow=4,ncol=5)


datas <- c(100,5000,10000,10000)	# number of datapoints
knots <- c(50,2500,5000,9000)		# number of spline knots

for (i in 1:2){

	bspace <- new("BSpline",knots[i],degree=3,ub=10,lb=0)
	btime <- system.time(mymat  <- as.matrix(bspace$calc(seq(0,10,length=datas[i]))))
	dims[i,] <- dim(mymat)
	btime
	MB[i] <- object.size(mymat)
	y      <- rnorm(n=dim(mymat)[1])

	if (i<3) {
		lmtimes[i,]<- round(system.time(my.lm <- lm.fit(mymat,y)),3)
	} else {
		lmtimes[i,] <- "crash"
	}
	sptimes[i,] <- round(system.time(my.sp <- sparseLS(mymat,y,nmethod)),3)

	if (i<3) stopifnot( all.equal(unname(my.lm$coefficients),my.sp$betaSpChol))
}

# last measurements carry huge overhead from loop. do without loop:

i <- 3
	bspace <- new("BSpline",knots[i],degree=3,ub=10,lb=0)
	btime <- system.time(mymat  <- as.matrix(bspace$calc(seq(0,10,length=datas[i]))))
	dims[i,] <- dim(mymat)
	btime
	MB[i] <- object.size(mymat)
	y      <- rnorm(n=dim(mymat)[1])

	if (i<3) {
		lmtimes[i,]<- round(system.time(my.lm <- lm.fit(mymat,y)),3)
	} else {
		lmtimes[i,] <- "crash"
	}
	sptimes[i,] <- round(system.time(my.sp <- sparseLS(mymat,y,nmethod)),3)

i <- 4
	bspace <- new("BSpline",knots[i],degree=3,ub=10,lb=0)
	btime <- system.time(mymat  <- as.matrix(bspace$calc(seq(0,10,length=datas[i]))))
	dims[i,] <- dim(mymat)
	btime
	MB[i] <- object.size(mymat)
	y      <- rnorm(n=dim(mymat)[1])

	if (i<3) {
		lmtimes[i,]<- round(system.time(my.lm <- lm.fit(mymat,y)),3)
	} else {
		lmtimes[i,] <- "crash"
	}
	sptimes[i,] <- round(system.time(my.sp <- sparseLS(mymat,y,nmethod)),3)



outdf <- data.frame(run=1:4,num.data=datas,num.knots=knots,dim=dims,MB=round(MB/1048576,2),lm.time=lmtimes[ ,3],sp.time=sptimes[,3])
outdf

