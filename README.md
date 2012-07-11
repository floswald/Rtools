Rtools
======

some tools for R in tools.r
some random examples in howto.r (mostly things from my own work, but also some nice ones from R-bloggers.)


tools.r content
---------------

1. grid creator: function to place points in an interval according to various rules.
2. linear mapping: function to map an arbitrary interval into [0,up] 
3. nonlinear mapping: function to map an arbitrary interval into [0,up] with some nonlinear rules
4. construct a Normal copula [tlamadon](https://github.com/tlamadon/Utils)
5. repmat [tlamadon](https://github.com/tlamadon/Utils)
6. spread [tlamadon](https://github.com/tlamadon/Utils)
7. kron.prod: an R-implementation of a very efficient fortran routine to compute kronecker products
8. knot.selector: a function that selects spline knots after a certain rule.


howto.r content
---------------

examples for

+ R base plot: random walks
+ ggplot2: random walks
+ grid.arrange: arrange ggplot2 plots in a panel
+ viewPort: arrange ggplot2 plots in a panel (more flexible). 
+ scrape data from HTML tables
+ data.table usage
+ manipulate strings
+ manipulate time series objects


eigen-sparse-lm content
-----------------------

worked example from RcppEigen vignette. compares performance of `lm()` with sparse matrix solver. 


kron-inline content
-------------------

inline test file for various kronecker product routines. uses Eigen sparse/dense.

