##' Log likelihood function for one evolutionary rate matrix.
##'
##' This function applies a math shortcut to avoid inverting the kronecker product of matrices. The trick is based on an equality described in Ho and Ané (2014). A similar equality also make the calculation of the determinant of the kronecker product faster. Thus, one need to provide the inverse of the phylogenetic covariance matrix as the argument 'C.prime' and the determinant of this same matrix as the argument 'det.C'.
##' @title Log likelihood for one evolutionary rate matrix fitted to the data and tree.
##' @param X matrix. The matrix of trait data.
##' @param phy ape phylo. The phylogenetic tree with branch lengths.
##' @param root numeric. The phylogenetic mean.
##' @param R matrix. The covariance matrix fitted to the tree and data.
##' @param n numeric. Number of tips in the phylogeny.
##' @param r numeric. Number of traits, the dimension of the R matrix.
##' @return numeric. The log likelhood.
##' @references Revell, L. J., and D. C. Collar. 2009. Phylogenetic Analysis of the Evolutionary Correlation Using Likelihood. Evolution 63:1090–1100.
##' @references Ho, L. si T., and C. Ané. 2014. A Linear-Time Algorithm for Gaussian and Non-Gaussian Trait Evolution Models. Syst Biol 63:397–408.
##' @export
singleR.loglik <- function(X, phy, root, R, r, n){
    ## X.c are the residuals. We can think of the root as the intercept of a linear regression.
    ## This standardize the data. So we can apply the 'three.point.compute'.
    X.c <- sapply(1:ncol(X), function(x) X[,x] - root[x])
    P <- data.frame(X.c)
    comp <- three.point.compute(phy, P=P, Q=NULL)
    xiVx <- sum( comp$PP * chol2inv(chol(R)) )
    detV <- r * comp$logd + n * determinant(R)$modulus[1]
    logl <- -0.5 * ( r * n * log(2 * pi) + detV + xiVx )
    return(logl)
}

## Old implementation of the same function:
## singleR.loglik <- function(X, root, R, C.prime, det.C, D, n, r){
##     ## (y - Da)
##     Xvec <- as.vector( as.matrix(X) )
##     residuals <- D %*% root - Xvec
##     ## solve( R ) -> Use the Cholesky decomposition to compute the inverse. This works with square positive-definite matrices. Should be twice as fast as the general solve.
##     R.prime <- chol2inv( chol( R ) )
##     ## solve( V ) -> This is the inverse of the kronecker product.
##     V.prime <- kronecker(R.prime, C.prime)
##     ## det.C is in log scale. So we need to transform all the expression to the log.
##     ## This is log( ( det(R) )^n * ( det(C) )^r )
##     det.V <- ( determinant(R)$modulus[1] ) * n + ( det.C ) * r
##     ## The log.lik function:
##     ## Guide code:
##     ## - ((t(y - b) %*% solve(V) %*% (y - b))/2) - determinant(V)$modulus[1]/2 - (n*r*log(2*pi)/2)
##     logl <- -0.5 * ( ( t(residuals) %*% V.prime %*% residuals) + det.V + (r * n * log(2 * pi) ) )
##     return(logl)
## }
