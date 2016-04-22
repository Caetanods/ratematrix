##' Log likelihood function for one evolutionary rate matrix.
##'
##' This function applies a math shortcut to avoid inverting the kronecker product of matrices. The trick is based on an equality described in Ho and Ané (2014). A similar equality also make the calculation of the determinant of the kronecker product faster. Thus, one need to provide the inverse of the phylogenetic covariance matrix as the argument 'C.prime' and the determinant of this same matrix as the argument 'det.C'.
##' @title Log likelihood for one evolutionary rate matrix fitted to the data and tree.
##' @param X matrix. The matrix of trait data.
##' @param root numeric. The phylogenetic mean.
##' @param R matrix. The covariance matrix fitted to the tree and data.
##' @param C.prime matrix. The inverse of the phylogenetic covariance matrix. This is a fixed term.
##' @param det.C numeric. The determinant (log scale) of the phylogenetic covariance matrix. This is a fixed term.
##' @param D matrix. The design matrix. This is a fixed term, dependent on the phylogeny.
##' @param n numeric. Number of tips in the phylogeny.
##' @param r numeric. Number of traits, the dimension of the R matrix.
##' @param method string. Method to be used to calculate the likelihood. Right now only the option 'rpf' is working. More options will be implemented in the future.
##' @return numeric. The log likelhood.
##' @references Revell, L. J., and D. C. Collar. 2009. Phylogenetic Analysis of the Evolutionary Correlation Using Likelihood. Evolution 63:1090–1100.
##' @references Ho, L. si T., and C. Ané. 2014. A Linear-Time Algorithm for Gaussian and Non-Gaussian Trait Evolution Models. Syst Biol 63:397–408.
singleR.loglik <- function(X, root, R, C.prime, det.C, D, n, r){

    ## (y - Da)
    Xvec <- as.vector( as.matrix(X) )
    residuals <- D %*% root - Xvec
    ## solve( R ) -> Use the Cholesky decomposition to compute the inverse. This works with square positive-definite matrices. Should be twice as fast as the general solve.
    R.prime <- chol2inv( chol( R ) )
    ## solve( V ) -> This is the inverse of the kronecker product.
    V.prime <- kronecker(R.prime, C.prime)
    ## det( V ) -> This is the determinant of the kronecker product. NOT in log scale.
    ## det.V <- kronecker( (det(R))^n, ( det.C )^r )
    det.V <- ( determinant(R, logarithm=TRUE)$modulus[1] )^n * ( det.C )^r
    ## The log.lik function:
    ## Guide code:
    ## - ((t(y - b) %*% solve(V) %*% (y - b))/2) - determinant(V)$modulus[1]/2 - (n*r*log(2*pi)/2)
    logl <- -0.5 * ( t(residuals) %*% V.prime %*% residuals) -0.5 * det.V -0.5 * (r * n * log(2 * pi) )
    
    return(logl)
}
