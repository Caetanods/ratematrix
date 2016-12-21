##' Log likelihood function for one evolutionary rate matrix.
##'
##' This function applies a math shortcut to avoid inverting the kronecker product of matrices. The trick is based on an equality described in Ho and Ané (2014). A similar equality also make the calculation of the determinant of the kronecker product faster. Thus, one need to provide the inverse of the phylogenetic covariance matrix as the argument 'C.prime' and the determinant of this same matrix as the argument 'det.C'.
##' @title Log likelihood for one evolutionary rate matrix fitted to the data and tree.
##' @param data the data cache.
##' @param chain the chain cache.
##' @param root numeric. The phylogenetic mean.
##' @param R matrix. The covariance matrix fitted to the tree and data.
##' @return numeric. The log likelhood.
##' @references Revell, L. J., and D. C. Collar. 2009. Phylogenetic Analysis of the Evolutionary Correlation Using Likelihood. Evolution 63:1090–1100.
##' @references Ho, L. si T., and C. Ané. 2014. A Linear-Time Algorithm for Gaussian and Non-Gaussian Trait Evolution Models. Syst Biol 63:397–408.
##' @importFrom phylolm three.point.compute
##' @noRd
logLikSingleRegime <- function(data, chain, phy, root, R){
    ## X.c are the residuals. We can think of the root as the intercept of a linear regression.
    ## This standardize the data. So we can apply the 'three.point.compute'.
    X.c <- sapply(1:ncol(data$X), function(x) data$X[,x] - root[x])
    P <- data.frame(X.c)
    comp <- three.point.compute(phy, P=P, Q=NULL)
    xiVx <- sum( comp$PP * chol2inv(chol(R)) )
    detV <- data$k * comp$logd + data$n * determinant(R)$modulus[1]
    logl <- -0.5 * ( data$k * data$n * log(2 * pi) + detV + xiVx )
    return(logl)
}
