##' Log likelihood function for two or more evolutionary rate matrices.
##'
##' Internal function to be used by the MCMC. This function uses the 'rpf' method to calculate the approximation of the inverse of matrices as impplemented in mvMORPH.
##' @title Log likelihood for multiple matrices.
##' @param data cache data.
##' @param chain cache chain
##' @param root phylogenetic mean
##' @param R list. List of the R matrices.
##' @return numeric. The log likelhood.
##' @references Revell, L. J., and D. C. Collar. 2009. Phylogenetic Analysis of the Evolutionary Correlation Using Likelihood. Evolution 63:1090–1100.
##' @noRd
logLikMultRegime <- function(data, chain, root, R){
    ## Try to keep the same parameters for the likelihood function. Make it easier for the MCMC to run as it is now.
    ## Default method is "rpf"
    C.m <- data$C.m
    X <- data$X
    D <- data$D
    n <- data$n
    r <- data$k
    p <- data$p
    V.m <- lapply(1:p, function(x) kronecker(R[[x]], C.m[[x]]) )
    V <- Reduce("+", V.m)
    ntot <- n*r
    ## This method will not accept matrices with elements equal to zero.
    ## Unclear if this would be enough to break the whole MCMC.
    ## Maybe check if elements are zero prior to calculating the loglik.
    ## But this might make the process slower. Maybe not slower than using always the
    ##     full inverse.
    ms <- 0 ## Here we do not allow for measurement error. Need to implement in the future.
    ## k <- r
    ## This is the relevant part. The calculation using C.
    error <- NULL
    cholres <- .Call(mvMORPH:::Chol_RPF, V, D, X, as.integer(r), as.integer(ntot), mserr=error, ismserr=as.integer(ms))
    det <- cholres[[2]]
    ## Adjusting the format of the data.
    Xvec <- as.vector(as.matrix(X))
    residus <- D %*% root - Xvec
    quad <- .Call(mvMORPH:::Chol_RPF_quadprod, cholres[[1]], residus, as.integer(ntot))
    logl <- -.5*quad-.5*as.numeric(det)-.5*(ntot*log(2*pi))
    return(logl)
}

## This is the old implementation. In case one need to backtrack to it.
## Snapshot taken July 6, 2016
## logLikMultRegime <- function(data, chain, root, R){
##     ## Default method is "rpf"
##     C.m <- data$C.m
##     X <- data$X
##     D <- data$D
##     n <- data$n
##     r <- data$k
##     p <- data$p
##     V.m <- lapply(1:p, function(x) kronecker(R[[x]], C.m[[x]]) )
##     V <- Reduce("+", V.m)
##     ntot <- n*r
##     ## This method will not accept matrices with elements equal to zero.
##     ## Unclear if this would be enough to break the whole MCMC.
##     ## Maybe check if elements are zero prior to calculating the loglik.
##     ## But this might make the process slower. Maybe not slower than using always the
##     ##     full inverse.
##     ms <- 0 ## Here we do not allow for measurement error. Need to implement in the future.
##     ## k <- r
##     ## This is the relevant part. The calculation using C.
##     error <- NULL
##     cholres <- .Call(mvMORPH:::Chol_RPF, V, D, X, as.integer(r), as.integer(ntot), mserr=error, ismserr=as.integer(ms))
##     det <- cholres[[2]]
##     ## Adjusting the format of the data.
##     Xvec <- as.vector(as.matrix(X))
##     residus <- D %*% root - Xvec
##     quad <- .Call(mvMORPH:::Chol_RPF_quadprod, cholres[[1]], residus, as.integer(ntot))
##     logl <- -.5*quad-.5*as.numeric(det)-.5*(ntot*log(2*pi))
##     return(logl)
## }
