##' Log likelihood function for two or more evolutionary rate matrices.
##'
##' Internal function to be used by the MCMC. This function uses the 'rpf' method to calculate the approximation of the inverse of matrices as impplemented in mvMORPH.
##' @title Log likelihood for multiple matrices.
##' @param X matrix. The matrix of trait data.
##' @param root numeric. The phylogenetic mean.
##' @param R.m list. List of matrices fitted to the tree with length equal to 'cache.data$p'.
##' @param C.m list. List of C matrices with length equal to 'cache.data$p'. This is calculated by the function 'multiC' of the 'phytools' package.
##' @param D matrix. The design matrix. This is a fixed term, dependent on the phylogeny.
##' @param n numeric. Number of tips in the phylogeny.
##' @param r numeric. Number of traits, the dimension of the R matrix.
##' @param p numeric. Number of R matrices fitted to the tree.
##' @param method string. Method to be used to calculate the likelihood. Right now only the option 'rpf' is working. More options will be implemented in the future.
##' @return numeric. The log likelhood.
##' @references Revell, L. J., and D. C. Collar. 2009. Phylogenetic Analysis of the Evolutionary Correlation Using Likelihood. Evolution 63:1090–1100.
multiR.loglik <- function(X, root, R.m, C.m, D, n, r, p){
    ## Default method is "rpf"
    method <- method[1]
    V.m <- lapply(1:p, function(x) kronecker(R.m[[x]], C.m[[x]]) )
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
