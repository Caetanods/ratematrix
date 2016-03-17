##' Log likelihood function for one evolutionary rate matrices.
##'
##' Internal function to be used by the MCMC.
##' @title Log likelihood for one evolutionary rate matrix fitted to the data and tree.
##' @param X matrix. The matrix of trait data.
##' @param root numeric. The phylogenetic mean.
##' @param R matrix. The covariance matrix fitted to the tree and data.
##' @param C matrix. The phylogenetic covariance matrix. This is a fixed term.
##' @param D matrix. The design matrix. This is a fixed term, dependent on the phylogeny.
##' @param n numeric. Number of tips in the phylogeny.
##' @param r numeric. Number of traits, the dimension of the R matrix.
##' @param method string. Method to be used to calculate the likelihood. Right now only the option 'rpf' is working. More options will be implemented in the future.
##' @return numeric. The log likelhood.
##' @references Revell, L. J., and D. C. Collar. 2009. Phylogenetic Analysis of the Evolutionary Correlation Using Likelihood. Evolution 63:1090–1100.
singleR.loglik <- function(X, root, R, C, D, n, r, method=c("rpf","sparse")){
    ## The loglik function for a single R matrix in a tree.
    ## This is based on the faster options for calculating the likelihood implemented in mvMORPH.
    ## X = The matrix of the trait data.
    ## root = The phylogenetic mean.
    ## R = The R matrix.
    ## C = The phylogenetic covariance matrix. This is a fixed term.
    ## D = The desing matrix. This is a fixed term.
    ## y = Column vector of trait values.
    ## b = Design matrix * the phylogenetic mean vector. This is a nuisance parameter.
    ##     I can both set the phylogenetic mean to be a fixed value or try to make proposals
    ##     to it.
    ## n = Number of tips of the phylogeny.
    ## r = Number of traits. Dimension of the R matrix.

    ## Default method is "rpf"
    method <- method[1]
    V <- kronecker(R, C) ## In mvMORPH this is the parameter 'tree'.
    ntot <- n*r
    if(method=="rpf"){
        ## This method will not accept matrices with elements equal to zero.
        ## Unclear if this would be enough to break the whole MCMC.
        ## Maybe check if elements are zero prior to calculating the loglik.
        ## But this might make the processo slower. Maybe not slower than using always the
        ##     full inverse.
        ms <- 0 ## Here we do not allow for measurement error. Need to implement in the future.
        ## This is the relevant part. The calculation using C.
        error <- NULL
        cholres <- .Call("Chol_RPF", V, D, X, as.integer(r), as.integer(ntot), mserr=error, ismserr=as.integer(ms))
        det <- cholres[[2]]
        ## Adjusting the format of the data.
        Xvec <- as.vector(as.matrix(X))        
        residus <- D %*% root - Xvec
        quad <- .Call("Chol_RPF_quadprod", cholres[[1]], residus, as.integer(ntot))
        logl <- -.5*quad-.5*as.numeric(det)-.5*(ntot*log(2*pi))
    }
    if(method=="sparse"){
        ## Sparce method will take more time than the normal inverse if the matrix is
        ##        not sparse enough. Problem is that the matrix change its format in the MCMC.
        ## This part is still not funtional. Need to test with more care.
        spambig <- as.spam(V)
        U <- chol(spambig)
        beta <- root
        Xvec <- as.vector(as.matrix(X))
        res <- D %*% root - Xvec
        vec1 <- forwardsolve(U,res)
        a <- sum(vec1^2)
        DET <- determinant(U)
        logl <- -.5*(a)-.5*as.numeric(DET$modulus*2)-.5*(ntot*log(2*pi))
    }
    return(logl)
}
