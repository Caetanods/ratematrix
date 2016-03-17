##' Calculate the effective variance of a list of matrices.
##'
##' This funtion is not optimized to be called by the user.
##' @title Effective variance.
##' @param mat list. list of covariance matrices.
##' @param burn numeric. proportion of the burnin.
##' @param thin numeric. the thinning to be applied in the posterior distribution. Same of the argument 'by' of the function 'seq'.
##' @return numeric. Return a vector with the effective variance index of the matrices.
##' @author Daniel S. Caetano, \email{caetanods1@@gmail.com}
effective.var <- function(mat, burn=0.5, thin=1){
    ## Calculates the effective variance of a list of matrices.
    ## mat = list of covariance matrices.
    ## burn = proportion of burnin.
    ## thin = thinning, same of 'by' in function 'seq'.
    ll <- length(mat)
    if(burn == 0){
        out <- round( ll / thin )
        ss <- round(seq(from = 1, to = ll, length.out = out))
    } else{
        bb <- round( ll * burn )
        out <- round( ( ll - bb ) / thin )
        ss <- round(seq(from = bb, to = ll, length.out = out))
    }
    k <- nrow(mat[[1]])
    detMat <- unlist(lapply(mat[ss], det))
    return(detMat^(1/k))
}
