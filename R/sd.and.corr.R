##' Decomposes a covariance matrix in the correlation and standard deviation components.
##'
##' The function applies 'cov2cor' to calculate the correlation matrix from the covariance matrix and Takes the 'sqrt' of the diagonal elements to extract the standard deviations.
##' @title Matrix to correlation and standard deviation.
##' @param mat list. A list of posterior distribution of evolutionary rate matrices (R) or covariance matrices.
##' @param burn numeric. The proportion of the burnin to take from the chain.
##' @param thin thinning. The thinning to be applied to the posterior distribution of matrices. This is the same as the argument 'by' of the function 'seq'.
##' @return Returns a list with two components. 'corr' is a list of the correlation matrices and 'sd' is a matrix with the standard deviations in the rows in the same order as the correlation matrices in 'corr'.
sd.and.corr <- function(mat, burn=0.5, thin=1){
    ## Calculate the corr matrix and vector of standard deviations.
    ## Return list with two elements. First is list of correlation matrices.
    ##      Second is matrix with standard deviations as rows.
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
    mcor <- lapply(mat[ss], cov2cor)
    sd <- sapply(mat[ss], function(x) sqrt( diag(x) ) )
    return( list(corr = mcor, sd = t(sd) ) )
}
