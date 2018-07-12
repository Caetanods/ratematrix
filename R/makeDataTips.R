##' Transform data vector into a data matrix for the discrete model
##'
##' Transform data vector into a data matrix for the discrete model
##' @title Transform data vector into a data matrix for the discrete model
##' @param X 
##' @return A data matrix.
##' @author daniel
##' @noRd
makeDataTips <- function(X, states){
    ## Function to correct the format of the data to pass to the likelihood.
    X.mat <- sapply(states, function(y) as.numeric(X == y) )
    rownames(X.mat) <- names(X)
    colnames(X.mat) <- states
    return( X.mat )
}
