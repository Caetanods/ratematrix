##' Function to sample a starting point for the MCMC. This function sample the vector of root values, the correlation matrix and the vector of standard deviations in a format that can be used to start the MCMC chain ('ratematrixMCMC' function).
##'
##' The functionality of this function is already integraded within 'ratematrixMCMC'. Sampling a starting point using this function is equivalent to setting start="prior_sample" and prior="unif" in the 'ratematrixMCMC' function.\cr
##' \cr
##' The format for the start of the MCMC is a list. The first element of the list is a vector with the root values for each trait. The second element is a correlation matrix and the third element is a vector of standard deviations. The root values and the vector of standard deviations need to have the same length, the correlation matrix needs to have ncol and nrow equal to the length of the root values.
##' @title Sample a starting point for the MCMC.
##' @param k number of traits in the model.
##' @param p number of evolutionary rate matrices fitted to the tree.
##' @return A list with the starting point for a MCMC chain.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' data(centrarchidae)
##' ## Sample starting point for two traits and two regimes.
##' start <- makeStart(k=2, p=2)
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=2000, start=start)
##' posterior <- readMCMC(handle)
##' plotRatematrix(posterior)
##' @importFrom corpcor decompose.cov
makeStart <- function(k, p){
    mn <- stats::runif(k, min=-100, max=100)
    if( p == 1 ){
        corr <- riwish(v=k+1, S=diag(nrow=k))
        sd <- sqrt( decompose.cov(corr)$v )
    } else{
        corr <- list()
        sd <- list()
        for(i in 1:p){
            corr[[i]] <- riwish(v=k+1, S=diag(nrow=k))
            sd[[i]] <- sqrt( decompose.cov(corr[[i]])$v )
        }
    }
    start <- list(mn, corr, sd)
    return(start)
}
