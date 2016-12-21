##' Check the convergence of the variables using a single chain.
##'
##' Function applies the 'heidel.diag' test from the package 'coda'. For a better test of convergence need to use more than one chain and apply Gelman's R diagnostic.
##' @title Heigel diagnostic of convergence.
##' @param mcmc An 'mcmc' class object as defined in the 'coda' package.
##' @return A vector of length 2 with TRUE or FALSE dependent if the test passed.
##' @importFrom coda heidel.diag
##' @noRd
checkHeidelTest <- function(mcmc){
    ## Check if the variables in the mcmc passed the convergence test.
	## Convergence test for only one chain. For more than one chain use
	##      Gelman's R test.
    ## See the function herdel.diag for more information.
    ## Return a vector with TRUE or FALSE with length 2.
    ## stest = Stationarity test.
    ## htest = Half-width test.
    ## mcmc = An 'mcmc' class object from coda.
    diag <- heidel.diag(mcmc)
    nvar <- dim(mcmc)[2]
    stest <- sum( as.numeric( diag[,1] ) )
    htest <- sum( as.numeric( diag[,4] ) )
    res <- c(stest == nvar, htest == nvar)
    names(res) <- c("stest","htest")
    return(res)
}
