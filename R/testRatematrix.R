##' Function uses summary statistics to test for differences betwwen the posterior distribution of parameter estimates for the evolutionary rate matrix regimes.
##'
##' This functions performs a test to check whether the posterior distribution of the fitted matrices are different. When the test is significant this means that the posterior distribution of the elements of the evolutionary rate matrices does not overlap more than 5\%. This test is not a p value, since this is not an estimate of a probability of deviance from a given null distribution. The test relies on the median overlap of the posterior distribution. It assumes that when the posterior distribution of two or more paramaters do not overlap, then there is a significative difference between the parameters. \cr
##' \cr
##' When a posterior distribution with more than two rate regimes is fitted to the data, the function performs all pairwise combinations tests.
##' @title Test for difference between evolutionary rate matrix estimates
##' @param chain the posterior distribution of parameter estimates as output of the function 'readMCMC'.
##' @param par sets which of the attributes of the rate matrices that are checked by the test. One of 'all', 'correlation', or 'rates' (default is 'all'). Chose 'all' to test the median overlap among all cells of the covariance matrix. Chose 'rates' to check for the median overlap among the evolutionary rates of each trait (the variance vector). Chose 'correlation' to test for the median overlap among the correlations (the distribution of correlation matrices).
##' @return Return a matrix with the value of the test statistics. Values bellow 0.05 supports a difference between the posterior distribution of R matrices.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' \donttest{
##' ## Load data
##' data(centrarchidae)
##' ## Run MCMC. This is just a very short chain.
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=1000)
##' ## Load posterior distribution, make plots and check the log.
##' posterior <- readMCMC(handle, burn=0.25, thin=1)
##' testRatematrix(chain=posterior, par="all")
##' testRatematrix(chain=posterior, par="rates")
##' testRatematrix(chain=posterior, par="correlation")
##' }
testRatematrix <- function(chain, par=c("all","correlation","rates")){

    par <- match.arg(par)

    if( inherits(chain, what=c("ratematrix_single_chain", "ratematrix_multi_chain")) ){
        if( inherits(chain, what=c("ratematrix_single_chain")) ){
            p <- 1
        } else{
            p <- length( chain$matrix )
        }
    } else{
        stop("Arguments need to be the posterior distribution of class 'ratematrix_single_chain' or 'ratematrix_multi_chain'. See 'readMCMC'. \n")
    }
    
    comb <- combn(1:p, 2)
    median.diff <- list()
    if(par == "all"){
        for(i in 1:ncol(comb)){
            mat1 <- t( sapply(chain$matrix[[comb[1,i]]], function(x) c(x) ) )
            mat2 <- t( sapply(chain$matrix[[comb[2,i]]], function(x) c(x) ) )
            mat.diff <- mat1 - mat2
            median.diff[[i]] <- apply(mat.diff, 1, median)
        }
    }
    if(par == "rates"){
        for(i in 1:ncol(comb)){
            mat1 <- t( sapply(chain$matrix[[comb[1,i]]], function(x) c( diag(x) ) ) )
            mat2 <- t( sapply(chain$matrix[[comb[2,i]]], function(x) c( diag(x) ) ) )
            mat.diff <- mat1 - mat2
            median.diff[[i]] <- apply(mat.diff, 1, median)
        }
    }
    if(par == "correlation"){
        ## This part will break if the rate matrix is 2x2. The median do not need to be calculated.
        if( ncol( chain$matrix[[1]][[1]] ) > 2 ){
            upper <- upper.tri( chain$matrix[[1]][[1]] )
            for(i in 1:ncol(comb)){
                mat1 <- t( sapply(chain$matrix[[comb[1,i]]], function(x) c( cov2cor(x)[upper] ) ) )
                mat2 <- t( sapply(chain$matrix[[comb[2,i]]], function(x) c( cov2cor(x)[upper] ) ) )
                mat.diff <- mat1 - mat2
                median.diff[[i]] <- apply(mat.diff, 1, median)
            }
        } else{
            upper <- upper.tri( chain$matrix[[1]][[1]] )
            for(i in 1:ncol(comb)){
                mat1 <- t( sapply(chain$matrix[[comb[1,i]]], function(x) c( cov2cor(x)[upper] ) ) )
                mat2 <- t( sapply(chain$matrix[[comb[2,i]]], function(x) c( cov2cor(x)[upper] ) ) )
                median.diff[[i]] <- mat1 - mat2 ## Not a real 'median', because there is only one correlation.
            }
        }
    }
    cdf.list <- lapply(median.diff, FUN = ecdf)
    qq.list <- lapply(cdf.list, FUN = function(x) x(0) )
    test <- lapply(qq.list, FUN = function(x) 2*apply(cbind(x, 1-x), 1, min) )
    test.dt <- do.call(cbind, test)
    colnames(test.dt) <- paste("mat #",comb[1,], " x #", comb[2,], sep="")
    rownames(test.dt) <- "test value"
    return(test.dt)
}
