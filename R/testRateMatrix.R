##' Test difference between fitted evolutionary rate matrices.
##'
##' This functions performs a test to check whether the posterior distribution of the fitted matrices are
##'    different. When the test is significant this means that the posterior distribution of the elements
##'    of the evolutionary rate matrices does not overlap more than 5%. This test is not a p value, since
##'    this is not an estimate of a probability of deviance from a given null distribution. On the other hand,
##'    the test relies on the median overlap of the posterior distribution. It assumes that when the posterior
##'    distribution of two or more paramaters do not overlap, then there is a significative difference
##'    between the parameters. Note that the test will applied to all pairwise combinations of the
##'    matrices fitted to the tree. You can use the 'par' parameter to set which are the elements of the R matrix
##'    you are interested in comparing.
##' @title Test for difference between evolutionary rate matrix estimates
##' @param chain The MCMC chain produced by the read function.
##' @param par This parameter sets the attribute of the R matrices that is checked by the test. Chose 'all' to test the median overlap among all cells of the covariance matrix. Chose 'rates' to check for the median overlap among the evolutionary rates of each trait (variance vector). Chose 'correlation' to test for the median overlap among the correlations (correlation matrix).
##' @return Return a matrix with the value of the test statistics. Values bellow 0.05 supports a difference between the posterior distribution of R matrices.
##' @export
testRatematrix <- function(chain, par="all"){
    ## Right now the function offer for a plot. Better than that would be to print the partial results, the difference between each cell and then the median of the difference. Also return this as the elements of a list.
    ## Right now the function is doing the tests but it is not returning the differences trait by trait. Need to find a way to do that.

    ## Need to find the number of regimes fitted to the tree automatically. Also perform some tests to make sure that the imput data is of the correct format.

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
        for(i in ncol(comb)){
            mat1 <- t( sapply(chain$matrix[[comb[1,i]]], function(x) c(x) ) )
            mat2 <- t( sapply(chain$matrix[[comb[2,i]]], function(x) c(x) ) )
            mat.diff <- mat1 - mat2
            median.diff[[i]] <- apply(mat.diff, 1, median)
        }
    }
    if(par == "rates"){
        for(i in ncol(comb)){
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
            for(i in ncol(comb)){
                mat1 <- t( sapply(chain$matrix[[comb[1,i]]], function(x) c( cov2cor(x)[upper] ) ) )
                mat2 <- t( sapply(chain$matrix[[comb[2,i]]], function(x) c( cov2cor(x)[upper] ) ) )
                mat.diff <- mat1 - mat2
                median.diff[[i]] <- apply(mat.diff, 1, median)
            }
        } else{
            upper <- upper.tri( chain$matrix[[1]][[1]] )
            for(i in ncol(comb)){
                mat1 <- t( sapply(chain$matrix[[comb[1,i]]], function(x) c( cov2cor(x)[upper] ) ) )
                mat2 <- t( sapply(chain$matrix[[comb[2,i]]], function(x) c( cov2cor(x)[upper] ) ) )
                median.diff[[1]] <- mat1 - mat2 ## Not a real 'median', because there is only one correlation.
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
