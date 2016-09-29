##' Test difference between fitted evolutionary rate matrices.
##'
##' This functions performs a test to check whether the posterior distribution of the fitted matrices are
##'    different. When the test is significant this means that the posterior distribution of the elements
##'    of the evolutionary rate matrices does not overlap more than 5%. This test is not a p value, since
##'    this is not an estimate of a probability of deviance from a given null distribution. On the other hand,
##'    the test relies on the overlap of the posterior distribution. It assumes that when the posterior
##'    distribution of two or more paramaters do not overlap, then there is a significative difference
##'    between the parameters. Note that the test will applied to all pairwise combinations of the
##'    matrices fitted to the tree.
##' @title Test different of posterior distributions of R matrices.
##' @param mcmc.chain The MCMC chain produced by the read function.
##' @param out The output of the MCMC chain. This is a list with a series of important informations about the MCMC chain.
##' @param plot logical. Whether plots should be saved as .pdf files in the current directory. The number of the plot increases in function of the number of matrices fitted to the tree. All plots will have a pattern like "post_diff_mat_test_X.pdf", where X is a single digit number.
##' @param file string. The name of the file to save the plot. Default is 'NULL' and will plot in the R standard plotting device. If a string is passed to it, then it will plot to a '.pdf' file. This is going to be used as a the start of the string. The number of plots is equal to the combination 2 by 2 of the number of the matrices fitted to the phylogentic tree.
##' @return Return a matrix with the value of the test statistics. Values bellow 0.05 supports a difference between the posterior distribution of R matrices.
##' @export
testRateMatrix <- function(mcmc.chain, out, plot = TRUE, file = NULL){
    comb <- combn(1:out$p, 2)
    mean.diff <- list()
    for(i in ncol(comb)){
        mat1 <- t( sapply(mcmc.chain$matrix[[comb[1,i]]], function(x) c(x) ) )
        mat2 <- t( sapply(mcmc.chain$matrix[[comb[2,i]]], function(x) c(x) ) )
        mat.diff <- mat1 - mat2
        mean.diff[[i]] <- apply(mat.diff, 1, mean)
    }
    cdf.list <- lapply(mean.diff, FUN = ecdf)
    qq.list <- lapply(cdf.list, FUN = function(x) x(0) )
    test <- lapply(qq.list, FUN = function(x) 2*apply(cbind(x, 1-x), 1, min) )
    if(plot == TRUE){
        if( is.null(file) ){
            hist(mean.diff[[i]], main = paste("mat #",comb[1,i]," vs. mat #",comb[2,i]," test value ="
                                              , round(test[[i]], digits = 3))
                   , freq = FALSE, xlab = "Posterior mean difference between matrices"
               , col = "grey", border = "white")
            abline(v = 0, col = "red", lty = 3, lwd = 3)
            lines( density(mean.diff[[i]]) )
        } else {
            for(i in ncol(comb)){
                pdf( paste(file, i, ".pdf", sep="") )
                hist(mean.diff[[i]], main = paste("mat #",comb[1,i]," vs. mat #",comb[2,i]," test value ="
                                                , round(test[[i]], digits = 3))
                   , freq = FALSE, xlab = "Posterior mean difference between matrices"
                   , col = "grey", border = "white")
                abline(v = 0, col = "red", lty = 3, lwd = 3)
                lines( density(mean.diff[[i]]) )
                dev.off()
            }
        }
    }
    test.dt <- do.call(cbind, test)
    colnames(test.dt) <- paste("mat #",comb[1,], " x #", comb[2,], sep="")
    rownames(test.dt) <- "test value"
    return(test.dt)
}
