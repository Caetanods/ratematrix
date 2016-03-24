##' Make convergence test for the MCMC chain.
##'
##' This function makes convergence test for the MCMC chains. It uses the Heidel test from the package 'coda'.
##'    A better test of convergence is the Gelman's R, but this test requires more than one MCMC chain,
##'    recommended to start for distinct start values of the chain for better accessing the convergence.
##'    This other test will be implemented soon.
##' To calculate the convergence this function will check if each value of each cell of the matrices have
##'    converged. It is hard to find a better way to check for convergence. This is the same approach used by
##'    MCMCglmm and the function to test the significance of the difference between two posterior distributions
##'    of matrices. This approach is most likely more conservative than would be looking to the matrix as
##'    a whole.
##' @title Convergence test.
##' @param mcmc.chain The MCMC chain. Output of the function to read the chain. 
##' @param out The MCMC result output. This is a list with several details about the MCMC run.
##' @return A matrix with colums varying dependent of the number of matrices fitted to the tree. The matrix
##'    has two lines, each is a component of the Heidel diagnostic. TRUE is when the convergence test passed and
##'    FALSE is when it falied.
##' @export
checkConvergence <- function(mcmc.chain, out){
    ## mcmc.chain = A MCMC chain read by the read function.
    ## out = The correspondent chain description. Output of the MCMC function.
    if(out$k == 1){
        root.mcmc <- coda::mcmc( mcmc.chain[[1]] )
        matrix.mcmc <- coda::mcmc( do.call(rbind, lapply(mcmc.chain[[2]], c) ) )
        ## Using the Heidelberger and Welch diagnostic for convergence.
        ## Note that better would be to use Gelman's R. But we only have one chain.
        hei.root <- ratematrix:::passed.heidel(root.mcmc)
        hei.matrix <- ratematrix:::passed.heidel(matrix.mcmc)
        diag <- data.frame(hei.root, hei.matrix)
        colnames(diag) <- c("root","matrix")
        return(diag)
    }
    if(out$k > 1){
        root.mcmc <- coda::mcmc( mcmc.chain[[1]] )
        hei.root <- ratematrix:::passed.heidel(root.mcmc)
        hei.root <- data.frame(hei.root)
        diag <- list()
        for(i in 1:out$k){            
            matrix.mcmc <- coda::mcmc( do.call(rbind, lapply(mcmc.chain[[2]][[i]], c) ) )
            hei.matrix <- ratematrix:::passed.heidel(matrix.mcmc)
            diag[[i]] <- data.frame(hei.matrix)
        }
        res.mat <- do.call(cbind, diag)
        diag.res <- cbind(hei.root, res.mat)
        colnames(diag.res) <- c("root", paste("matrix_", 1:out$k, sep="") )
        return(diag.res)
    }
}
