##' Make convergence test for the MCMC chain.
##'
##' This function makes convergence test for the MCMC chains. It uses the Heidel test or the potential scale reduction factor (Gelman's R), both from the package 'coda'.
##' To calculate the convergence this function will check if each value of each cell of the matrices have
##'    converged.
##' @title Convergence test.
##' @param mcmc.chain The MCMC chain. Output of the function to read the chain. It works with any number of matrices fitted to the tree. If 'method="gelman"' the argument need to be a list of the MCMC chain results. Each element of the list will be the result of an independent MCMC run starting from a different starting position.
##' @param p numeric. The number of evolutionary rate matrices fitted to the tree.
##' @param method The method to be used to check convergence. If 'method="heidel"' then the Heidelberger and Welch diagnostic for convergence will be conducted. This test requires a single MCMC chain. If 'method="gelman"' then the potential scale reduction factor will be calculated. This method requires at least two MCMC chains for the same model starting from different starting points. The later is recommended.
##' @param multivariate logical. Whether the 'gelman.diag' function from the package 'coda' should perform the diagnostics using the multivariate version. The default value is TRUE. However, for some matrix configurations the function can return an error that the '(...) leading (...) of order (...) is not positive definite'. Then, setting this argument as FALSE will solve the issue.
##' @return Function returns a list with two elements. If the 'method="heidel"' the first element of the list is a matrix with colums varying dependent of the number of matrices fitted to the tree. This matrix has two lines, each is a component of the Heidel diagnostic. TRUE is when the convergence test passed and FALSE is when it falied. If the 'method="gelman"' the first element of the list is another list with the results of the 'gelman.diag' function. This list will have length equal to 1+p: first element is the result for the root value and the following elements are results for the fitted R matrices. The second element of the list returned by the function is a matrix with the effective sample sizes for the root and each cell of the fitted R matrices.
##' @export
##' @importFrom coda mcmc effectiveSize mcmc.list gelman.diag
checkConvergence <- function(mcmc.chain, p, method=c("heidel","gelman"), multivariate=TRUE){

    if(p == 1){
        if(method == "heidel"){
            root.mcmc <- coda::mcmc( mcmc.chain[[1]] )
            matrix.mcmc <- coda::mcmc( do.call(rbind, lapply(mcmc.chain[[2]], c) ) )
            ## Using the Heidelberger and Welch diagnostic for convergence.
            ## Note that better would be to use Gelman's R. But we only have one chain.
            hei.root <- .checkHeidelTest(root.mcmc)
            ess.root <- coda::effectiveSize(root.mcmc)
            hei.matrix <- .checkHeidelTest(matrix.mcmc)
            ess.matrix <- coda::effectiveSize(matrix.mcmc)
            ess <- c(ess.root, ess.matrix)
            names(ess) <- c( paste("root_",1:length(ess.root),sep=""), paste("matrix_cel_",1:length(ess.matrix),sep=""))
            diag <- data.frame(hei.root, hei.matrix)
            colnames(diag) <- c("root","matrix")
            return( list(heidel=diag,ess=ess) )
        }
        if(method == "gelman"){
            root.mcmc <- lapply(mcmc.chain, function(x) coda::mcmc( x[[1]] ) )
            root.mcmc.list <- coda::mcmc.list(root.mcmc)
            matrix.mcmc <- lapply(mcmc.chain, function(x) coda::mcmc( do.call(rbind, lapply(x[[2]], c) ) ) )
            matrix.mcmc.list <- coda::mcmc.list(matrix.mcmc)
            diag.root <- coda::gelman.diag(root.mcmc.list, autoburnin=FALSE, multivariate=multivariate)
            diag.matrix <- coda::gelman.diag(matrix.mcmc.list, autoburnin=FALSE, multivariate=multivariate)
            ess.root <- coda::effectiveSize(root.mcmc.list)
            ess.matrix <- coda::effectiveSize(matrix.mcmc.list)
            ess <- c(ess.root, ess.matrix)
            names(ess) <- c( paste("root_",1:length(ess.root),sep=""), paste("matrix_cel_",1:length(ess.matrix),sep=""))
            return( list(diag_root=diag.root, diag_matrix=diag.matrix, ess=ess) )
        } else {
            stop("'method' need to be one of 'heidel' or 'gelman'.")
        }
    }
    if(p > 1){
        if(method == "heidel"){
            root.mcmc <- coda::mcmc( mcmc.chain[[1]] )
            hei.root <- .checkHeidelTest(root.mcmc)
            ess.root <- coda::effectiveSize(root.mcmc)
            hei.root <- data.frame(hei.root)
            diag <- list()
            ess.matrix <- list()
            for(i in 1:p){
                matrix.mcmc <- coda::mcmc( do.call(rbind, lapply(mcmc.chain[[2]][[i]], c) ) )
                hei.matrix <- .checkHeidelTest(matrix.mcmc)
                ess.matrix[[i]] <- coda::effectiveSize(matrix.mcmc)
                diag[[i]] <- data.frame(hei.matrix)
            }
            res.mat <- do.call(cbind, diag)
            diag.res <- cbind(hei.root, res.mat)
            colnames(diag.res) <- c("root", paste("matrix_", 1:p, sep="") )
            nn <- rep(NA, times= length(ess.matrix[[1]]) )
            ess.matrix.mat <- do.call(cbind, ess.matrix)
            nn[1:length(ess.root)] <- ess.root
            ess <- cbind(nn, ess.matrix.mat)
            colnames(ess) <- c("root", paste("matrix_", 1:p, sep="") )
            return(list(heidel=diag.res,ess=ess))
        }
        if(method == "gelman"){
            root.mcmc <- lapply(mcmc.chain, function(x) coda::mcmc( x[[1]] ) )
            root.mcmc.list <- coda::mcmc.list(root.mcmc)
            diag.root <- coda::gelman.diag(root.mcmc.list, autoburnin=FALSE, multivariate=multivariate)
            ess.root <- coda::effectiveSize(root.mcmc.list)
            diag.matrix <- list()
            ess.matrix <- list()
            for(i in 1:p){
                matrix.mcmc <- lapply(mcmc.chain, function(x) coda::mcmc( do.call(rbind, lapply(x[[2]][[i]], c))))
                matrix.mcmc.list <- coda::mcmc.list(matrix.mcmc)
                ess.matrix[[i]] <- coda::effectiveSize(matrix.mcmc.list)
                diag.matrix[[i]] <- coda::gelman.diag(matrix.mcmc.list,autoburnin=FALSE,multivariate=multivariate)
            }
            names(diag.matrix) <- paste("diag_matrix_",1:p,sep="")
            res.root <- list(diag_root=diag.root)
            res <- c(res.root, diag.matrix)
            nn <- rep(NA, times= length(ess.matrix[[1]]) )
            ess.matrix.mat <- do.call(cbind, ess.matrix)
            nn[1:length(ess.root)] <- ess.root
            ess <- cbind(nn, ess.matrix.mat)
            colnames(ess) <- c("root", paste("matrix_", 1:p, sep="") )            
            return( list( gelman=res, ess=ess ) )
            } else {
            stop("'method' need to be one of 'heidel' or 'gelman'.")
        }
    }
}
