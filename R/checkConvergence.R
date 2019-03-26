##' Make convergence test for the MCMC chain. We STRONGLY recommend doing at least two independent searches to test convergence. Please see 'Details' for more information.
##'
##' Function performs convergence tests using the potential scale reduction factor (Gelman's R) or the Heidelberger test. The Gelman's R test will be performed if two or more MCMC chains are provided as input. If only one MCMC chain is provided, then the function will perform the Heidelberger test (and print a message about it). \cr
##' \cr
##' Multiple chains need to be replicates of the same analysis (e.g., multiple runs of the 'ratematrixMCMC' function with the same set of arguments and, in the best scenario, with varying starting points). We recommend users to perform the Gelman's R test by providing two or more independent MCMC chains with different starting points. This test is more robust than the Heidelberger test. The advantage of the Heidelberger test is that it can be used with a single MCMC chain, so it can be useful for a preliminary test prior to running a full convergence analysis with multiple chains. (Our experience shows that performing the Heidelberger test alone can return false convergence results.) Convergence can also be investigated using the 'logAnalizer' and 'computeESS'.\cr
##' \cr
##' The 'Gelman's R' test is based on the potential scale reduction factor which is expected to be equal to 1 when convergence is achieved. If you see values close to 1 (e.g., ~1.01 to 1.05) it means that you just need to get more samples from the MCMC (see 'continueMCMC' function). See more information about each of these tests in the references below and in the documentation for the functions 'coda::gelman.diag' and 'coda::heidel.diag', both from the package 'coda'.
##' @title Performs convergence tests
##' @param ... posterior(s) distribution(s) of parameter estimates. This can be a single MCMC chain or multiple independent chains from the same model. The type of convergence analysis will be dependent on the number of MCMC chains provided as input. See 'Details'.
##' @return The format of the output depends of the type of test performed. The Gelman's R test will return a list with two elements. The first element is a list with the results for the potential scale reduction factor for the root values and the evolutionary rate matrices. The test for the R matrices is performed element by element, the names of the columns show the number of the row and column for each element. The length of this list will depend on the number of rate regimes fitted to the phylogenetic tree. The Heidelberger test also returns a list with two elements, the first element is a table with one column for the root values and each evolutionary rate matrix regime fitted to the tree. The colnames show the type of diagnostic used, the values are whether the test passed or not. The second element of the list, independent of the type of convergence test, is a estimate of the Effective Sample Size for each parameter of the model.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @importFrom coda mcmc effectiveSize mcmc.list gelman.diag
##' @examples
##' \donttest{
##' data(centrarchidae)
##' handle1 <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=10000
##'                           , dir=tempdir())
##' posterior1 <- readMCMC(handle1, burn=0.25, thin=10)
##' handle2 <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=10000
##'                           , dir=tempdir())
##' posterior2 <- readMCMC(handle2, burn=0.25, thin=1)
##' ## Note that these are short chains used here as example only.
##' ## A convergence test using 'Gelman's R' calculated from two independent MCMC chains.
##' checkConvergence(posterior1, posterior2)
##' }
##' @references
##' \describe{
##'   \item{}{Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.}
##'   \item{}{Heidelberger, P and Welch, PD (1983) Simulation run length control in the presence of an initial transient. Opns. Res., 31, 1109-44.}
##' }
checkConvergence <- function(...){

    chains <- list(...)

    if( inherits(chains[[1]], what = c("ratematrix_multi_chain", "ratematrix_single_chain", "ratematrix_poly_chain")) ){
        ## Probably a single posterior as input.
        if( length( chains[[1]] ) == 1 ) stop( "Need two or more posterior distributions to merge." )
    } else if( inherits(chains[[1]], what = "list") ){
        ## Input likely a list of posteriors.
        if( inherits(chains[[1]][[1]], what = c("ratematrix_multi_chain","ratematrix_single_chain","ratematrix_poly_chain") ) ){
            ## Yes, this is a list of posteriors:
            if( length( chains[[1]] ) == 1 ) stop( "Need two or more posterior distributions to merge." )
            class.post <- sapply(chains[[1]], function(x) inherits(x, what = c("ratematrix_multi_chain","ratematrix_single_chain","ratematrix_poly_chain") ))
            if( !all( class.post ) ) stop( "All input objects need to be posterior distributions." )
            chains <- chains[[1]] ## Decrease one level in the listing structure.
        } else{
            stop( "All input objects need to be posterior distributions." )
        }
    } else{
        ## Don't know what type this is. Returns error.
        stop( "All input objects need to be posterior distributions." )
    }
    
    ## Test if all the posteriors provided are of the same class.
    ## Then calculate the value of p, if necessary and proceed.
    ## The type of convergence test will be given by the class of the posterior.
    
    if( sum( sapply(chains, function(x) inherits(x, what=c("ratematrix_single_chain", "ratematrix_multi_chain","ratematrix_poly_chain")) ) ) == length(chains) ){
        
        if( length( unique( sapply(chains, function(x) class(x) ) ) ) > 1 ) stop("Posterior chains need to belong to the same model. \n")

        ## Get the number of rate regimes:
        if( inherits(chains[[1]], what=c("ratematrix_single_chain")) ){
            p <- 1
        }
        if( inherits(chains[[1]], what=c("ratematrix_multi_chain")) ){
            p <- length( chains[[1]]$matrix )
        }
        if( inherits(chains[[1]], what=c("ratematrix_poly_chain")) ){
            if( is.matrix(chains[[1]]$matrix[[1]]) ){
                p <- 1
            } else{
                p <- length( chains[[1]]$matrix )
            }                
        }

        ## Define if the root is included in the model:
        include.root <- ifelse( is.null(chains[[1]]$root), FALSE, TRUE)
        
        if( length(chains) == 1 ){
            method <- "heidel"
        }
        
        if( length(chains) > 1 ){
            method <- "gelman"

            ## Define function to reduce the length of the chains.
            ## This is necessary because the coda functions require chains with the same length.
            set.size <- function(mcmc, n, p, include.root = TRUE){
                res <- mcmc
                if( include.root ){
                    res$root <- mcmc$root[1:n,]
                }
                ## Need to check the number of regimes in order to manipulate the R matrices on the posterior distribution.
                if( p > 1 ){
                    for(i in 1:p){
                        res$matrix[[i]] <- mcmc$matrix[[i]][1:n]
                    }
                } else{
                    ## A single regime analysis.
                        res$matrix <- mcmc$matrix[1:n]
                }
                return( res )
            }

            ## Get the length of each chain:
            if( p == 1 ){
                ll <- sapply(chains, function(x) length(x$matrix) )
            } else{
                ll <- sapply(chains, function(x) length(x$matrix[[1]]) )
            }
            
            if( as.logical(max(ll) - min(ll)) ){
                ## Need to set the length of all the chains to be the same.
                warning( "MCMC chains have different lengths. Pruning all chains to the minimum length observed." )
                chains <- lapply(chains, function(x) set.size(mcmc=x, n=min(ll), p=p, include.root=include.root) )
            }
        }
        
    } else{
        stop("Arguments need to be output of 'readMCMC' function. Of class 'ratematrix_single_chain', 'ratematrix_multi_chain', or 'ratematrix_poly_chain'. \n")
    }

    if(p == 1){
        k <- ncol( chains[[1]]$matrix[[1]] )
        if(method == "heidel"){
            ## There is only a single chain.
            matrix.mcmc <- coda::mcmc( do.call(rbind, lapply(chains[[1]]$matrix, c) ) )
            ## Using the Heidelberger and Welch diagnostic for convergence.
            ## Note that better would be to use Gelman's R. But we only have one chain.
            hei.matrix <- checkHeidelTest(matrix.mcmc)
            ess.matrix <- effectiveSize(matrix.mcmc)
            
            if( include.root ){
                ## This returns only the ratematrix results and ignores the tip and anc samples (even if they are present).
                root.mcmc <- coda::mcmc( chains[[1]]$root )
                hei.root <- checkHeidelTest(root.mcmc)
                ess.root <- effectiveSize(root.mcmc)
                ess <- c(ess.root, ess.matrix)
                names(ess) <- c( paste("root_",1:length(ess.root),sep=""), paste("matrix_cel_",1:length(ess.matrix),sep=""))
                diag <- data.frame(hei.root, hei.matrix)
                colnames(diag) <- c("root","matrix")
            }

            if( inherits(chains[[1]], what=c("ratematrix_poly_chain")) ){
                ## Need to get the test for the tip and node samples.
                tip.samples.mcmc <- coda::mcmc( chains[[1]]$tip_samples )
                hei.tip.samples <- checkHeidelTest(tip.samples.mcmc)
                ess.tip.samples <- effectiveSize(tip.samples.mcmc)
                anc.samples.mcmc <- coda::mcmc( chains[[1]]$anc_samples )
                hei.anc.samples <- checkHeidelTest(anc.samples.mcmc)
                ess.anc.samples <- effectiveSize(anc.samples.mcmc)
                ess <- c(ess.matrix, ess.tip.samples, ess.anc.samples)
                names(ess) <- c( paste("matrix_cel_",1:length(ess.matrix),sep="")
                              , colnames( chains[[1]]$tip_samples )
                              , colnames( chains[[1]]$anc_samples ) )
                diag <- data.frame(hei.matrix, hei.tip.samples, hei.anc.samples)
                colnames(diag) <- c("matrix", "tip_samples", "anc_samples")
            }

            return( list(heidel=diag, ess=ess) )
            
        }
        if(method == "gelman"){
            ## Single rate regime but multiple chains.
            matrix.mcmc <- lapply(chains, function(x) coda::mcmc( do.call(rbind, lapply(x$matrix, c) ) ) )
            matrix.mcmc.list <- mcmc.list(matrix.mcmc)
            diag.matrix <- gelman.diag(matrix.mcmc.list, autoburnin=FALSE, multivariate=FALSE)
            nm <- expand.grid(1:k, 1:k)
            rownames(diag.matrix$psrf) <- sprintf('%s,%s', nm[,2], nm[,1])
            ess.matrix <- effectiveSize(matrix.mcmc.list)

            if( include.root ){
                ## Include the Gelman diag for the root and ignore the tip and anc samples.
                root.mcmc <- lapply(chains, function(x) coda::mcmc( x$root ) )
                root.mcmc.list <- mcmc.list(root.mcmc)
                ess.root <- effectiveSize(root.mcmc.list)
                diag.root <- gelman.diag(root.mcmc.list, autoburnin=FALSE, multivariate=FALSE)
                ess <- c(ess.root, ess.matrix)
                names(ess) <- c( colnames(chains[[1]]$root), paste("matrix_cel_", sprintf('%s,%s', nm[,2], nm[,1]), sep=""))
                gelman.diag <- list(root=diag.root, ratematrix=diag.matrix)
            }

            if( inherits(chains[[1]], what=c("ratematrix_poly_chain")) ){
                ## Need to get the test for the tip and node samples.
                tip.samples.mcmc <- lapply(chains, function(x) coda::mcmc( x$tip_samples ) )
                tip.samples.mcmc.list <- mcmc.list( tip.samples.mcmc )
                diag.tip.samples <- gelman.diag(tip.samples.mcmc.list, autoburnin=FALSE, multivariate=FALSE)
                ess.tip.samples <- effectiveSize(tip.samples.mcmc.list)

                anc.samples.mcmc <- lapply(chains, function(x) coda::mcmc( x$anc_samples ) )
                anc.samples.mcmc.list <- mcmc.list( anc.samples.mcmc )
                diag.anc.samples <- gelman.diag(anc.samples.mcmc.list, autoburnin=FALSE, multivariate=FALSE)
                ess.anc.samples <- effectiveSize(anc.samples.mcmc.list)
                ess <- c(ess.matrix, ess.tip.samples, ess.anc.samples)
                names(ess) <- c( paste("matrix_cel_",1:length(ess.matrix),sep="")
                              , colnames( chains[[1]]$tip_samples )
                              , colnames( chains[[1]]$anc_samples ) )
                gelman.diag <- list(ratematrix=diag.matrix, tip_samples=diag.tip.samples, anc_samples=diag.anc.samples)
            }
            
            return( list(gelman=gelman.diag, ess=ess) )
        }
    }
    
    if(p > 1){
        k <- ncol( chains[[1]]$matrix[[1]][[1]] )
        if(method == "heidel"){
            diag <- list()
            ess.matrix <- list()
            for(i in 1:p){
                matrix.mcmc <- coda::mcmc( do.call(rbind, lapply(chains[[1]]$matrix[[i]], c) ) )
                hei.matrix <- checkHeidelTest(matrix.mcmc)
                ess.matrix[[i]] <- effectiveSize(matrix.mcmc)
                diag[[i]] <- data.frame(hei.matrix)
            }
            res.mat <- do.call(cbind, diag)
            
            if( include.root ){
                ## Include the root and ignore samples for the tip and ancestral states.
                root.mcmc <- coda::mcmc( chains[[1]][[1]] )
                hei.root <- checkHeidelTest(root.mcmc)
                ess.root <- effectiveSize(root.mcmc)
                hei.root <- data.frame(hei.root)
                diag.res <- cbind(hei.root, res.mat)
                colnames(diag.res) <- c("root", names( chains[[1]]$matrix ) )
                nn <- rep(NA, times=length(ess.matrix[[1]]) )
                ess.matrix.mat <- do.call(cbind, ess.matrix)
                nn[1:length(ess.root)] <- ess.root
                ess <- cbind(nn, ess.matrix.mat)
                colnames(ess) <- c("root", names( chains[[1]]$matrix ) )
                return(list(heidel=diag.res,ess=ess))
            }

            if( inherits(chains[[1]], what=c("ratematrix_poly_chain")) ){
                ## Need to get the test for the tip and node samples.
                tip.samples.mcmc <- coda::mcmc( chains[[1]]$tip_samples )
                hei.tip.samples <- data.frame( checkHeidelTest(tip.samples.mcmc) )
                ess.tip.samples <- effectiveSize(tip.samples.mcmc)
                anc.samples.mcmc <- coda::mcmc( chains[[1]]$anc_samples )
                hei.anc.samples <- data.frame( checkHeidelTest(anc.samples.mcmc) )
                ess.anc.samples <- effectiveSize(anc.samples.mcmc)
                diag.res <- cbind(res.mat, hei.tip.samples, hei.anc.samples)
                colnames(diag.res) <- c(paste0( "regime_", names( chains[[1]]$matrix ) ), "tip_samples", "anc_samples")
                ess.matrix.mat <- do.call(cbind, ess.matrix)
                colnames(ess.matrix.mat) <- names( chains[[1]]$matrix )
                return( list(heidel=diag.res, ess_rates=ess.matrix.mat) )
            }
        }
        
        if(method == "gelman"){
            diag.matrix <- list()
            ess.matrix <- list()
            for(i in 1:p){
                matrix.mcmc <- lapply(chains, function(x) coda::mcmc( do.call(rbind, lapply(x$matrix[[i]], c))))
                matrix.mcmc.list <- mcmc.list(matrix.mcmc)
                ess.matrix[[i]] <- effectiveSize(matrix.mcmc.list)
                diag.matrix[[i]] <- gelman.diag(matrix.mcmc.list, autoburnin=FALSE, multivariate=FALSE)
                nm <- expand.grid(1:k, 1:k)
                nm.comb <- sprintf('%s,%s', nm[,2], nm[,1])
                rownames(diag.matrix[[i]]$psrf) <- nm.comb
            }
            ## Keep the names to use again below.
            nm.comb <- sprintf('%s,%s', nm[,2], nm[,1])
            names(diag.matrix) <- names( chains[[1]]$matrix )
            ess.matrix.mat <- do.call(cbind, ess.matrix)

            if( include.root ){
                root.mcmc <- lapply(chains, function(x) coda::mcmc( x[[1]] ) )
                root.mcmc.list <- coda::mcmc.list(root.mcmc)
                diag.root <- gelman.diag(root.mcmc.list, autoburnin=FALSE, multivariate=FALSE)
                ess.root <- effectiveSize(root.mcmc.list)
                res.root <- list(diag_root=diag.root)
                res <- c(res.root, diag.matrix)
                nn <- rep(NA, times=length(ess.matrix[[1]]) )
                nn[1:length(ess.root)] <- ess.root
                ess <- cbind(nn, ess.matrix.mat)
                colnames(ess) <- c("root", names( chains[[1]]$matrix ) )            
                return( list( gelman=res, ess=ess ) )
            }

            if( inherits(chains[[1]], what=c("ratematrix_poly_chain")) ){
                ## Need to get the test for the tip and node samples.
                tip.samples.mcmc <- lapply(chains, function(x) coda::mcmc( x$tip_samples ) )
                tip.samples.mcmc.list <- mcmc.list( tip.samples.mcmc )
                diag.tip.samples <- gelman.diag(tip.samples.mcmc.list, autoburnin=FALSE, multivariate=FALSE)
                ess.tip.samples <- effectiveSize(tip.samples.mcmc.list)

                anc.samples.mcmc <- lapply(chains, function(x) coda::mcmc( x$anc_samples ) )
                anc.samples.mcmc.list <- mcmc.list( anc.samples.mcmc )
                diag.anc.samples <- gelman.diag(anc.samples.mcmc.list, autoburnin=FALSE, multivariate=FALSE)
                ess.anc.samples <- effectiveSize(anc.samples.mcmc.list)
                ess <- ess.matrix.mat                
                rownames(ess) <- paste("rate_", nm.comb, sep="")
                colnames(ess) <- paste0("regime_", names( chains[[1]]$matrix ))
                
                gelman.diag <- list(ratematrix=diag.matrix, tip_samples=diag.tip.samples, anc_samples=diag.anc.samples)
                return( list( gelman = gelman.diag, ess = ess ) )
            }
            
        }
    }
}
