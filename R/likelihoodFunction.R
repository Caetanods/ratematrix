##' Function calculates the log-likelihood for the multivariate Brownian motion model given a phylogenetic tree with rate regimes mapped to it, a data matrix with the trait value for the species, a vector with the phylogenetic means (root values) and a matrix (or list of matrices) with the evolutionary rate matrix for each regime.
##'
##' If more than one evolutionary rate matrix is used, then the function calculates the likelihood using the multirates prunning algorithm. Otherwise the function uses the three point algorithm to make calculations for the single regime case.
##' @title Likelihood function for the multivariate Brownian motion model
##' @param data a matrix with the data. Species names need to be provided as rownames (rownames(data) == phy$tip.label).
##' @param phy a phylogeny of the class "simmap" with the mapped regimes. The number of evolutionary rate matrices fitted to the phylogeny need to be equal to the number of regimes in phy.
##' @param root a numeric vector with the root value (phylogenetic mean).
##' @param R a matrix or a list of matrices. If 'R' is a matrix then the likelihood for a single regime is calculated. If 'R' is a list of matrices, then each matrix will be fitted to a regime in 'phy'. The R matrices will be fitted to the regimes in the same order as the columns of 'phy$mapped.edge'.
##' @return The log likelihood for the multivariate Brownian motion model.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' ## Set the prior, take a sample from it and compute the log-likelihood.
##' par.mu <- rbind( c(-10, 10), c(-10, 10) )
##' par.sd <- rbind( c(0, 10), c(0, 10) )
##' prior <- makePrior(r=2, p=2, par.mu=par.mu, par.sd=par.sd)
##' sample.prior <- samplePrior(n=1, prior=prior)
##' ## Reconstruct the variance-covariance matrix from the correlation matrix and the variances.
##' ## Note that the model has two evolutionary rate matrix regimes. So we need a list.
##' R1 <- diag(sample.prior$sd[[1]]^2) %*% sample.prior$matrix[[1]] %*% diag(sample.prior$sd[[1]]^2)
##' R2 <- diag(sample.prior$sd[[2]]^2) %*% sample.prior$matrix[[2]] %*% diag(sample.prior$sd[[2]]^2)
##' R <- list( R1, R2 )
##' ## Compute the log-likelihood.
##' likelihoodFunction(data=centrarchidae$data, phy=centrarchidae$phy.map, root=sample.prior$mu, R=R)
likelihoodFunction <- function(data, phy, root, R){
    ## Check if R is a matrix or a list.
    if( is.list(phy[[1]]) ) stop( "Function does not accept a list of phylo." )
    
    if( is.list(R) ){ ## The case with multiple regimes.
        
        if( !inherits(phy, what="simmap") ) stop( "R is a list or matrices but phy is not of type 'simmap'." )
        k <- ncol(data) ## Number of traits.

        ## Make the precalculation based on the tree. Here two blocks, depending of whether there is only one or several trees.
        ord.id <- reorder.phylo(phy, order="postorder", index.only = TRUE) ## Order for traversal.
        mapped.edge <- phy$mapped.edge[ord.id,] ## The regimes.
        ## Need to take care how to match the regimes and the R matrices.
        anc <- phy$edge[ord.id,1] ## Ancestral edges.
        des <- phy$edge[ord.id,2] ## Descendent edges.
        nodes <- unique(anc) ## The internal nodes we will traverse.

        ## Set the types for each of the nodes that are going to be visited.
        node.to.tip <- which( tabulate( anc[which(des <= length(phy$tip.label))] ) == 2 )
        node.to.node <- which( tabulate( anc[which(des > length(phy$tip.label))] ) == 2 )
        node.to.tip.node <- unique( anc )[!unique( anc ) %in% c(node.to.node, node.to.tip)]
        ## 1) nodes to tips: nodes that lead only to tips, 2) nodes to nodes: nodes that lead only to nodes, 3) nodes to tips and nodes: nodes that lead to both nodes and tips.
        names(anc) <- rep(1, times=length(anc))
        names(anc)[which(anc %in% node.to.node)] <- 2
        names(anc)[which(anc %in% node.to.tip.node)] <- 3

        loglik <- logLikPrunningMCMC(data, k, nodes, des, anc, mapped.edge, R=R, mu=root)
        return(loglik)
        
    } else{ ## The case of a single regime.
        
        ## Creates data and chain cache:
        cache.data <- list()
        cache.data$k <- ncol(data) ## Number of traits.
        cache.data$X <- data
        cache.data$n <- length(phy$tip.label)
        loglik <- logLikSingleRegime(data=cache.data, chain=NULL, phy=phy, root=root, R=R ) ## Lik start value.
        return(loglik)
    }
}
