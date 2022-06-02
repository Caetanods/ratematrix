##' Returns the log-likelihood for the multivariate Brownian motion model with 1 or more rate regimes mapped to the tree.
##'
##' If two or more rate regimes are mapped to the phylogenetic tree, then the function calculates the likelihood using the new prunning algorithm adapted to fit multiple rate regimes. The prunning algorithm is implemented in C++ using 'Rcpp' and 'RcppArmadillo'. Otherwise the function uses the three point algorithm (Ho and Ané, 2014) to make calculations for the single regime case.
##' @title Likelihood function for the multivariate Brownian motion model
##' @param data a matrix with the data. Species names need to be provided as rownames (rownames(data) == phy$tip.label).
##' @param phy a phylogeny of the class "simmap" with the mapped regimes or "phylo" for a single rate model.
##' @param root a numeric vector with the root value (phylogenetic mean).
##' @param R a matrix or a list of matrices. If 'R' is a matrix then the likelihood for a single regime is calculated. If 'R' is a list of matrices, then each matrix will be fitted to a regime in 'phy' and the length of the list need to match the number of regimes fitted to the tree.
##' @return The log likelihood for the multivariate Brownian motion model.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @references Ho, L. S. T. and Ané, C. (2014). "A linear-time algorithm for Gaussian and non-Gaussian trait evolution models". Systematic Biology *63*(3):397-408.
##' @examples
##' \donttest{
##' data( centrarchidae )
##' root <- colMeans( centrarchidae$data )
##' Rlist <- list( rbind(c(0.5, 0.1),c(0.1,0.5)), rbind(c(0.5, 0),c(0,0.5)) )
##' likelihoodFunction(data = centrarchidae$data, phy = centrarchidae$phy.map, root = root
##'                    , R = Rlist)
##' ## Get the likelihood for a single regime model:
##' phy.single <- mergeSimmap(phy = centrarchidae$phy.map, drop.regimes = TRUE)
##' Rsingle <- rbind(c(0.5, 0.1),c(0.1,0.5))
##' likelihoodFunction(data = centrarchidae$data, phy = phy.single, root = root, R = Rsingle)
##' }
likelihoodFunction <- function(data, phy, root, R){

    ## Need to transpose the data matrix.
    X <- t( as.matrix(data) )
    k <- nrow(X) ## Number of traits.

    ## Correct the format for the root.
    mu <- as.numeric( root )
    
    ## Check if R is a matrix or a list.
    phy_type <- check_phy_list(phy)
    if( phy_type ) stop( "Function does not accept a list of phylo." )
    
    if( is.list(R) ){ ## The case with multiple regimes.

        ## Correct format to matrix:
        R <- lapply(R, function(x) as.matrix(x) )
        p <- length(R) ## The number of regimes fitted to the tree.
        
        if( !inherits(phy, what="simmap") ) stop( "R is a list or matrices but phy is not of type 'simmap'." )
        k <- ncol(data) ## Number of traits.
        Rarray <- array(dim=c(k,k,p))
        for( i in 1:p){
            Rarray[,,i] <- R[[i]]
        }
        
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
        names_anc <- names(anc)
        names_anc <- as.numeric( names_anc ) ## Type integer.

        lik <- logLikPrunningMCMC_C(X=X, k=k, p=p, nodes=nodes, des=des, anc=anc,
                                    names_anc=names_anc
                                  , mapped_edge=mapped.edge, R=Rarray, mu=mu)
        return(lik)
        
    } else{ ## The case of a single regime.

        ## Correct the format for the matrix:
        R <- as.matrix( R )
        
        ## Creates data and chain cache:
        cache.data <- list()
        cache.data$k <- ncol(data) ## Number of traits.
        cache.data$X <- as.matrix(data)
        cache.data$n <- length(phy$tip.label)
        loglik <- logLikSingleRegime(data=cache.data, chain=NULL, phy=phy, root=root, R=R ) ## Lik start value.
        return(loglik)
    }
}
