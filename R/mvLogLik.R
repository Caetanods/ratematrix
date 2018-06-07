##' Calculates the log likelihood of a tree, rates and phylogenetic mean. Applies to both a single and multiple regimes. Right now it only works with a constant regime or with a regime with two rate matrices.
##'
##' This function uses the pruning algorithm. This avoids the calculation of the inverse and the determinant of large matrices. Also it makes calculation much more stable, so the likelihood can be calculated even for very large matrices. This function is to be used outside of the MCMC, since one can precaulculate some quantities that are constant given the tree and the data.
##' The regimes in the 'simmap' tree need to be named to assume the proper match between the rate matrices and the tree.
##' @title Multivariate Brownian motion log likelihood
##' @param X A matrix with the trait data.
##' @param phy The phylogeny with "named" regimes in the 'simmap' format. See 'phytools:make.simmap' function.
##' @param R A named list of rate matrices. Names must match the regimes of the 'simmap' tree.
##' @param mu A vector of phylogenetic means.
##' @return Returns the log likelihood.
##' @noRd
mvLogLik <- function(X, phy, R, mu){

    k <- ncol(X) ## Number of traits.
    ord.id <- reorder.phylo(phy, order="postorder", index.only = TRUE)
    mapped.edge <- phy$mapped.edge[ord.id,] ## The regimes. Need to take care how to match the regimes and the R matrices. Right it seems inverted.
    anc <- phy$edge[ord.id,1] ## Ancestral edges.
    des <- phy$edge[ord.id,2] ## Descendent edges.
    nodes <- unique(anc) ## The internal nodes we will traverse.

    node.to.tip <- which( tabulate( anc[which(des <= length(phy$tip.label))] ) == 2 )
    node.to.node <- which( tabulate( anc[which(des > length(phy$tip.label))] ) == 2 )
    node.to.tip.node <- unique( anc )[!unique( anc ) %in% c(node.to.node, node.to.tip)]
    ## 1) nodes to tips, 2) nodes to nodes, 3) nodes to tips and nodes.
    names(anc) <- rep(1, times=length(anc))
    names(anc)[which(anc %in% node.to.node)] <- 2
    names(anc)[which(anc %in% node.to.tip.node)] <- 3

    ## Create vectors to store results.
    X0 <- matrix(nrow = k, ncol = length(nodes)+1)
    V0 <- array(dim=c(k,k,length(nodes)+1))
    key <- vector(mode="numeric")
    ll <- vector(mode="numeric")
    cache <- list(X0=X0, V0=V0, key=key, ll=ll)

    ## Create the list of functions:
    ## This is the list that will have the functions to get the loglik in each case.
    node_calc <- list(calcNodeToTip, calcNodeToNode, calcNodeToTipNode)

    ## Traverse the tree.
    for (i in nodes) { ## Will visit all the internal nodes including the ROOT.
        node.id <- which(anc == i) ## The index for the 'des', 'anc', and 'mapped.edge (lines)'.
        type <- as.numeric( names(anc[node.id[1]]) )
        ## Next is a list of functions. Will return the updated cache for the tree traversal.
        cache <- node_calc[[type]](X, k, des, anc, mapped.edge, R, node.id, cache)
    }

    cache$ll <- c(cache$ll, logLikNode(cache$X0[,length(cache$key)] - mu, cache$V0[,,length(cache$key)], solve(cache$V0[,,length(cache$key)]), k))

    return( sum(cache$ll) )
}
