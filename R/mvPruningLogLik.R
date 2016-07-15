loglik.node <- function(ss, sigma_len ,sigma_len_inv, k){
    ## Internal function to calculate the log likelihood at each node.
    ## ss = contrast.
    ## sigma_len = Sigma * branch length
    ## sigma_len = solve( Sigma * branch length )
    ll <- -0.5 * ( k*log(2*pi) + log(det(sigma_len)) + (ss %*% sigma_len_inv %*% ss) )
    return(ll)
}

## Define the functions that will work in each of the branches.
## Those are going to be used for the pruning. One for each of the types of nodes.
calc_node_to_tip <- function(X, k, des, anc, mapped.edge, R, node.id, cache){
    ## Applies for every node that only has tip descendents.
    ## For description of the arguments, see code for the function 'mvloglik'.
    des.node <- des[ node.id ] ## The descendents of this node.
    ss <- X[des.node[1],] - X[des.node[2],] ## The contrast.
    Rs1 <- ( R[[1]] * mapped.edge[node.id[1],1] ) + ( R[[2]] * mapped.edge[node.id[1],2] )
    Rs2 <- ( R[[1]] * mapped.edge[node.id[2],1] ) + ( R[[2]] * mapped.edge[node.id[2],2] )
    Rinv <- chol2inv(chol(Rs1+Rs2))
    cache$ll <- c(cache$ll, loglik.node(ss, Rs1+Rs2, Rinv, k))
    cache$key <- c(cache$key, anc[node.id[1]]) ## 'key' for both X0 and V0. This is the node number.
    cache$X0[,length(cache$key)] <- ((Rs2 %*% Rinv) %*% X[des.node[1],]) + ((Rs1 %*% Rinv) %*% X[des.node[2],])
    cache$V0[,,length(cache$key)] <- chol2inv(chol( chol2inv(chol(Rs1)) + chol2inv(chol(Rs2)) ))
    return( cache )
}
calc_node_to_node <- function(X, k, des, anc, mapped.edge, R, node.id, cache){
    ## Applies for every node that only has node descendents.
    ## For description of the arguments, see code for the function 'mvloglik'.
    des.node <- des[ node.id ] ## The descendents of this node.
    key.id <- sapply(des.node, function(x) which(cache$key == x) )
    ss <- cache$X0[,key.id[1]] - cache$X0[,key.id[2]] ## The contrast for the nodes.
    Rs1 <- ( R[[1]] * mapped.edge[node.id[1],1] ) + ( R[[2]] * mapped.edge[node.id[1],2] ) ## The regime (Sigma * t) for the internal branch.
    Rs2 <- ( R[[1]] * mapped.edge[node.id[2],1] ) + ( R[[2]] * mapped.edge[node.id[2],2] ) ## The regime (Sigma * t) for the internal branch.
    Rs1 <- Rs1 + cache$V0[,,key.id[1]] ## The additional variance.
    Rs2 <- Rs2 + cache$V0[,,key.id[2]] ## The additional variance.
    Rinv <- chol2inv(chol(Rs1+Rs2))
    cache$ll <- c(cache$ll, loglik.node(ss, Rs1+Rs2, Rinv, k))
    cache$key <- c(cache$key, anc[node.id[1]]) ## 'key' for both X0 and V0. This is the node number.
    cache$X0[,length(cache$key)] <- ((Rs2 %*% Rinv) %*% cache$X0[,key.id[1]]) + ((Rs1 %*% Rinv) %*% cache$X0[,key.id[2]])
    cache$V0[,,length(cache$key)] <- chol2inv(chol( chol2inv(chol(Rs1)) + chol2inv(chol(Rs2)) ))
    return( cache )
}
calc_node_to_tip_node <- function(X, k, des, anc, mapped.edge, R, node.id, cache){
    ## Applies for every node that has BOTH tip and node descendents.
    ## For description of the arguments, see code for the function 'mvloglik'.
    des.node <- des[ node.id ] ## The descendents of this node.
    nn <- which(des.node > (min(anc)-1) )
    nd <- des.node[nn] ## This is a node.
    nd.id <- node.id[nn] ## This is the id position for this node.
    key.id <- which(cache$key == nd) ## Find the key for the node.
    tt <- which(des.node <= (min(anc)-1) )
    tip <- des.node[tt] ## This is a tip.
    tip.id <- node.id[tt] ## This is the id position for this tip.
    ss <- X[tip,] - cache$X0[,key.id] ## The contrast for the nodes.
    Rs1 <- ( R[[1]] * mapped.edge[tip.id,1] ) + ( R[[2]] * mapped.edge[tip.id,2] ) ## No additional variance. This is a tip.
    Rs2 <- ( R[[1]] * mapped.edge[nd.id,1] ) + ( R[[2]] * mapped.edge[nd.id,2] ) ## Need additional variance.
    Rs2 <- Rs2 + cache$V0[,,key.id] ## The additional variance.
    Rinv <- chol2inv(chol(Rs1+Rs2))
    cache$ll <- c(cache$ll, loglik.node(ss, Rs1+Rs2, Rinv, k))
    cache$key <- c(cache$key, anc[node.id[1]]) ## 'key' for both X0 and V0. This is the node number.
    cache$X0[,length(cache$key)] <- ((Rs2 %*% Rinv) %*% X[tip,]) + ((Rs1 %*% Rinv) %*% cache$X0[,key.id])
    cache$V0[,,length(cache$key)] <- chol2inv(chol( chol2inv(chol(Rs1)) + chol2inv(chol(Rs2)) ))
    return( cache )
}
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
##' @export
mvloglik <- function(X, phy, R, mu){

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
    node_calc <- list(calc_node_to_tip, calc_node_to_node, calc_node_to_tip_node)

    ## Traverse the tree.
    for (i in nodes) { ## Will visit all the internal nodes including the ROOT.
        node.id <- which(anc == i) ## The index for the 'des', 'anc', and 'mapped.edge (lines)'.
        type <- as.numeric( names(anc[node.id[1]]) )
        ## Next is a list of functions. Will return the updated cache for the tree traversal.
        cache <- node_calc[[type]](X, k, des, anc, mapped.edge, R, node.id, cache)
    }

    cache$ll <- c(cache$ll, loglik.node(cache$X0[,length(cache$key)] - mu, cache$V0[,,length(cache$key)], solve(cache$V0[,,length(cache$key)]), k))

    return( sum(cache$ll) )
}

loglikMCMC <- function(X, k, nodes, des, anc, mapped.edge, R, mu){
    ## This is a function to be used with the MCMC. The processing of the tree into the objects that have the information needed for traversing the tree will be computed prior to the call of this function.
    ## Create vectors to store results.
    X0 <- matrix(nrow = k, ncol = length(nodes)+1)
    V0 <- array(dim=c(k,k,length(nodes)+1))
    key <- vector(mode="numeric")
    ll <- vector(mode="numeric")
    cache <- list(X0=X0, V0=V0, key=key, ll=ll)

    ## Create the list of functions:
    ## This is the list that will have the functions to get the loglik in each case.
    ## Note that when you pass the as an argument to a function, R will make a copy of this object to a independent environment where the function will run. This means that making this list here and passing it as an argument might be the same thing.
    node_calc <- list(calc_node_to_tip, calc_node_to_node, calc_node_to_tip_node)

    ## Traverse the tree.
    for (i in nodes) { ## Will visit all the internal nodes including the ROOT.
        node.id <- which(anc == i) ## The index for the 'des', 'anc', and 'mapped.edge (lines)'.
        type <- as.numeric( names(anc[node.id[1]]) )
        ## Next is a list of functions. Will return the updated cache for the tree traversal.
        cache <- node_calc[[type]](X, k, des, anc, mapped.edge, R, node.id, cache)
    }

    cache$ll <- c(cache$ll, loglik.node(cache$X0[,length(cache$key)] - mu, cache$V0[,,length(cache$key)], solve(cache$V0[,,length(cache$key)]), k))

    return( sum(cache$ll) )
}
