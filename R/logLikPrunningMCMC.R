logLikPrunningMCMC <- function(X, k, nodes, des, anc, mapped.edge, R, mu){
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
