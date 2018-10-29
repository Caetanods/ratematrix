##' Make stochastic map simulations conditioned on a Markov matrix 'Q' and a vector of root probabilities 'pi'.
##'
##' This function is a simplification of Revell's 'phytools::make.simmap' function. Here the stochastic mapping is performed conditioned on a given Markov matrix and a vector of probabilities for the root node. This allows users to fit the Mk model using any preferred method and use this function to perform the stochastic mapping on the tree.
##'
##' The prior probabilities at the root can be set to "equal" (i.e., all states have the same probability to be observed at the root) or to "madfitz" (i.e., state probabilities follow the likelihood of the Mk model).
##' 
##' The argument 'max_nshifts' controls the max number of state changes in any given branch of the phylogeny. This sets the size of the "buffer" that records the events that happen on the branches. It has no influence in the model, it is only a contraint of the fast implementation of the stochastic mapping algorithm. Set this value to a high enough number (i.e., more changes that can possibly happen at any given branch).
##' 
##' The reduced time is accomplished by using compiled code to perfom simulations ( C++ ). All calculations are the same as Revell's original function.
##' @title Fast implementation of stochastic mapping.
##' @param tree a phylogenetic tree of class 'phylo'.
##' @param x a named vector with the states observed at the tips of the tree.
##' @param Q a Markov transition matrix for the Markov Model. This needs to be provided and the user can estimate such matrix from the observed data using any of a multitude of methods.
##' @param pi one of 'equal' or 'madfitz'.
##' @param nsim number of stochastic mappings to be performed conditioned on Q. NOT IMPLEMENTED YET!
##' @param mc.cores same as in 'parallel::mclapply'. This is used to make multiple simulations (controlled with the argument 'nsim') by calling 'parallel::mclapply'.
##' @param max_nshifts allocate the max number of events in any given branch. See 'Details'.
##' @param silence if function should stop printing messages. This will also stop checks for data format and some informative errors.
##' @return A stochastic mapped phylogeny of class 'simmap'.
##' @author Daniel Caetano
##' @export
##' @importFrom ape reorder.phylo Nnode Ntip
##' 
fastSimmap <- function(tree, x, Q, pi = "equal" , nsim = 1, mc.cores = 1, max_nshifts = 100, silence = FALSE){
    if( !silence ){
        ## Make essential tests on the data:
        if( inherits(tree,"multiPhylo") ) stop( "Don't work with 'multiPhylo'. See examples." )
        if( !inherits(tree,"phylo") ) stop("tree should be object of class \"phylo\".")
        if( is.null( colnames( Q ) ) ) stop( "Q needs colnames equal to the states in the data." )
        if( any( !colnames( Q ) %in% x ) ) stop(" colnames of Q need to be the states in the data." )
        if( !ncol( Q ) == length( unique(x) ) ) stop(" number of states in x and Q need to match." )
        if( any( !unique(x) %in% colnames(Q) ) ) stop(" all states in the data need to ne present in Q." )
        if( any(!rowSums( Q ) == 0) ) stop( " rowSums(Q) need to be 0 " )
        cat("MAKE SURE THE DATA AND THE TREE MATCH! \n")
    }

    ## NEED TO IMPLEMENT THIS:
    if( nsim > 1 ) stop( "Sorry. Multiple simmaps coming soon!" )
    
    ## Compute the elements to use the function:
    prun.tree <- reorder.phylo(x = tree, order = "postorder")
    n_nodes <- Nnode(prun.tree)
    n_tips <- Ntip(prun.tree)
    n_states <- length(unique(x))
    edge_len <- prun.tree$edge.length
    edge_mat <- prun.tree$edge
    parents <- unique( prun.tree$edge[,1] )
    ## The order of the columns in X need to be the same as the columns in Q.
    X <- makeDataTips(X = x, states = colnames(Q) )
    root_node <- n_tips + 1
    root_type <- switch(pi, "equal" = 0, "madfitz" = 1)

    ## Get the maps list:
    get.maps <- makeSimmapMaps(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                          , edge_len=edge_len, edge_mat=edge_mat
                                          , parents=parents, X=X, Q=Q, root_node=root_node
                                          , root_type=root_type, max_nshifts = max_nshifts)

    ## Process the results
    nrow <- nrow( get.maps ) / 2
    states <- colnames(Q)
    states.mat <- get.maps[1:nrow,]
    states.mat.ref <- get.maps[1:nrow,]
    ## Translate the states.
    states.mat[states.mat.ref == length(states)] <- NA ## Need to delete the dummies first.
    for( i in 1:length(states) ){
        states.mat[states.mat.ref == i-1] <- states[i]
    }

    ## Now can create the 'maps' object with the correct state names.
    maps <- lapply(1:nrow, function(x) setNames(get.maps[x+nrow,], states.mat[x,])[get.maps[x+nrow,] > 0] )

    ## Get the mapped edge.
    mapped.edge <- matrix(data = 0, nrow = length(maps), ncol = ncol(Q))
    colnames(mapped.edge) <- colnames(Q)
    for( j in 1:length(maps) ){
        for( i in 1:length(states)){
            mapped.edge[j,i] <- sum( maps[[j]][ names(maps[[j]]) == states[i] ] )
        }
    }

    ## Make the 'simmap' object and exit.
    prun.tree$maps <- maps
    prun.tree$mapped.edge <- mapped.edge
    rownames( prun.tree$mapped.edge ) <- apply(prun.tree$edge, 1, function(x) paste(x, collapse=",") )
    prun.tree$Q <- Q
    ## Get the log likelihood for the Q. This is an element of the simmap object.
    prun.tree$logL <- logLikMk_C(n_nodes=n_nodes, n_tips=n_tips
                                            , n_states=n_states, edge_len=edge_len
                                            , edge_mat=edge_mat, parents=parents, X=X, Q=Q
                                            , root_node=root_node, root_type=root_type)
    
    ## Cannot append the class here. Strangely the class need to be in this particular order.
    class(prun.tree) <- c("simmap", "phylo")

    return( prun.tree )
}
