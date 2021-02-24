##' Make a stochastic map simulation conditioned on a Markov matrix 'Q' and a vector of root probabilities 'pi'.
##'
##' This function is a simplification of Revell's 'phytools::make.simmap' function. Here the stochastic mapping is performed conditioned on a given Markov matrix and a vector of probabilities for the root node. This allows users to fit the Mk model using any preferred method and use this function to perform the stochastic mapping on the tree.
##'
##' The function returns a single stochastic map in the 'simmap' format. In order to get multiple simulations, simply call this function multiple times using 'lapply', see 'Examples'.
##'
##' The prior probabilities at the root can be set to "equal" (i.e., all states have the same probability to be observed at the root) or to "madfitz" (i.e., state probabilities follow the likelihood of the Mk model).
##' 
##' The argument 'max_nshifts' controls the size of the "memory buffer" that records the number of state changes in any given branch of the phylogeny. It DOES NOT influence the outcome of the stochastic character map simulations. Set this value to a high enough number (i.e., more changes that can possibly happen at any given branch). If the limit is reached the function will print a message and return a value of 0.0 instead of the stochastic map. If that happens, simply increase the number of 'max_nshifts' and run again. This is only a limitation of the computer algorithm used to speed up the simulation and DOES NOT affect the results in any way.
##' 
##' The reduced time is accomplished by using compiled code to perfom simulations ( C++ ). All calculations are the same as Revell's original function.
##' 
##' If some of the states in the transition matrix "Q" are not present among the observed tips of the phylogeny the function will return some warning messages. The stochastic mapping will work properly however. Please check that ALL states among the tips of the phylogeny are represented on some of the columns (and rows) of the transition matrix "Q".
##' @title Fast implementation of stochastic mapping.
##' @param tree a phylogenetic tree of class 'phylo'.
##' @param x a named vector with the states observed at the tips of the tree.
##' @param Q a Markov transition matrix for the Markov Model. This needs to be provided and the user can estimate such matrix from the observed data using any of a multitude of methods.
##' @param pi one of 'equal' or 'madfitz'.
##' @param mc.cores same as in 'parallel::mclapply'. This is used to make multiple simulations (controlled with the argument 'nsim') by calling 'parallel::mclapply'.
##' @param max_nshifts allocate the max number of events in any given branch. See 'Details'.
##' @param silence if function should skip data format checks tests and stop printing messages.
##' @return A stochastic mapped phylogeny of class 'simmap' or a value of 0 if 'max_nshifts' is reached. Please see 'Details'.
##' @author Daniel Caetano
##' @export
##' @importFrom ape reorder.phylo Nnode Ntip
##' @examples
##' \donttest{
##' ## Load data
##' data(anoles)
##' area <- setNames(object = as.character(anoles$data$area), nm = rownames(anoles$data))
##' phy <- mergeSimmap(phy = anoles$phy[[1]], drop.regimes = TRUE)
##' ## Define a transition matrix. This can be estimated using MLE or MCMC.
##' ## Building one as an example.
##' Q <- matrix(0.0007, nrow = 2, ncol = 2)
##' diag(Q) <- diag(Q) * -1
##' colnames(Q) <- unique(area)
##' ## Generate 10 stochastic mappings using lapply:
##' maps <- lapply(1:10, function(x) fastSimmap(tree = phy, x = area, Q = Q))
##' ## Now using a simple for loop.
##' maps <- vector(mode = "list", length = 10)
##' for( i in 1:10 ) maps[[i]] <- fastSimmap(tree = phy, x = area, Q = Q)
##' }
fastSimmap <- function(tree, x, Q, pi = "equal", mc.cores = 1, max_nshifts = 200, silence = FALSE){
    if( !silence ){
        ## Make essential tests on the data:
        if( inherits(tree,"multiPhylo") ) stop( "Don't work with 'multiPhylo'. See examples." )
        if( !inherits(tree,"phylo") ) stop("Tree should be object of class \"phylo\".")
        if( is.null( colnames( Q ) ) ) stop( "Q needs colnames equal to the states in the data." )
        if( any( !unique(x) %in% colnames(Q) ) ) stop("All states in the data need to ne present in Q." )
        ## if( any( !colnames( Q ) %in% x ) ) stop(" colnames of Q need to be the states in the data." )
        if( any( !colnames( Q ) %in% x ) ) warning( "Some states in Q are not present among the tips. See Details." )
        if( any(abs(rowSums(Q)) > 1e-8) ) stop( "Rows of Q need to sum to 0." )
        if( is.null( names(x) ) ) stop("Data need to have names matching the tips of the phylogeny.")
        if( !Ntip(tree) == length(x) ) stop("Number of species in the data and the tree must be the same.")
        if( any( is.na(x) ) ) stop("Data contains NAs.")
    }
    
    ## Always make sure that the tree and the data match:
    mm <- match(tree$tip.label, names(x))
    x <- x[mm]
    
    ## Check for zero length branches:
    zb <- any( tree$edge.length == 0 )
    ## These quantities need to be global variables.
    short_branch <- min( tree$edge.length[tree$edge.length > min(tree$edge.length)] )
    almost_zero <- short_branch / 100000 ## A very small fraction of the shortest branch!
    
    if( zb ){
        ## A fix when the phylogeny includes branch lengths of zero length.
        zero_branches <- tree$edge.length == 0
        tree$edge.length[ tree$edge.length == 0 ] <- almost_zero
    }
    
    ## Compute the elements to use the function:
    prun.tree <- reorder.phylo(x = tree, order = "postorder")
    n_nodes <- Nnode(prun.tree)
    n_tips <- Ntip(prun.tree)
    n_states <- ncol(Q)
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
    if( nrow(get.maps) == 2 ){
        ## Number of maximum events on the branch has been reached.
        ## Will just return a value of 0.
        return( 0.0 )
    }

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
    ## Some application in phytools try to get the states for the model from the rownames of Q:
    rownames( Q ) <- colnames( Q )
    prun.tree$Q <- Q
    ## Get the log likelihood for the Q. This is an element of the simmap object.
    prun.tree$logL <- logLikMk_C(n_nodes=n_nodes, n_tips=n_tips
                                            , n_states=n_states, edge_len=edge_len
                                            , edge_mat=edge_mat, parents=parents, X=X, Q=Q
                                            , root_node=root_node, root_type=root_type)
    
    ## Cannot append the class here. Strangely the class need to be in this particular order.
    class(prun.tree) <- c("simmap", "phylo")
    
    if( zb ){
        ## Before returning the tree we need to return the zero length branches to their original state.
        prun.tree$edge.length[ prun.tree$edge.length < short_branch ] <- 0.0
        for( i in 1:length(prun.tree$maps) ){
            if( length( prun.tree$maps[[i]] ) == 1 ){
                if( prun.tree$maps[[i]] < short_branch ){
                    prun.tree$maps[[i]] <- setNames(object = 0.0, nm = names(prun.tree$maps[[i]]))
                }
            }
        }
        prun.tree$mapped.edge[ prun.tree$mapped.edge == almost_zero ] <- 0.0
    }

    return( prun.tree )
}
