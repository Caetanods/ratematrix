##' Function perform marginal maximum likelihood estimate for the state at the nodes and map the branches to the respective regimes. Transitions between regimes are made in the mid point of the branches.
##'
##' If the tree shows politomies the function will solve the politomies randomly and add a branch of length equal to l/100, where l is the smallest branch length observed in the tree. This procedure will guarantee that there will be no zero branch length. The arbitrary branch added to the solved politomies is much smaller than the ones observed in the empirical tree, as a result it is very unlikely that the ancestral state reconstruction will be biased.
##'
##' This analysis does not substitute the simulations of transitions between states provided by a stochastic map. We suggest the user to contrast results with a distribution of stochastic maps and, if possible, make analyses based on this map and samples from the stochastic map to incorporate uncertainty in ancestral state estimates.
##' @title Map ML ancestral estimates to the tree.
##' @param tree A phylogenetic tree of the 'phylo' format.
##' @param state A vector with states to map the regimes to the tree. Names need to be species and match tree$tip.label .
##' @return A 'Simmap' format tree. This can be used in subsequent functions such as multi.R.iwish.mcmc.
##' @export
makeMLEmap <- function(tree, state){
    
    if( !is.binary.tree(tree) ){
        tree <- multi2di( tree ) ## This will generate zero branches.
        zero <- which(tree$edge.length == 0)
        tiny <- min(tree$edge.length[-zero])/100 ## 100 times smaller than the smallest branch.
        tree$edge.length[zero] <- tiny ## Very tiny branch. But not zero.
    }
    anc <- ace(x=state, phy=tree, type="discrete", marginal=TRUE)
    ## This will assume that the node state is equal to the larger likelihood. Strict.
    ## Need to make it clear in the function and give the user another option too.
    node.state <- as.numeric( anc$lik.anc[,1] < anc$lik.anc[,2] )

    state.code <- colnames(anc$lik.anc)
    node <- (length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode)
    map <- rep(0, times=length(tree$edge.length) )
    for(i in 1:length(tree$edge.length)) map[i] <- node.state[ which(node == tree$edge[i,1]) ]
    name.map <- map
    ## Line bellow assumes that 'map' are sequential numbers starting from 0.
    for(i in unique(name.map) ) name.map[which(name.map == i)] <- state.code[i+1]
    value.map <- tree$edge.length
    names(value.map) <- name.map
    maps <- list()
    for(i in 1:length(value.map)) maps[[i]] <- value.map[i]
    mapped.edge <- cbind( as.numeric(!map), map ) * tree$edge.length
    
    ## This section will set the change for one state to the other to occur in the middle of the edge.
    nodes <- seq(from=Ntip(tree)+1, length.out = length(node.state) )
    for(j in 1:length(nodes) ){
        desc <- geiger:::.get.desc.of.node(nodes[j], tree)
        for(i in desc){
            if( i <= Ntip(tree) ) next
            st.desc <- node.state[ which( i == nodes ) ]
            st.node <- node.state[ j ]
            if( st.desc == st.node ) next
            node.change <- nodes[ which( i == nodes ) ]
            edge.change <- which( node.change == tree$edge[,2] ) ## branch to change.
            edge.mid <- sum( mapped.edge[ edge.change,] )/2
            mapped.edge[edge.change, c(st.node+1, st.desc+1)] <- edge.mid
            maps.mid <- c(edge.mid, edge.mid)
            names(maps.mid) <- state.code[ c(st.node+1,st.desc+1) ]
            maps[[edge.change]] <- maps.mid
        }
    }
    
    colnames(mapped.edge) <- state.code
    tree$mapped.edge <- mapped.edge
    tree$maps <- maps
    class(tree) <- c("simmap", "phylo")
    return(tree)
}
