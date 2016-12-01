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
