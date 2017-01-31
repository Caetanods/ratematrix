calcNodeToNode <- function(X, k, p, des, anc, mapped.edge, R, node.id, cache){
    ## Applies for every node that only has node descendents.
    ## For description of the arguments, see code for the function 'mvloglik'.
    des.node <- des[ node.id ] ## The descendents of this node.
    key.id <- sapply(des.node, function(x) which(cache$key == x) )
    ss <- cache$X0[,key.id[1]] - cache$X0[,key.id[2]] ## The contrast for the nodes.
    ## 'node.id' always have length 2 given that the tree is bifurcating.
    ## 'Rs1' is relative to the branch 'node.id[1]' and 'Rs2' to the node 'node.id[2]'.
    Rs1 <- Reduce( "+", lapply(1:p, function(x) R[[x]] * mapped.edge[node.id[1],x]) ) ## The regime (Sigma * t) for the internal branch.
    Rs2 <- Reduce( "+", lapply(1:p, function(x) R[[x]] * mapped.edge[node.id[2],x]) ) ## The regime (Sigma * t) for the internal branch.
    Rs1 <- Rs1 + cache$V0[,,key.id[1]] ## The additional variance.
    Rs2 <- Rs2 + cache$V0[,,key.id[2]] ## The additional variance.
    Rinv <- chol2inv(chol(Rs1+Rs2))

    cache$ll <- c(cache$ll, logLikNode(ss, Rs1+Rs2, Rinv, k)) ## Incremental sum. Do not need to store.

    cache$key <- c(cache$key, anc[node.id[1]]) ## 'key' for both X0 and V0. This is the node number.
    ## Given the formula for the Harmonic mean, you just add more elements for this sum.
    cache$X0[,length(cache$key)] <- ((Rs2 %*% Rinv) %*% cache$X0[,key.id[1]]) + ((Rs1 %*% Rinv) %*% cache$X0[,key.id[2]])
    cache$V0[,,length(cache$key)] <- chol2inv(chol( chol2inv(chol(Rs1)) + chol2inv(chol(Rs2)) ))
    return( cache )
}
