calcNodeToTipNode <- function(X, k, p, des, anc, mapped.edge, R, node.id, cache){
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
    ## 'node.id' always have length 2 given that the tree is bifurcating.
    ## 'Rs1' is relative to the branch leading to a tip 'tip.id' and 'Rs2' to the branch leading to a node 'nd.id'.
    Rs1 <- Reduce( "+", lapply(1:p, function(x) R[[x]] * mapped.edge[tip.id,x]) ) ## The regime (Sigma * t) for the internal branch. No additional variance. This is a tip.
    Rs2 <- Reduce( "+", lapply(1:p, function(x) R[[x]] * mapped.edge[nd.id,x]) ) ## The regime (Sigma * t) for the internal branch.
    Rs2 <- Rs2 + cache$V0[,,key.id] ## The additional variance for the node.
    Rinv <- chol2inv(chol(Rs1+Rs2))
    cache$ll <- c(cache$ll, logLikNode(ss, Rs1+Rs2, Rinv, k))
    cache$key <- c(cache$key, anc[node.id[1]]) ## 'key' for both X0 and V0. This is the node number.
    cache$X0[,length(cache$key)] <- ((Rs2 %*% Rinv) %*% X[tip,]) + ((Rs1 %*% Rinv) %*% cache$X0[,key.id])
    cache$V0[,,length(cache$key)] <- chol2inv(chol( chol2inv(chol(Rs1)) + chol2inv(chol(Rs2)) ))
    return( cache )
}
