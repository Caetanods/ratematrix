calcNodeToTip <- function(X, k, p, des, anc, mapped.edge, R, node.id, cache){
    ## Applies for every node that only has tip descendents.
    ## For description of the arguments, see code for the function 'mvloglik'.
    des.node <- des[ node.id ] ## The descendents of this node.
    ss <- X[des.node[1],] - X[des.node[2],] ## The contrast.
    ## 'node.id' always have length 2 given that the tree is bifurcating.
    ## 'Rs1' is relative to the branch 'node.id[1]' and 'Rs2' to the node 'node.id[2]'.
    Rs1 <- Reduce( "+", lapply(1:p, function(x) R[[x]] * mapped.edge[node.id[1],x]) ) ## The regime (Sigma * t) for the internal branch.
    Rs2 <- Reduce( "+", lapply(1:p, function(x) R[[x]] * mapped.edge[node.id[2],x]) ) ## The regime (Sigma * t) for the internal branch.
    Rinv <- chol2inv(chol(Rs1+Rs2))

    cache$ll <- c(cache$ll, logLikNode(ss, Rs1+Rs2, Rinv, k))

    cache$key <- c(cache$key, anc[node.id[1]]) ## 'key' for both X0 and V0. This is the node number.
    cache$X0[,length(cache$key)] <- ((Rs2 %*% Rinv) %*% X[des.node[1],]) + ((Rs1 %*% Rinv) %*% X[des.node[2],])
    cache$V0[,,length(cache$key)] <- chol2inv(chol( chol2inv(chol(Rs1)) + chol2inv(chol(Rs2)) ))
    return( cache )
}
