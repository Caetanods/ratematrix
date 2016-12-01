calcNodeToTip <- function(X, k, des, anc, mapped.edge, R, node.id, cache){
    ## Applies for every node that only has tip descendents.
    ## For description of the arguments, see code for the function 'mvloglik'.
    des.node <- des[ node.id ] ## The descendents of this node.
    ss <- X[des.node[1],] - X[des.node[2],] ## The contrast.
    Rs1 <- ( R[[1]] * mapped.edge[node.id[1],1] ) + ( R[[2]] * mapped.edge[node.id[1],2] )
    Rs2 <- ( R[[1]] * mapped.edge[node.id[2],1] ) + ( R[[2]] * mapped.edge[node.id[2],2] )
    Rinv <- chol2inv(chol(Rs1+Rs2))
    cache$ll <- c(cache$ll, logLikNode(ss, Rs1+Rs2, Rinv, k))
    cache$key <- c(cache$key, anc[node.id[1]]) ## 'key' for both X0 and V0. This is the node number.
    cache$X0[,length(cache$key)] <- ((Rs2 %*% Rinv) %*% X[des.node[1],]) + ((Rs1 %*% Rinv) %*% X[des.node[2],])
    cache$V0[,,length(cache$key)] <- chol2inv(chol( chol2inv(chol(Rs1)) + chol2inv(chol(Rs2)) ))
    return( cache )
}
