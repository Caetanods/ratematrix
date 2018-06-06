##' Simulates correlated traits under a multivariate Brownian motion model. The function uses a covariance matrix (evolutionary rate matrix-R) to indicate the rates of the traits.
##'
##' This is a function derived from 'sim.corrs' in the package 'phytools'. This version has some edits to make the simulations more efficient for this particular use. For all other applications please refer to the original implementation of 'sim.corrs' in the package 'phytools' wrote by Liam Revell.
##' @title Simulates multivariate trait evolution using a Brownian motion model
##' @param tree a phylogenetic tree of 'phylo' format.
##' @param vcv a variance covariance matrix (the evolutionary rate matrix).
##' @param anc a vector of the same length as the number of traits to be simulated (same as the dimension of the 'vcv' matrix). This is used as the values for the root in the simulations. If 'NULL' then all traits have root value of 0.
##' @param internal whether to return the values simulated for the nodes in the phylogeny. If FALSE (default), then only returns the values simulated for the tips.
##' @return Returns a matrix with each trait values for the tips. Traits are distributed in the rows and tips are distributed in the columns.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @references
##' \describe{
##'   \item{}{Revell, L. J. 2012. phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution 3:217–223.}
##' }
simRatematrix <- function(tree, vcv, anc = NULL, internal = FALSE){
    if (!inherits(tree, "phylo")) 
        stop("tree should be an object of class \"phylo\".")
    if (!is.list(vcv)) {
        p <- nrow(vcv)
        if (is.null(anc)) 
            anc <- rep(0, p)
        cholvcv <- chol(vcv)
        X <- matrix(stats::rnorm(p * nrow(tree$edge), sd = rep(sqrt(tree$edge.length), 
            p)), nrow(tree$edge), p)
        X <- X %*% cholvcv
    }
    else {
        p <- nrow(vcv[[1]])
        if (is.null(anc)) 
            anc <- rep(0, p)
        if (is.null(names(vcv))) {
            names(vcv) <- colnames(tree$mapped.edge)
            message("names absent from vcv: assuming same order as $mapped.edge")
        }
        vcv <- vcv[colnames(tree$mapped.edge)]
        cholvcv <- lapply(vcv, chol)
        X <- matrix(0, nrow(tree$edge), p)
        for (i in 1:length(vcv)) {
            Y <- matrix(stats::rnorm(p * nrow(tree$edge), sd = rep(sqrt(tree$mapped.edge[, 
                i]), p)), nrow(tree$edge), p)
            X <- X + Y %*% cholvcv[[i]]
        }
    }
    Y <- array(0, dim = c(nrow(tree$edge), ncol(tree$edge), p))
    n <- length(tree$tip)
    for (i in 1:nrow(X)) {
        if (tree$edge[i, 1] == (n + 1)) 
            Y[i, 1, ] <- anc
        else Y[i, 1, ] <- Y[match(tree$edge[i, 1], tree$edge[, 
            2]), 2, ]
        Y[i, 2, ] <- Y[i, 1, ] + X[i, ]
    }
    X <- matrix(data = rbind(Y[1, 1, ], as.matrix(Y[, 2, ])), 
        length(tree$edge.length) + 1, p)
    rownames(X) <- c(n + 1, tree$edge[, 2])
    X <- as.matrix(X[as.character(1:(n + tree$Nnode)), ])
    rownames(X)[1:n] <- tree$tip.label
    if (internal == TRUE) 
        return(X)
    else return(X[1:length(tree$tip.label), ])
}
