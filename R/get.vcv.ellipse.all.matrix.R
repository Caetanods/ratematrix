##' Calculates points for ellipses.
##'
##' This funtion is not optimized to be called by the user. Function calculates the points to draw ellipses from a list of covariance matrices.
##' @title Calculate all ellipses.
##' @param mat list. List of covariance matrices.
##' @param center numeric. Coordinates for the center.
##' @param traits numeric (?). The traits (?).
##' @param n.points numeric. The number of points used to approximate the shape of the ellipse.
##' @return List of coordinates to plot the ellipses.
##' @author Daniel S. Caetano, \email{caetanods1@@gmail.com}
get.vcv.ellipse.all.matrix <- function(mat, center=c(0,0), traits, n.points=200){
    ## Function is just like 'get.vcv.ellipse.matrix'. However there is no sample of the matrices.
    ## Points are calculated for all the matrices.

    mat.list <- mat
    
    RR <- lapply(mat.list, chol.default) ## Cholesky decomposition.
    angles <- seq(0, 2*pi, length.out = n.points) ## angles for ellipse
    ## ellipse scaled with factor 1
    ell <- lapply(RR, function(x) 1 * cbind(cos(angles), sin(angles)) %*% x[traits,traits] )
    ## center ellipse to the data centroid
    ellCtr <- lapply(ell, function(x) sweep(x, 2, center, "+") )
    
    all.points <- do.call(rbind, ellCtr)
    xlim <- range(all.points[,1])
    ylim <- range(all.points[,2])
    return( list(limits=data.frame("xlim"=xlim, "ylim"=ylim), ellCtr=ellCtr) )
}
