##' Calculates points for ellipses.
##'
##' This funtion is not optimized to be called by the user. Function calculates the points to draw ellipses from a list of covariance matrices.\cr
##' \cr
##' If 'mat' is only a matrix then only the coordinates for the plot is returned. Otherwise, a list with the coordinates and limits calculated for
##' the series of ellipses is returned.
##' @title Calculate sample of ellipses.
##' @param mat list. List of covariance matrices.
##' @param center numeric. Coordinates for the center.
##' @param traits numeric (?). The traits (?).
##' @param sample.line numeric. Sample size to be used to calculate the ellipses.
##' @param n.points numeric. The number of points used to approximate the shape of the ellipse.
##' @return List of coordinates to plot the ellipses.
##' @author Daniel S. Caetano, \email{caetanods1@@gmail.com}
##' @noRd
getEllipseMatrix <- function(mat, center=c(0,0), traits, sample.line=50, n.points=200){

    if( is.matrix(mat) ){
        ## Make calculations for the case when one single matrix is provided.
        ## Same of the end of the function. But do not calculates the limits for the
        ##      plot.
        RR <- chol.default(mat) ## Cholesky decomposition.
        angles <- seq(0, 2*pi, length.out = n.points) ## angles for ellipse
        ## ellipse scaled with factor 1
        ell <- 1 * cbind(cos(angles), sin(angles)) %*% RR[traits,traits]
        ## center ellipse to the data centroid
        ellCtr <- sweep(ell, 2, center, "+")
        return( ellCtr )
    }

    ss <- 1:length(mat)
    ii <- sample(ss, size = sample.line)
    qq <- 1:length(ii)
    mat.list <- mat[ii]
    
    RR <- lapply(mat.list, chol.default) ## Cholesky decomposition.
    angles <- seq(0, 2*pi, length.out = n.points) ## angles for ellipse
    ## ellipse scaled with factor 1
    ell <- lapply(qq, function(x) 1 * cbind(cos(angles), sin(angles)) %*% RR[[x]][traits,traits] )
    ## center ellipse to the data centroid
    ellCtr <- lapply(qq, function(x) sweep(ell[[x]], 2, center, "+") )
    
    all.points <- do.call(rbind, ellCtr)
    xlim <- range(all.points[,1])
    ylim <- range(all.points[,2])
    return( list(limits=data.frame("xlim"=xlim, "ylim"=ylim), ellCtr=ellCtr) )
}
