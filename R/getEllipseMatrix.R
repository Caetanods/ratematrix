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
##' @param burn numeric. The proportion of the burnin.
##' @param thin numeric. The thinnig to be applied to the posterior distribution. Same as the 'by' argument of the 'seq' function.
##' @param sample.line numeric. Sample size to be used to calculate the ellipses.
##' @param n.points numeric. The number of points used to approximate the shape of the ellipse.
##' @return List of coordinates to plot the ellipses.
##' @author Daniel S. Caetano, \email{caetanods1@@gmail.com}
getEllipseMatrix <- function(mat, center=c(0,0), traits, burn=0.5, thin=1, sample.line=50, n.points=200){
    ## Function return objects to plot the matrix on the new visualization method.
    ## mat = list of covariance matrices. If 'mat' is only a matrix then only the coordinates
    ##     for the plot is returned. Otherwise, a list with the coordinates and limits calculated
    ##     for the series of ellipses is returned.
    ## center = vector with centroid for the two selected traits.
    ## traits = vector with two indexes for the vcv matrix indicating the traits to be plotted.
    ## burn = proportion of burnin.
    ## thin = thinning, same of 'by' in function 'seq'.
    ## sample.line = number of lines to appear in the final plot.
    ## n.points = Number of points to approximate the ellipse.

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

    ll <- length(mat)
    if(burn == 0){
        out <- round( ll / thin )
        ss <- round(seq(from = 1, to = ll, length.out = out))
    } else{
        bb <- round( ll * burn )
        out <- round( ( ll - bb ) / thin )
        ss <- round(seq(from = bb, to = ll, length.out = out))
    }
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
