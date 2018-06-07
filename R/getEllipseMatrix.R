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
getEllipseMatrix <- function(mat, center=c(0,0), traits, sample.line=NULL, n.points=200){

    if( is.matrix(mat) ){
        ellCtr <- ellipse::ellipse(x = mat, centre = center, which = traits, n.points = n.points)
        return( ellCtr )
    } else{
        
        if( !is.null( sample.line ) ){
            ## Take a sample from the lines.
            ss <- 1:length(mat)
            ii <- sample(ss, size = sample.line)
            qq <- 1:length(ii)
            mat.list <- mat[ii]
            ellCtr <- lapply(qq, function(x) ellipse::ellipse(x = mat[[x]], centre = center, which = traits, n.points = n.points) )
        }

        if( is.null( sample.line ) ){
            ## Use all the matrices in the list
            ss <- 1:length(mat)
            ellCtr <- lapply(ss, function(x) ellipse::ellipse(x = mat[[x]], centre = center, which = traits, n.points = n.points) )
        }   
        
        all.points <- do.call(rbind, ellCtr)
        xlim <- range(all.points[,1])
        ylim <- range(all.points[,2])
        return( list(limits=data.frame("xlim"=xlim, "ylim"=ylim), ellCtr=ellCtr) )
    }
}
