##' Make plots for all elements in a list (posterior distribution) of evolutionary covariance matrices (R).
##'
##' The function make a grid of plots with ellipses representing the distribution of each variance and covariance element
##'     of a list of evolutionary covariance matrices (R). Disposition of the plot is calculated to prefer spread over the
##'     columns while minimizing number of rows.\cr
##' \cr
##' When a matrix for 'true.R' is provided, function use this matrix to create red lines and superpose them above the distribution
##'     of ellipses.
##' @title Ellipse plot of covariance matrices.
##' @param mat list. List of covariance matrices. 
##' @param mu matrix. Matrix of mean values to use as centroid for each ellipse.
##' @param burn numeric. Proportion of burnin.
##' @param thin numeric. The thinning of the posterior distribution. Same as the argument 'by' of the function 'seq'.
##' @param sample.line numeric. Number of ellipse lines to be plotted in each plot.
##' @param n.points numeric. Number of points to be used to calculate the shape of the ellipses.
##' @param center.mu logical. Whether the center of the ellipses should be guven by the argument 'mu'. Otherwise all ellipses are centered in the '0,0' coordinate.
##' @param true.R matrix. A single evolutinary to be plotted as a red line. Usually the true value of a set of simulations.
##' @param file string. When provided the functions creates a '.jpeg' figure instead of plotting in the standard device.
##' @return Plot in standard device or a '.jpeg' file. See argument 'file'.
##' @author Daniel S. Caetano, \email{caetanods1@@gmail.com}
ellipse.chart <- function(mat, mu, burn=0.5, thin=1, sample.line=50, n.points=200, center.mu=TRUE, true.R=NULL, file=NULL){
    ## Make plots for all correlations of a list of variance-covariance matrices.
    ## Number of plots equal to the unordered combination of traits in the matrix.
    ## Disposition of the plots is calculated to prefer spread in the columns while minimizing
    ##     number of rows.
    ## When the 'true.R' is provided, function use the info to make a red line as the true value.
    ##     The centroid of the true ellipse is equal to the mean of the centroids of the posterior.
    ##     The root value is badly estimated and, thus, most of the time the true centroid is out
    ##     of the plot.
    ## mat = list of covariance matrices.
    ## mu = matrix of mu values to use as centroid for each ellipse.
    ## burn = proportion of burnin.
    ## thin = thinning, same of 'by' in function 'seq'.
    ## sample.line = number of lines to appear in the final plot.
    ## n.points = Number of points to approximate the ellipse.
    ## center.mu = Logical. Whether the center of the ellipses should be given by mu. All ellipses
    ##             centered in 0,0 otherwise.
    ## true.R = The true evolutionary matrix used to generate the data.
    ## file = When a file string is provided the function creates a jpeg figure file istead of
    ##        creating the plot in the standard device.

    dd <- dim(mat[[1]])[1]

    if(center.mu == FALSE){
        mu <- matrix(0, nrow = dim(mu)[1], ncol = dim(mu)[2])
    }

    ## Special case when matrix is 2x2:
    if(dd == 2){
        ## If file has filename then save plot in jpeg format.
        if( !is.null(file) ){
            jpeg(file)
            vcv.ellipse(mat, mu, traits = c(1,2), burn, thin, sample.line, n.points)
            if( is.matrix(true.R) ){
                RR <- chol(true.R)
                angles <- seq(0, 2*pi, length.out=n.points)
                ell <- cbind(cos(angles), sin(angles)) %*% RR
                root <- colMeans(mu)
                ellCtr <- sweep(ell, 2, root, "+")
                points(ellCtr, col = "red", type = "l", lwd = 1.5)
            }
            dev.off()
            cat("Ellipse plot chart saved to ", file, "\n", sep="")
        } else{
            ## Otherwise, print to default device.
            vcv.ellipse(mat, mu, traits = c(1,2), burn, thin, sample.line, n.points)
            if( is.matrix(true.R) ){
                RR <- chol(true.R)
                angles <- seq(0, 2*pi, length.out=n.points)
                ell <- cbind(cos(angles), sin(angles)) %*% RR
                root <- colMeans(mu)
                ellCtr <- sweep(ell, 2, root, "+")
                points(ellCtr, col = "red", type = "l", lwd = 1.5)
            }
        }
        return()
    }
        
    ## Create list of combinations of traits.
    cc <- combn(x=dd, m=2, simplify=FALSE)
    chart <- length(cc)
    cat("Number of plots = ", chart, "\n", sep = "")

    ## Calculates size of the plot chart.
    ## This will optimize for less rows in function of columns.
    col <- row <- ceiling( sqrt( chart ) )
    while(col * row > chart){
        row <- row-1
        if(col * row < chart){
                row <- row+1
                break()
            }
    }

    ## If file has filename then save plot in jpeg format.
    if( !is.null(file) ){
        jpeg(file, width=480*col , height=480*row, pointsize = 6*col)
        par(mfrow = c(row,col))
        for(i in 1:chart){
            vcv.ellipse(mat, mu, traits = cc[[i]], burn, thin, sample.line, n.points)
            if( is.matrix(true.R) ){
                RR <- chol(true.R[ cc[[i]], cc[[i]] ])
                angles <- seq(0, 2*pi, length.out=n.points)
                ell <- cbind(cos(angles), sin(angles)) %*% RR
                root <- colMeans(mu)
                ellCtr <- sweep(ell, 2, root[ cc[[i]] ], "+")
                points(ellCtr, col = "red", type = "l", lwd = 1.5)
            }
        }
        dev.off()
        cat("Ellipse plot chart saved to ", file, "\n", sep="")
    } else{
        ## Otherwise, print to default device.
        par(mfrow = c(row,col))
        for(i in 1:chart){
            vcv.ellipse(mat, mu, traits = cc[[i]], burn, thin, sample.line, n.points)
            if( is.matrix(true.R) ){
                RR <- chol(true.R[ cc[[i]], cc[[i]] ])
                angles <- seq(0, 2*pi, length.out=n.points)
                ell <- cbind(cos(angles), sin(angles)) %*% RR
                root <- colMeans(mu)
                ellCtr <- sweep(ell, 2, root[ cc[[i]] ], "+")
                points(ellCtr, col = "red", type = "l", lwd = 1.5)
            }
        }
    }
}
