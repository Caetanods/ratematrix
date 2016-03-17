##' This creates and draws a ellipse from a covariance matrix.
##'
##' Creates ellipse plots from a covariance matrix.
##' @title Ellipse plots.
##' @param mat list. A list of covariance matrices. The posterior distribution of evolutionary rate matrices fitted to the tree.
##' @param mu matrix. Matrix of values to be used as centroid to the ellipses. This allows the user to supply a different centroid for each trait.
##' @param traits vector. This is the index of the traits to be plotted. It is a length 2 numeric vector, the first indicates the row of the matrix and the second indicates the column.
##' @param burn numeric. Proportion of burnin to be excluded from the chain.
##' @param thin numeric. Thinning of the posterior chain. This is the same as the argument 'by' of the function 'seq'.
##' @param sample.line numeric. The number of lines to appear in the final plot. This is equal to the number of samples taken from the posterior chain.
##' @param n.points numeric. Number of points used to calculate the ellipse. 
##' @return Plot the ellipse for a combination of two traits. This function does not plot the combination of all plots in a grid. For this see 'make.grid.plot'.
vcv.ellipse <- function(mat, mu, traits=c(1,2), burn=0.5, thin=1, sample.line=50, n.points=200){
    ## Drawing an ellipse from a variance-covariance matrix.
    ## mat = list of covariance matrices.
    ## mu = matrix of mu values to use as centroid for each ellipse.
    ## traits = vector with two indexes for the vcv matrix indicating the traits to be plotted.
    ## burn = proportion of burnin.
    ## thin = thinning, same of 'by' in function 'seq'.
    ## sample.line = number of lines to appear in the final plot.
    ## n.points = Number of points to approximate the ellipse.

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
    mu.mm <- mu[ii,]
    
    RR <- lapply(mat.list, chol.default) ## Cholesky decomposition.
    angles <- seq(0, 2*pi, length.out = n.points) ## angles for ellipse
    ## ellipse scaled with factor 1
    ## 'traits' define which traits to be plotted.
    ell <- lapply(qq, function(x) 1 * cbind(cos(angles), sin(angles)) %*% RR[[x]][traits,traits] )
    ## center ellipse to the data centroid
    ## 'traits' define which mu values to be used.
    ellCtr <- lapply(qq, function(x) sweep(ell[[x]], 2, mu.mm[x,traits], "+") )
    
    all.points <- do.call(rbind, ellCtr)
    xlim <- range(all.points[,1])
    ylim <- range(all.points[,2])
    plot(axes = FALSE,x = 0, y = 0, type = "n", xlim = xlim, ylim = ylim
       , asp = 1, xlab = NA, ylab = NA, pty = "s", main = paste( "Ellipse plot -", paste(traits, collapse = ","), sep = " ") )
    axis(1, padj = -0.8, c(roundoff(xlim[1], 1, 1), round(mean(xlim),  digits = 1),
        roundoff(xlim[2], 1, -1)))
    axis(2, padj = 0.5, c(roundoff(ylim[1], 1, 1), round(mean(ylim),  digits = 1),
        roundoff(ylim[2], 1, -1)))
    box()
    title(xlab = paste("Trait",traits[1],sep=" "), ylab = paste("Trait",traits[2],sep=" ")
        , cex.lab = 1, line = 1.5)
    invisible( lapply(ellCtr, points, col = "grey35", type = "l", lwd = 0.5) )
}
