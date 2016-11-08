##' Make a grid plot with the posterior distribution of the evolutionary rate matrix (R). Function works with one or two matrices fitted to the same tree and will superpose the results.
##'
##' The upper-tri plots are the posterior density for the off-diagonal elements of the
##'     R evolutionary matrix. The diagonal plots are the posterior elements of the evolutionary
##'     rates for each of the traits. The lower-tri plots are ellipse plots for the rates of correlated
##'     evolution.\cr
##' \cr
##' Two posterior distributions of matrices can be plotted in the same grid if a distribution for
##'     'mat2' is also provided (optional). Otherwise only 'mat1' is plotted.\cr
##' \cr
##' If a second distribution of matrices is provided ('mat2') the function will check whether the number
##'      of matrices in each list is the same. If not, the first elements of the lengthier list will be pruned
##'      to match the length of the shorter list.\cr
##' \cr
##' Function use 'par()' to set graphical parameters. It may take some time to create all the plots if
##'      the number of traits is large. Thinning the posterior distribution of matrices before plotting is
##'      recommended.
##' @title Plot posterior distribution of rate matrices.
##' @param chain The MCMC chain loaded from the files.
##' @param p A vector with the regimes to be plotted.
##' @param colors A vector with colors. Same length as p.
##' @param alphaOff Transparency of the off-diagonals plots. Numeric between 0 and 1.
##' @param alphaDiag Transparency of the diagonals plots. Numeric between 0 and 1.
##' @param alphaEll Transparency of the ellipses lines. Numeric between 0 and 1.
##' @param ell.wd The width of the ellipses lines.
##' @param point.matrix Optional. A list with matrices of length equal to p. This will be plotted as lines in the grid.plot.
##' @param point.color Optional. A vector of colors with length equal to p. This is a alternative color vector to the 'colors' argument to be used with the point matrices. If not provided then the colors of the point matrices will be equal to 'colors'.
##' @param point.wd Optional. The width of the lines for the vertical and ellipses of the point matrices.
##' @param leg string. Legend for the traits. Vector need to have length equal to the dimension of the R matrix.
##' @param l.cex numeric. 'cex' parameter for 'leg'. See 'help(par)' for more information on 'cex'.
##' @param hpd numeric. Set the proportion of the highest posterior density (HPD) to be highlighted in the plot. If set .95 the distributions and the ellipses within the 95\% HPD will be highlighted in the plot. Default value does not show highlighted region.
##' @param show.zero logical. Whether to plot a 'blue' vertical line at 0.
##' @param set.xlim numeric. Two elements vector to set the 'xlim' [see 'help(hist)'] for the density plots manually. If 'NULL', the default value, then the limits are calculated from the data.
##' @param n.lines numeric. Number of lines to be displayed in the ellipsed plots.
##' @return Plot a grid of charts.
##' @export
make.grid.plot <- function(chain, p, colors, alphaOff=1, alphaDiag=1, alphaEll=1, ell.wd=0.5, point.matrix=NULL, point.color=NULL, point.wd=0.5, leg=NULL, l.cex=0.7, hpd=100, show.zero=FALSE, set.xlim=NULL, n.lines=50){
    
    ## First do a batch of tests:
    check.mat <- vector()
    check.length <- vector()
    for(i in p){
        check.mat[i] <- ncol( chain$matrix[[i]][[1]] ) ## First element of each regime.
        check.length[i] <- length( chain$matrix[[i]] ) ## The length of each chain regime.
    }
    equal.size <- sapply(2:length(p), function(x) check.mat[1] == check.mat[x] )
    equal.length <- sapply(2:length(p), function(x) check.length[1] == check.length[x] )

    if( !sum(equal.size) == (length(p)-1) ) stop("Matrix regimes do not have the same number of traits.")
    if( !sum(equal.length) == (length(p)-1) ) stop("Chain for regimes do not have the same length.")
    
    ll <- check.length[1]
    dd <- check.mat[1]

    ## Create default titles:
    ## This need to be changed to get the name of the traits from the MCMC chain.
    ## Also need to be able to relate the color of the plot with the regime of the R matrix.
    if(is.null(leg)){
        leg <- paste("trait_", 1:dd, sep="")
    }

########################################################################################
    ## BLOCK FOR HIGHLIGHT THE HPD:

    if(hpd < 100){
        ## Calculate the probabilities for the hpd:
        frac <- (1-(hpd/100))/2
        prob <- c(frac,1-frac)

        ## Create a list where each element is a line of the grid plot:
        ## One list for each matrix to be plotted.
        LL <- lapply(p, function(x) lapply(1:dd, function(y) data.frame( t( sapply(chain$matrix[[x]], function(x) x[y,] ) ) ) ) )
        
        ## Calculate the limits for the HPD for each cell in the matrix:
        ## Now this is a list of lists and have the info for all the matrices.
        qq.list <- list()
        for( w in 1:length(LL) ){
            qq.list[[w]] <- list()
            qq.count <- 1
            for(i in 1:dd){
                for(j in i:dd){
                    qq.list[[w]][[qq.count]] <- quantile(x=LL[[w]][[i]][,j], probs=prob)
                    qq.count <- qq.count + 1
                }
            }
        }

        ## Make Sample Of Matrices To The ellipse plots:
        ## Here I am going to make sure that each sample is within the HPD using rejection sampling.        
        if( ll < n.lines ) stop(" n.lines is larger than number of matrices in the chain.")
        if( n.lines < 1 ) stop(" n.lines need to be > 1.")
        ss.list <- list() ## list for the samples, each element for each regime.
        sampled.LL <- list() ## The sampled matrices.
        for( i in 1:length(p) ){ ## For each regime.
            count <- 1
            ss.list[[i]] <- vector()
            while( count < n.lines ){
                ss <- sample(1:ll, size = 1) ## Take a sample
                if( ss %in% ss.list[[i]] ) next ## Check if sample is already in the vector.
                if( check.hpd(chain$matrix[[ p[i] ]][[ss]], qq.list[[i]], dd) ) next ## Check if ss is within HPD.
                ss.list[[i]][count] <- ss
                count <- count + 1
            }
            sampled.LL[[i]] <- lapply(1:dd, function(y) data.frame( t( sapply(chain$matrix[[ p[i] ]][ ss.list[[i]] ], function(x) x[y,] ) ) ) )
        }
        
        ## Get the data for the ellipse plots. The data is going to be in the same order of
        ##         the graphs in the par grid plot. Only the 'lower tri' part of the grid.
        ## Here i is the x axis and j is the y axis.
        ## In this implementation all samples are inside the HPD, so there are no ellipses to be printed in gray color.        

        ell.data <- list() ## Ellipses. This is a list with all of the regimes.
        for( w in 1:length(p) ){
            ell.data.count <- 1 ## A counter.
            ell.data[[w]] <- list()
            for(i in 2:dd){
                for(j in 1:(i-1)){                    
                    ell.data[[w]][[ell.data.count]] <- get.vcv.ellipse.all.matrix( mat=chain$matrix[[ p[w] ]][ ss.list[[w]] ], traits=c(i,j) )
                    ## ell.data[[w]][[ell.data.count]] <- get.vcv.ellipse.all.matrix( mat=sampled.LL[[w]], traits=c(i,j) )
                    ell.data.count <- ell.data.count + 1 ## Update the counter.
                }
            }
        }

        ## Create ellipse for the MLE estimate.
        ## The mle argument now is a list of matrices. The color of each line will match the color
        ##     respective regime in the same order. One more parameter will be added to control the width
        ##     of the line that is plotted as the mle. The name of the parameter will be changed to
        ##     drop the assumption that the line need to be the mle estimate.
        if( !is.null(point.matrix) ){ ## Check if the value of the argument is what we want:
            if( !class( point.matrix ) == "list" ) stop( " point.matrix need to be a list of matrices." )
            if( !length( point.matrix ) == length(p) ) stop( "Lenght of point.matrix need to be equal to the number of regimes fitted to the tree." )
            if( !class( point.matrix[[1]] ) == "matrix" ) stop(" point.matrix need to be a list of matrices." )
            ell.point <- list()
            for( w in 1:length(p) ){
                ell.point[[w]] <- list()
                ell.point.count <- 1
                for(i in 2:dd){
                    for(j in 1:(i-1)){
                        ell.point[[w]][[ell.point.count]] <- get.vcv.ellipse.matrix( mat=point.matrix[[w]], traits=c(i,j) )
                        ell.point.count <- ell.point.count + 1
                    }
                }
            }
        }
            
        ## Get the data for the histogram plots. Also in the same order as the plotting.
        ## Here using the data to calculate the range for the x and y axes.
        ## Here i is the x axis and j is the y axis.
        ## Note that we need to check for each regime whether the x and y limits need to be amplified or not.
        ## We want to fit all the matrices at once.
        y.hist <- vector()
        x.hist <- vector()
        ## First need to find the xlim of the plots. Will need this quantity to set the breaks for
        ##      the histograms.
        for(w in 1:length(p) ){
            for(i in 1:dd){
                for(j in i:dd){
                    x.hist <- c( min(x.hist[1], LL[[w]][[i]][,j], na.rm=TRUE) , max(x.hist[2], LL[[w]][[i]][,j], na.rm=TRUE) )
                }
            }
        }
        
        if( is.null(set.xlim) ){ ## Check for a user defined limits.
            ## Calculate the window of each bar:
            wd <- (x.hist[2] - mean(x.hist))/20
            ## Make a sequence for the breaks. Note that we add one window (bar) to the borders.
            brk <- seq(from=x.hist[1]-wd, to=x.hist[2]+wd, by=wd)
            xlim.hist <- x.hist
        } else{
            wd <- (set.xlim[2] - mean(set.xlim))/20
            brk <- seq(from=x.hist[1]-wd, to=x.hist[2]+wd, by=wd)
            xlim.hist <- set.xlim
        }

        ## Note that we have here one list for each matrix. Now the number of matrices is a variable.
        hists <- list() ## The histograms. Using the plot function and setting plot=FALSE
        ccat <- list() ## The categories, the breaks to be used in the histograms.
        for(w in 1:length(p) ){
            hists[[w]] <- list()
            ccat[[w]] <- list()
            hist.count <- 1
            for(i in 1:dd[1]){
                for(j in i:dd[1]){
                    ## Remember that the list 'LL' were made over the lines of the grid.
                    ## This is to calculate the x and y limits of the histogram.
                    hists[[w]][[hist.count]] <- hist(LL[[w]][[i]][,j], plot=FALSE, breaks=brk)
                    ## Create the cuts for the hpd:
                    ccat[[w]][[hist.count]] <- cut(hists[[w]][[hist.count]]$breaks
                                                 , c(-Inf, qq.list[[w]][[hist.count]][1], qq.list[[w]][[hist.count]][2], Inf))             
                    y.hist <- max(y.hist, hists[[w]][[hist.count]]$density)
                    hist.count <- hist.count + 1
                }
            }
        }
        ylim.hist <- c(0, y.hist)

        ## Get the xlim and ylim for ellipse plots:
        ## This is just like the original case. But here the function will take into account the limits of both
        ##       mat1 and mat2 in the case that the user have two matrices.
        ## Important to note that the new implementation will not plot ellipses out of the HPD, also the width
        ##       of the lines and the transparency can be modified.
        ell.lim <- lapply( do.call(c, ell.data), function(x) x[[1]] )
        ell.lim <- do.call(rbind, ell.lim)
        ell.lim <- apply(ell.lim, 2, range)

        ## Save the original par:
        old.par <- par(no.readonly = TRUE)

        ## Set par block:
        par(mfrow = c(dd, dd))
        par(cex = 0.6)
        par(mar = c(0, 0, 0, 0), oma = c(2, 2, 2, 2))
        par(tcl = -0.25)
        par(mgp = c(2, 0.6, 0))
        ## par(xpd = NA) ## To make plots outside plotting area.

        ## Plot the graphs in the grid:
        ## The if and else make sure that the upper.tri and the lower.tri get the correct plots.
        ## Both the ellipse lines and the densities for the diagonal and off-diagonal elements
        ##     are in lists objects with the same order of the plots in the grid. Then at each
        ##     plot the function will call the next element of the correspondent list.
        ## The ell and dens.plot.count are here to help to call the correct data to put in each plot.

        ## Produce some of the colors.
        ## Need to create the off-diagonal colors based on the alpha for the colors in the diagonal.
        ## Repeat the process for all elements of the graph, and now it will need to be a vector.
        ## colors has the color for each regime.
        colOff <- adjustcolor(col=colors, alpha.f=alphaOff)
        colDiag <- adjustcolor(col=colors, alpha.f=alphaDiag)
        colEll <- adjustcolor(col=colors, alpha.f=alphaEll)
        if( is.null(point.color) ){ ## Set the colors for the point matrices equal to the colors if not provided.
            point.color <- colors
        }

        ell.plot.count <- 1
        hist.plot.count <- 1
        for(i in 1:dd){
            for(j in 1:dd){
                if(j >= i){
                    plot(1, xlim=xlim.hist, ylim=ylim.hist, axes=FALSE, type="n", xlab="", ylab="")
                    mid <- mean(xlim.hist)
                    first.quart <- xlim.hist[1] + (mid - xlim.hist[1])/2
                    second.quart <- mid + (xlim.hist[2] - mid)/2
                    lines(x=c(mid,mid), y=ylim.hist, type="l", lty = 3, col="grey")
                    lines(x=c(first.quart, first.quart), y=ylim.hist, type="l", lty = 3, col="grey")
                    lines(x=c(second.quart, second.quart), y=ylim.hist, type="l", lty = 3, col="grey")
                    if(show.zero == TRUE){
                        lines(x=c(0,0), y=ylim.hist, type="l", lty = 3, col="blue")
                    }
                    box(col="grey")
                    if(j != i){
                        for( w in 1:length(p) ){
                            plot(hists[[w]][[hist.plot.count]], add=TRUE, freq=FALSE, border="gray"
                               , col=c("white", colOff[w], "white")[ccat[[w]][[hist.plot.count]]] )
                        }
                    } else{
                        for( w in 1:length(p) ){
                            plot(hists[[w]][[hist.plot.count]], add=TRUE, freq=FALSE, border="black"
                               , col=c("white", colDiag[w], "white")[ccat[[w]][[hist.plot.count]]] )
                        }
                        if(i == dd[1]){
                            axis(1, at=round(c(xlim.hist[1], mean(xlim.hist), xlim.hist[2]), digits = 2) )
                        }
                    }
                    if(i == 1){
                        mtext(text=leg[j], side=3, cex=l.cex)
                    }
                    if(j == 1){
                        mtext(text=leg[i], side=2, cex=l.cex)
                    }
                    if( !is.null(point.matrix) ){
                        for( w in 1:length(p) ){
                            lines(x=c(point.matrix[[w]][i,j], point.matrix[[w]][i,j]), y=ylim.hist, type="l", col=point.color[w], lwd=point.wd)
                        }
                    }
                    hist.plot.count <- hist.plot.count + 1
                } else{
                    plot(1, xlim=ell.lim[,1], ylim=ell.lim[,2], axes=FALSE, type="n", xlab="", ylab="")
                    box(col="grey")
                    for( w in 1:length(p) ){
                        invisible( lapply(ell.data[[w]][[ell.plot.count]][[2]], points, col = colEll[w]
                                        , type = "l", lwd = ell.wd) )
                    }
                    if( !is.null(point.matrix) ){
                        for( w in 1:length(p) ){
                            invisible( points(ell.point[[w]][[ell.plot.count]], col=point.color[w], type="l", lwd=point.wd) )
                        }
                    }
                    ell.plot.count <- ell.plot.count + 1
                    if(j == 1){
                        mtext(text=leg[i], side=2, cex=l.cex)
                    }
                }            
            }
        }
        
        ## Return the old parameters:
        par(old.par)
    } else{

        ########################################################################################
        ## BLOCK FOR POSTERIOR HPD IS NOT HIGHLIGHTED.
        ## This is used when hpd == 100.
        
        ## Create a list where each element is a line of the grid plot:
        LL <- lapply(p, function(x) lapply(1:dd, function(y) data.frame( t( sapply(chain$matrix[[x]], function(x) x[y,] ) ) ) ) )

        ## Get the data for the ellipse plots. The data is going to be in the same order of
        ##         the graphs in the par grid plot. Only the 'lower tri' part of the grid.
        ## Here i is the x axis and j is the y axis.
        ell.data <- list()
        for( w in 1:length(p) ){
            ell.data[[w]] <- list()
            ell.data.count <- 1    
            for(i in 2:dd){
                for(j in 1:(i-1)){
                    ell.data[[w]][[ell.data.count]] <- get.vcv.ellipse.matrix(mat=chain$matrix[[ p[w] ]], traits=c(i,j) )
                    ell.data.count <- ell.data.count + 1
                }
            }
        }
        
        ## Create ellipse for the MLE estimate.
        if( !is.null(point.matrix) ){ ## Check if the value of the argument is what we want:
            if( !class( point.matrix ) == "list" ) stop( " point.matrix need to be a list of matrices." )
            if( !length( point.matrix ) == length(p) ) stop( "Lenght of point.matrix need to be equal to the number of regimes fitted to the tree." )
            if( !class( point.matrix[[1]] ) == "matrix" ) stop(" point.matrix need to be a list of matrices." )
            ell.point <- list()
            for( w in 1:length(p) ){
                ell.point[[w]] <- list()
                ell.point.count <- 1
                for(i in 2:dd){
                    for(j in 1:(i-1)){
                        ell.point[[w]][[ell.point.count]] <- get.vcv.ellipse.matrix( mat=point.matrix[[w]], traits=c(i,j) )
                        ell.point.count <- ell.point.count + 1
                    }
                }
            }
        }

        ## Get the data for the histogram plots. Also in the same order as the plotting.
        ## Here using the data to calculate the range for the x and y axes.
        y.hist <- vector()
        x.hist <- vector()
        ## First need to find the xlim of the plots. Will need this quantity to set the breaks for
        ##      the histograms.
        for(w in 1:length(p) ){
            for(i in 1:dd){
                for(j in i:dd){
                    x.hist <- c( min(x.hist[1], LL[[w]][[i]][,j], na.rm=TRUE) , max(x.hist[2], LL[[w]][[i]][,j], na.rm=TRUE) )
                }
            }
        }
        
        if(is.null(set.xlim)){
            ## Calculate the window of each bar:
            wd <- (x.hist[2] - mean(x.hist))/20
            ## Make a sequence for the breaks. Note that we add one window to the borders.
            brk <- seq(from=x.hist[1]-wd, to=x.hist[2]+wd, by=wd)
            xlim.hist <- x.hist
        } else{
            wd <- (set.xlim[2] - mean(set.xlim))/20
            brk <- seq(from=x.hist[1]-wd, to=x.hist[2]+wd, by=wd)
            xlim.hist <- set.xlim
        }

        hists <- list() ## Here there is no ccat because there is no "white region" in the histograms.
        for( w in 1:length(p) ){
            hists[[w]] <- list()
            hist.count <- 1
            for(i in 1:dd){
                for(j in i:dd){
                    ## Remember that the list 'LL' where made over the lines of the grid.
                    ## This is to calculate the x and y limits of the histogram.
                    hists[[w]][[hist.count]] <- hist(LL[[w]][[i]][,j], plot=FALSE, breaks=brk)
                    y.hist <- max(y.hist, hists[[w]][[hist.count]]$density)
                    hist.count <- hist.count + 1
                }
            }
        }
        ylim.hist <- c(0,y.hist)

        ## Get the xlim and ylim for ellipse plots:
        ell.lim <- lapply( do.call(c, ell.data), function(x) x[[1]] )
        ell.lim <- do.call(rbind, ell.lim)
        ell.lim <- apply(ell.lim, 2, range)
        
        ## Save the original par:
        old.par <- par(no.readonly = TRUE)

        ## Set par block:
        par(mfrow = c(dd[1], dd[1]))
        par(cex = 0.6)
        par(mar = c(0, 0, 0, 0), oma = c(2, 2, 2, 2))
        par(tcl = -0.25)
        par(mgp = c(2, 0.6, 0))
        ## par(xpd = NA) ## To make plots outside plotting area.

        ## Create the color for the plot:
        colOff <- adjustcolor(col=colors, alpha.f=alphaOff)
        colDiag <- adjustcolor(col=colors, alpha.f=alphaDiag)
        colEll <- adjustcolor(col=colors, alpha.f=alphaEll)
        if( is.null(point.color) ){ ## Set the colors for the point matrices equal to the colors if not provided.
            point.color <- colors
        }

        ## Plot the graphs in the grid:
        ## The if and else make sure that the upper.tri and the lower.tri get the correct plots.
        ## Both the ellipse lines and the densities for the diagonal and off-diagonal elements
        ##     are in lists objects with the same order of the plots in the grid. Then at each
        ##     plot the function will call the next element of the correspondent list.
        ## The ell and dens.plot.count are here to help to call the correct data to put in each plot.
        
        ell.plot.count <- 1
        hist.plot.count <- 1
        for(i in 1:dd){
            for(j in 1:dd){
                if(j >= i){
                    plot(1, xlim=xlim.hist, ylim=ylim.hist, axes=FALSE, type="n", xlab="", ylab="")
                    mid <- mean(xlim.hist)
                    first.quart <- xlim.hist[1] + (mid - xlim.hist[1])/2
                    second.quart <- mid + (xlim.hist[2] - mid)/2
                    lines(x=c(mid,mid), y=ylim.hist, type="l", lty = 3, col="grey")
                    lines(x=c(first.quart, first.quart), y=ylim.hist, type="l", lty = 3, col="grey")
                    lines(x=c(second.quart, second.quart), y=ylim.hist, type="l", lty = 3, col="grey")
                    if(show.zero == TRUE){
                        lines(x=c(0,0), y=ylim.hist, type="l", lty = 3, col="blue")
                    }
                    box(col="grey")
                    if(j != i){
                        for( w in 1:length(p) ){
                            plot(hists[[w]][[hist.plot.count]], add=TRUE, freq=FALSE, col=colOff[w], border="gray")
                        }
                    } else{
                        for( w in 1:length(p) ){
                            plot(hists[[w]][[hist.plot.count]], add=TRUE, freq=FALSE, col=colDiag[w], border="black")
                        }
                        if(i == dd){ axis(1, at=round(c(xlim.hist[1],mean(xlim.hist),xlim.hist[2]), digits = 2) ) }
                    }
                    if(i == 1){
                        mtext(text=leg[j], side=3, cex=l.cex)
                    }
                    if(j == 1){
                        mtext(text=leg[i], side=2, cex=l.cex)
                    }
                    if( !is.null(point.matrix) ){
                        for( w in 1:length(p) ){
                            lines(x=c(point.matrix[[w]][i,j], point.matrix[[w]][i,j]), y=ylim.hist, type="l", col=point.color[w], lwd=point.wd)
                        }
                    }
                    hist.plot.count <- hist.plot.count + 1
                } else{
                    plot(1, xlim=ell.lim[,1], ylim=ell.lim[,2], axes=FALSE, type="n", xlab="", ylab="")
                    box(col="grey")
                    for( w in 1:length(p) ){
                    invisible( lapply(ell.data[[w]][[ell.plot.count]][[2]], points, col = colEll[w]
                                    , type = "l", lwd = ell.wd) )
                    }
                    if( !is.null(point.matrix) ){
                        for( w in 1:length(p) ){
                            invisible( points(ell.point[[w]][[ell.plot.count]], col=point.color[w], type="l", lwd=point.wd) )
                        }
                    }
                    ell.plot.count <- ell.plot.count+1
                    if(j == 1){
                        mtext(text=leg[i], side=2, cex=l.cex)
                    }
                }            
            }
        }
        ## Return the old parameters:
        par(old.par)
    }
}
