##' Plot the posterior distribution of root values sampled from the MCMC analysis or samples from the prior distribution.
##'
##' 
##' @title Plot posterior distribution of root values for the traits
##' @param chain the posterior distribution loaded from the files using 'readMCMC' or samples from the prior generated with the 'samplePrior' function.
##' @param color the color for the histograms.
##' @param set.xlab a vector with legends for the x axes. If 'NULL' (default), the names are 'trait_1' to 'trait_n". 
##' @param set.cex.lab the cex value for the labels (default is 1).
##' @param set.cex.axis the cex value for the axes numbers (default is 1.5).
##' @param set.xlim the xlim for the plot. Need to be a vector with the lower and higher bound.
##' @param hpd the Highest Posterior Density interval to highlight in the plot. Parameter values outside this interval will be colored in white. A numeric value between 0 and 100 (default is 100).
##' @param mfrow the number of rows to use in the figure (default is 1).
##' @param vline.values numeric values for plotting vertical lines. Can be a single value recycled for each of the plots or a vector with length equal to the number of traits.
##' @param vline.color character vector with colors for the vertical lines. Can be a single color if length of 'vline.values' is 1 otherwise need to have length equal to the number of traits.
##' @param vline.wd numeric value for the width of the vertical lines. Can be a single value if length of 'vline.values' is 1 otherwise need to have length equal to the number of traits.
##' @param show.zero whether a vertical line should be plotted showing the position of the value 0 in the plot.
##' @return A plot with the posterior density of root values or distribution of root values sampled from the prior.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' \donttest{
##' data( centrarchidae )
##' dt.range <- t( apply( centrarchidae$data, 2, range ) )
##' ## The step size for the root value can be set given the range we need to sample from:
##' w_mu <- ( dt.range[,2] - dt.range[,1] ) / 10
##' par.sd <- cbind(c(0,0), sqrt( c(10,10) ))
##' prior <- makePrior(r=2, p=2, den.mu="unif", par.mu=dt.range, den.sd="unif", par.sd=par.sd)
##' prior.samples <- samplePrior(n = 1000, prior = prior)
##' start.point <- samplePrior(n=1, prior=prior)
##' ## Plot the prior. Red line shows the sample from the prior that will set the starting 
##' ##      point for the MCMC.
##' plotRatematrix(prior.samples, point.matrix = start.point$matrix, point.color = "red"
##'                , point.wd = 2)
##' plotRootValue(prior.samples)
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, prior=prior
##'                          , gen=10000, w_mu=w_mu, dir=tempdir())
##' posterior <- readMCMC(handle, burn = 0.2, thin = 10)
##' ## Again, here the red line shows the starting point of the MCMC.
##' plotRatematrix( posterior, point.matrix = start.point$matrix, point.color = "red"
##'                , point.wd = 2)
##' plotRootValue(posterior)
##' }
plotRootValue <- function(chain, color="black", set.xlab=NULL, set.cex.lab=1, set.cex.axis=1.5, set.xlim=NULL, hpd=100, mfrow=1, vline.values=NULL, vline.color=NULL, vline.wd=NULL, show.zero=FALSE){
    ## Plot Root value. Will work just like the plotRatematrix.
    ntr <- ncol( chain$root )

    ## If custom legend is not provided, then use the names of the traits.
    if(is.null(set.xlab)){
       set.xlab <- colnames( chain$root )
    }

    ## Calculate the xlim from the data:
    x.hist <- c( min( apply(chain$root, 2, min) ), max( apply(chain$root, 2, max) ) )

    if( is.null(set.xlim) ){ ## Check for a user defined limits.
        ## Calculate the window of each bar:
        wd <- (x.hist[2] - mean(x.hist))/20
        ## Make a sequence for the breaks. Note that we add one window (bar) to the borders.
        brk <- seq(from=x.hist[1]-wd, to=x.hist[2]+wd, by=wd)
    } else{
        wd <- (set.xlim[2] - mean(set.xlim))/20
        brk <- seq(from=x.hist[1]-wd, to=x.hist[2]+wd, by=wd)
        x.hist <- set.xlim
    }

    ## Calculate the disposition of the plots.
    graph.col <- ceiling(ntr/mfrow)
    graph.dim <- c(mfrow, graph.col)
    
    ## Calculating region to highlight with the 'hpd':
    if(hpd < 100){
        ## Calculate the probabilities for the hpd:
        frac <- (1-(hpd/100))/2
        prob <- c(frac, 1-frac)
        qq.mat <- matrix(data=NA, nrow=ntr, ncol=2)
        for( i in 1:ntr ){
            qq.mat[,i] <- stats::quantile(x=chain$root[,i], probs=prob)
        }
    }

    y.hist <- vector()
    hists <- list() ## The histograms. Using the plot function and setting plot=FALSE
    ccat <- list() ## The categories, the breaks to be used in the histograms.
    for( i in 1:ntr ){
        ## This is to calculate the x and y limits of the histogram.
        hists[[i]] <- graphics::hist(chain$root[,i], plot=FALSE, breaks=brk)             
        y.hist <- max(y.hist, hists[[i]]$density, na.rm=TRUE)
    }
    ylim.hist <- c(0,y.hist) ## Limits for the y axis start from 0.

    ## Create the cuts for the hpd:
    if(hpd < 100){
        for( i in 1:ntr ) ccat[[i]] <- cut(hists[[i]]$breaks, c(-Inf, qq.mat[1,i], qq.mat[2,i], Inf) )
    }

    ## Calculate things for the vertical line.
    ## vline.values=NULL, vline.color=NULL, vline.wd=0.5
    if( !is.null( vline.values ) ){
        if( length( vline.values ) > 1 ){
            if( !length( vline.values ) == length( vline.color ) ) stop("'vline.color' has to have length equal to 'vline.values'.\n")
            if( is.null( vline.wd ) ){
                vline.wd <- rep(2, times = length( vline.values ) )
            } else{
                if( !length( vline.values ) == length( vline.wd ) ) stop("'vline.wd' has to have length equal to 'vline.values'.\n")
            }
            recycle.vline <- FALSE
        } else{
            if( is.null( vline.wd ) ){
                vline.wd <- 2
            }
            recycle.vline <- TRUE
        }
    }

    ## Save the original par:
    old.par <- graphics::par(no.readonly = TRUE)

    ## Set par block:
    graphics::par(mfrow = graph.dim )
    graphics::par(cex = 0.6)
    graphics::par(oma = c(3, 3, 3, 3))

    ## Now make the plots:
    for( i in 1:ntr ){
        graphics::plot(1, xlim=x.hist, ylim=ylim.hist, axes=FALSE, type="n", xlab="", ylab="")     
        mid <- mean(x.hist)
        first.quart <- x.hist[1] + (mid - x.hist[1])/2
        second.quart <- mid + (x.hist[2] - mid)/2
        graphics::lines(x=c(mid,mid), y=ylim.hist, type="l", lty = 3, col="grey")
        graphics::lines(x=c(first.quart, first.quart), y=ylim.hist, type="l", lty = 3, col="grey")
        graphics::lines(x=c(second.quart, second.quart), y=ylim.hist, type="l", lty = 3, col="grey")
        
        if(show.zero == TRUE){
            graphics::lines(x=c(0,0), y=ylim.hist, type="l", lty = 3, col="blue")
        }
        
        if( hpd < 100 ){
            graphics::plot(hists[[i]], add=TRUE, freq=FALSE, border="gray", col=c("white", color, "white")[ ccat[[i]] ] )
        } else{
            graphics::plot(hists[[i]], add=TRUE, freq=FALSE, border="gray", col=color)
        }

        if( !is.null( vline.values ) ){
            if( recycle.vline ){
                graphics::abline(v=vline.values[1], col=vline.color[1], lwd=vline.wd[1])
            } else{
                graphics::abline(v=vline.values[i], col=vline.color[i], lwd=vline.wd[i])
            }
        }

        
        graphics::axis(1, at=round(c(x.hist[1], mean(x.hist), x.hist[2]), digits = 2), cex.axis=set.cex.axis)
        graphics::axis(2, at=round(c(ylim.hist[1], mean(ylim.hist), ylim.hist[2]), digits = 2), cex.axis=set.cex.axis)
        graphics::mtext(text=set.xlab[i], side=1, cex=set.cex.lab, line=3)
        graphics::mtext(text="Density", side=2, cex=set.cex.lab, line=3)
    }

    ## Return the original parameters.
    graphics::par(old.par)    
}
