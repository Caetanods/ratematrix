##' Generates a plate with plots showing the posterior distribution of evolutionary rate matrices.
##' 
##' The function provides the option to plot a single evolutionary rate matrix on top of the posterior distribution of each regime as a vertical line on the upper-diagonal and diagonal histogram plots and as an ellipse on the lower-diagonal plots. This can be set using the argument 'point.matrix' (as well as the 'point.color' and 'point.wd' options). One can use this option to contrast the posterior with some point estimate or summary statistics. \cr
##' \cr
##' Colors can be provided either as color names recognized by R-base plot functions or in the HEX format. \cr
##' \cr
##' The lines showed by the ellipse plots (lower-diagonal) are a sample from the posterior distribution. The user can set the number of lines plotted using the argument 'n.lines'. Note that more lines will take more time to plot. \cr
##' \cr
##' The 'hpd' argument can be used to set some regions of the plot to be colored in white. For example, if 'hpd=95' the histograms will plot the region outside the 95\% HPD (Highest Posterior Density) in white and ellipse lines will only be showed if within this 95\% HPD interval. If the region chosen is too small (~10\% or lower), the plot might return an error. This happens because the function take random samples from the posterior distribution to plot as ellipse lines and exclude the samples that are outside the defined HPD interval. If this happens, try to choose a more inclusive percentage or increase the number of samples taken for the ellise lines (see argument 'n.lines') or repeat the plot until sucessful. [A better solution for this issue will be provided soon.] The default is 100 (no highlight is performed and ellipse lines are not restricted). \cr
##' \cr
##' The plots are divided into three groups. Upper-diagonal plots show histograms with the posterior distribution for the covariance values between each pairwise combination of the traits. Plot in the diagonal show histograms with the posterior distribution of evolutionary rates for each trait. Plots on the lower-diagonal slots show a collection of ellipses sampled from the posterior distribution of the model. Each ellipse line represents a bivariate distribution for the 95% CI given the pairwise correlation between the traits. \cr
##' \cr
##' Lower-diagonal plots are ideal to visualize the evolutionary correlation and variance between two traits. The orientation of the ellipses show whether there is a positive, negative or lack of correlation (horizontal or vertical orientation) between traits. The shape of the ellipses show the major axis of variation between traits. A 'cigar-shaped' ellipse indicates that one of the traits show faster evolutionary rates than the other, so one axis of variation is much larger than the other whereas a more circular (round) ellipse is a result of comparable rates of evolution between the two traits. A completely circular shape denotes lack of evolutionary correlation between two traits. It might help to understand the meaning of the ellipses lines by imagining each ellipse line marks the spread of the dots in a scatterplot with data generated with a particular covariance value (i.e., the covariance value the ellipse is representing).
##' @title Plot the distribution of evolutionary rate matrices
##' @param chain the posterior distribution of parameter estimates as produced by 'readMCMC' or samples from the prior using 'samplePrior'..
##' @param p a numeric vector with the regimes to be plotted. This parameter can be used to subset the rate regimes to be plotted as well as to control the order of the plotting. If 'NULL' (default), then all rate regimes are plotted in the same order as in the data.
##' @param colors a vector with colors for each rate regime with length equal to the number of regimes or to the number of regimes provided to the argument 'p'. If not provided the function will use pre-selected colors up to 8 regimes.
##' @param set.xlim user limits for the x axes. Need to be a vector with two elements, the minimum and the maximum.
##' @param set.leg user defined legends for the trait names. A character vector with same length as the number of traits in the model.
##' @param l.cex the 'cex' parameter for legends of the plot. See 'help(par)' for more information on 'cex'. Default is 0.7 .
##' @param ell.wd a number for the width of the ellipse lines. Default is 0.5 .
##' @param alphaOff a number between 0 and 1 with the transparency of the color used for the off-diagonal plots. Default is 1.
##' @param alphaDiag a number between 0 and 1 with the transparency of the color used for the diagonal plots. Default is 1.
##' @param alphaEll a number between 0 and 1 with the transparency of the color used for the lines of the ellipse plots. Using transparency in the lines might enhance the visualization of regions with more or less density of samples. Default is 1.
##' @param hpd a number between 0 and 100 to set the proportion of the highest posterior density (HPD) to be highlighted in the plot.
##' @param show.zero whether to plot a thin blue line showing the position of the 0 value on the histograms.
##' @param n.lines number of lines to be displayed in the ellipse plots. Default is 50 lines.
##' @param n.points number of points used to approximate the shape of the ellipses. 
##' @param point.matrix optional argument. A list of variance-covariance matrices with length equal to p. If p=NULL then length need to be equal to the number of rate regimes fitted to the data. Each element of the list will be plotted as a single line on top of the distribution of parameter estimates.
##' @param point.color optional argument. A vector with color names for the matrices set in 'point.matrix'. The vector need to have same length as 'point.matrix. If not provided, the colors of the lines will be equal to the colors of the distribution (argument 'colors').
##' @param point.wd optional argument. The width of the lines plotted using 'point.matrix'. Default is 0.5 .
##' @return A plate with a grid of plots with dimension equal to the number of traits fitted to the data.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' \donttest{
##' data( centrarchidae )
##' dt.range <- t( apply( centrarchidae$data, 2, range ) )
##' ## The step size for the root value can be set given the range we need to sample from:
##' w_mu <- ( dt.range[,2] - dt.range[,1] ) / 10
##' par.sd <- cbind(c(0,0), sqrt( c(1,1) ))
##' prior <- makePrior(r=2, p=2, par.mu=dt.range, par.sd=par.sd)
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, prior=prior
##'                          , gen=50000, w_mu=w_mu, dir=tempdir())
##' posterior <- readMCMC(handle, burn = 0.2, thin = 10)
##' plotRatematrix( posterior )
##' plotRatematrix( posterior, colors = c("black","red"))
##' plotRatematrix( posterior, colors = c("black","red"), alphaOff = 0.5)
##' plotRatematrix( posterior, colors = c("black","red"), alphaOff = 0.5, alphaDiag = 0.5)
##' plotRatematrix( posterior, colors = c("black","red"), alphaOff = 0.5, alphaDiag = 0.5)
##' ref.matrix <- list( rbind(c(0.5,0),c(0,0.5)), rbind(c(0.5,0),c(0,0.5)) )
##' plotRatematrix( posterior, colors = c("black","red"), alphaOff = 0.5, alphaDiag = 0.5
##'                , point.matrix = ref.matrix)
##' plotRatematrix( posterior, colors = c("black","red"), alphaOff = 0.5, alphaDiag = 0.5
##'                , point.matrix = ref.matrix, point.color = "orange", point.wd = 3)
##' plotRatematrix( posterior, colors = c("black","red"), alphaOff = 0.5, alphaDiag = 0.5
##'                , point.matrix = ref.matrix, point.color = "orange", point.wd = 3
##'                , alphaEll = 0.05, n.lines = 2000)
##' plotRatematrix( posterior, colors = c("black","red"), alphaOff = 0.5, alphaDiag = 0.5
##'                , point.matrix = ref.matrix, point.color = "orange", point.wd = 3
##'                , alphaEll = 0.5, n.lines = 200, hpd = 90)
##' }
plotRatematrix <- function(chain, p=NULL, colors=NULL, set.xlim=NULL, set.leg=NULL, l.cex=0.7, ell.wd=0.5, alphaOff=1, alphaDiag=1, alphaEll=1, hpd=100, show.zero=FALSE, n.lines=50, n.points=200, point.matrix=NULL, point.color=NULL, point.wd=0.5){

    ## Quickly check the correct format for chain.
    ll_class <- length( class( chain ) ) ## Usually a single class.
    if( ll_class == 1 ){
        correct_class <- grepl(pattern = "ratematrix", x = class( chain ) )
    } else{
        ## In this case one class might be "list" and the other a "ratematrix" type.
        correct_class <- any( sapply(class( chain ), function(x) grepl(pattern = "ratematrix", x = x ) ) )
    }
    if( ! correct_class ) stop("chain argument need to be a single MCMC chain. See 'mergePosterior' function to merge multiple MCMC chains together.")
   
    ## Check if the sample is larger than 'n.lines'. Otherwise, it needs to be equal to.
    ## Here using a new element added to the mcmc object type.
    ## But also protecting in case of legacy usage.
    if( is.null(chain$n_post_samples) ){
        ## Legacy mode:
        if( n.lines > nrow( chain[[1]] ) ){
            n.lines <- nrow( chain[[1]] )
        }        
    } else{
        ## New mode:
        if( n.lines > chain$n_post_samples ){
            n.lines <- chain$n_post_samples
        }
    }
    
    ## Set default values for p and colors if none is provided:
    if( is.null(p) ){
        ## p need to be all the regimes in the data.
        if( is.list(chain$matrix) & is.matrix(chain$matrix[[1]][[1]]) ){ ## More than one regime.
            np <- length(chain$matrix)
            p <- 1:np ## All the regimes in the same order.
        } else{ ## A single regime.
            p <- 1
        }
    }
    if( is.null(colors) ){
        np <- length(p)
        if( np > 9 ) stop("Unable to generate colors for more than 9 regimes. Please define color vector using 'colors' argument.") 
        if( np == 1 ){
            colors <- "black"
        } else{
            check <- c(np < 4, 3 < np && np < 6, np > 5)
            cols <- list( c("#002244", "#69BE28", "#A5ACAF"), c("#7fc97f", "#beaed4", "#fdc086", "#386cb0", "#ffff99"), c("#bc80bd", "#d9d9d9", "#fccde5", "#b3de69", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd") )
            colors <- unlist(cols[check])[1:np]
        }
    }

    ## If custom legend is not provided, then use the names of the traits.
    ## Here using a new method, but protecting for legacy usage:
    if( is.null(set.leg) ){
        if( is.null(chain$trait_names) ){
            ## Legacy usage.
            set.leg <- colnames( chain$root )
        } else{
            ## New method.
            set.leg <- chain$trait_names
        }
    }
            
    ## Check if there is only one regime to be plotted:
    if(length(p) == 1){ ## Nothing to check here, because of single regime.        
        if( is.list(chain$matrix) & is.matrix(chain$matrix[[p]][[1]]) ){ ## The first level is a list.
            ## No fix is needed in this case.
            cat("Plotting a single regime.","\n")
            dd <- ncol( chain$matrix[[p]][[1]] )
            ll <- length( chain$matrix[[p]] )
            ## Check if the chain is a sample from the prior and correct.
            if( class( chain ) == "ratematrix_prior_sample" ){
                ## This step will rebuild the R matrix from the samples taken by the prior.
                corr <- lapply(chain$matrix[[p]], decompose.cov ) ## Over the matrices.
                rb.matrix <- lapply(1:ll, function(x) rebuild.cov(r=corr[[x]]$r, v=chain$sd[[p]][x,]^2) )
                chain$matrix[[p]] <- rb.matrix
            }
        }
        if( is.list(chain$matrix) & is.matrix(chain$matrix[[1]]) ){ ## The first level is a matrix and not a list.
            if( !p == 1 ) stop("There is only one regime in the chain, then p need to be equal to 1.")
            cat("Plotting a single regime.","\n")
            dd <- ncol( chain$matrix[[1]] )
            ll <- length( chain$matrix )
            if( class( chain ) == "ratematrix_prior_sample" ){
                ## This step will rebuild the R matrix from the samples taken by the prior.
                corr <- lapply(chain$matrix, decompose.cov ) ## Over the matrices.
                rb.matrix <- lapply(1:ll, function(x) rebuild.cov(r=corr[[x]]$r, v=chain$sd[x,]^2) )
                chain$matrix <- rb.matrix
            }
            ## Need to fix the format for the rest of the function.
            temp <- chain$matrix
            rm( chain )
            chain <- list()
            chain$matrix <- list()
            chain$matrix[[1]] <- temp
        }
    }
    
    if(length(p) > 1){
        cat("Plotting multiple regimes.","\n")
        
        ## Print a text with the association between the colors and the regimes.
        name.table <- rbind(names(chain$matrix)[p], colors)
        cat("Table with regimes and colors (names or HEX):\n")
        utils::write.table(format(name.table, justify="right"), row.names=F, col.names=F, quote=F)
        
        ## First do a batch of tests:
        check.mat <- vector()
        check.length <- vector()
        for(i in 1:length(p)){
            check.mat[i] <- ncol( chain$matrix[[ p[i] ]][[1]] ) ## First element of each regime.
            check.length[i] <- length( chain$matrix[[ p[i] ]] ) ## The length of each chain regime.
        }
        equal.size <- sapply(2:length(p), function(x) check.mat[1] == check.mat[x] )
        equal.length <- sapply(2:length(p), function(x) check.length[1] == check.length[x] )

        if( !sum(equal.size) == (length(p)-1) ) stop("Matrix regimes do not have the same number of traits.")
        if( !sum(equal.length) == (length(p)-1) ) stop("Chain for regimes do not have the same length.")
        
        ll <- check.length[1]
        dd <- check.mat[1]

        if( class( chain ) == "ratematrix_prior_sample" ){
            ## This step will rebuild the R matrix from the samples taken by the prior.
            for(i in p){
                corr <- lapply(chain$matrix[[i]], decompose.cov ) ## Over the matrices.
                rb.matrix <- lapply(1:ll, function(x) rebuild.cov(r=corr[[x]]$r, v=chain$sd[[i]][x,]^2) )
                chain$matrix[[i]] <- rb.matrix
            }
        }
    }

    ## ######################################################################################
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
                    qq.list[[w]][[qq.count]] <- stats::quantile(x=LL[[w]][[i]][,j], probs=prob)
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
                if( checkHPD(chain$matrix[[ p[i] ]][[ss]], qq.list[[i]], dd) ) next ## Check if ss is within HPD.
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
            for(j in 2:dd){
                for(i in 1:(j-1)){                    
                    ell.data[[w]][[ell.data.count]] <- getEllipseMatrix( mat=chain$matrix[[ p[w] ]][ ss.list[[w]] ]
                                                                        , traits=c(i,j), n.points=n.points)
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
                for(j in 2:dd){
                    for(i in 1:(j-1)){
                        ell.point[[w]][[ell.point.count]] <- getEllipseMatrix( mat=point.matrix[[w]], traits=c(i,j)
                                                                            , sample.line=n.lines, n.points=n.points)
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
            brk.var <- seq(0, to=x.hist[2]+wd, by=wd)
            xlim.hist <- x.hist
        } else{
            wd <- (set.xlim[2] - mean(set.xlim))/20
            brk <- seq(from=x.hist[1]-wd, to=x.hist[2]+wd, by=wd)
            brk.var <- seq(from=0, to=x.hist[2]+wd, by=wd)
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
                    if( i == j ){
                        ## Remember that the list 'LL' were made over the lines of the grid.
                        ## This is to calculate the x and y limits of the histogram.
                        hists[[w]][[hist.count]] <- graphics::hist(LL[[w]][[i]][,j], plot=FALSE, breaks=brk.var)
                        ## Create the cuts for the hpd:
                        ccat[[w]][[hist.count]] <- cut(hists[[w]][[hist.count]]$breaks
                                                     , c(-Inf, qq.list[[w]][[hist.count]][1], qq.list[[w]][[hist.count]][2], Inf))             
                        y.hist <- max(y.hist, hists[[w]][[hist.count]]$density)
                        hist.count <- hist.count + 1
                    } else{
                        hists[[w]][[hist.count]] <- graphics::hist(LL[[w]][[i]][,j], plot=FALSE, breaks=brk)
                        ccat[[w]][[hist.count]] <- cut(hists[[w]][[hist.count]]$breaks
                                                     , c(-Inf, qq.list[[w]][[hist.count]][1], qq.list[[w]][[hist.count]][2], Inf))             
                        y.hist <- max(y.hist, hists[[w]][[hist.count]]$density)
                        hist.count <- hist.count + 1
                    }
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
        ell.lim <- apply(ell.lim, 2, range) ## One column for each of the plots.
        ## Need to make such that the plot uses the same limits for the x and y axis. Otherwise it is hard to
        ##    compare the major axis of variance of the plots.
        ell.lim.iso <- c(min(ell.lim[1,]),max(ell.lim[2,]))

        ## Save the original par:
        old.par <- graphics::par(no.readonly = TRUE)

        ## Set par block:
        graphics::par(mfrow = c(dd, dd))
        graphics::par(cex = 0.6)
        graphics::par(mar = c(0, 0, 0, 0), oma = c(2, 2, 2, 2))
        graphics::par(tcl = -0.25)
        graphics::par(mgp = c(2, 0.6, 0))
        ## graphics::par(xpd = NA) ## To make plots outside plotting area.

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
        ## Need to adapt this to accept both color names and HEX.
        colOff <- grDevices::adjustcolor(col=colors, alpha.f=alphaOff)
        colDiag <- grDevices::adjustcolor(col=colors, alpha.f=alphaDiag)
        colEll <- grDevices::adjustcolor(col=colors, alpha.f=alphaEll)
        if( is.null(point.color) ){ ## Set the colors for the point matrices equal to the colors if not provided.
            point.color <- colors
        }

        ell.plot.count <- 1
        hist.plot.count <- 1
        for(i in 1:dd){
            for(j in 1:dd){
                if(j >= i){
                    graphics::plot(1, xlim=xlim.hist, ylim=ylim.hist, axes=FALSE, type="n", xlab="", ylab="")
                    mid <- mean(xlim.hist)
                    first.quart <- xlim.hist[1] + (mid - xlim.hist[1])/2
                    second.quart <- mid + (xlim.hist[2] - mid)/2
                    graphics::lines(x=c(mid,mid), y=ylim.hist, type="l", lty = 3, col="grey")
                    graphics::lines(x=c(first.quart, first.quart), y=ylim.hist, type="l", lty = 3, col="grey")
                    graphics::lines(x=c(second.quart, second.quart), y=ylim.hist, type="l", lty = 3, col="grey")
                    if(show.zero == TRUE){
                        graphics::lines(x=c(0,0), y=ylim.hist, type="l", lty = 3, col="blue")
                    }
                    graphics::box(col="grey")
                    if(j != i){
                        for( w in 1:length(p) ){
                            graphics::plot(hists[[w]][[hist.plot.count]], add=TRUE, freq=FALSE, border="gray"
                               , col=c("white", colOff[w], "white")[ccat[[w]][[hist.plot.count]]] )
                        }
                    } else{
                        for( w in 1:length(p) ){
                            graphics::plot(hists[[w]][[hist.plot.count]], add=TRUE, freq=FALSE, border="black"
                               , col=c("white", colDiag[w], "white")[ccat[[w]][[hist.plot.count]]] )
                        }
                        if(i == dd[1]){
                            graphics::axis(1, at=round(c(xlim.hist[1], mean(xlim.hist), xlim.hist[2]), digits = 2) )
                        }
                    }
                    if(i == 1){
                        graphics::mtext(text=set.leg[j], side=3, cex=l.cex)
                    }
                    if(j == 1){
                        graphics::mtext(text=set.leg[i], side=2, cex=l.cex)
                    }
                    if( !is.null(point.matrix) ){
                        for( w in 1:length(p) ){
                            graphics::lines(x=c(point.matrix[[w]][i,j], point.matrix[[w]][i,j]), y=ylim.hist, type="l", col=point.color[w], lwd=point.wd)
                        }
                    }
                    hist.plot.count <- hist.plot.count + 1
                } else{
                    graphics::plot(1, xlim=ell.lim.iso, ylim=ell.lim.iso, axes=FALSE, type="n", xlab="", ylab="")
                    graphics::box(col="grey")
                    for( w in 1:length(p) ){
                        invisible( lapply(ell.data[[w]][[ell.plot.count]][[2]], function(x) graphics::points(x, col = colEll[w]
                                        , type = "l", lwd = ell.wd)) )
                    }
                    if( !is.null(point.matrix) ){
                        for( w in 1:length(p) ){
                            invisible( graphics::points(ell.point[[w]][[ell.plot.count]], col=point.color[w], type="l", lwd=point.wd) )
                        }
                    }
                    ell.plot.count <- ell.plot.count + 1
                    if(j == 1){
                        graphics::mtext(text=set.leg[i], side=2, cex=l.cex)
                    }
                }            
            }
        }
        
        ## Return the old parameters:
        graphics::par(old.par)
       
    } else{

        ## ######################################################################################
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
            for(j in 2:dd){
                for(i in 1:(j-1)){
                    ell.data[[w]][[ell.data.count]] <- getEllipseMatrix(mat=chain$matrix[[ p[w] ]], traits=c(i,j)
                                                                      , sample.line=n.lines, n.points=n.points )
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
                for(j in 2:dd){
                    for(i in 1:(j-1)){
                        ell.point[[w]][[ell.point.count]] <- getEllipseMatrix( mat=point.matrix[[w]], traits=c(i,j)
                                                                            , sample.line=n.lines, n.points=n.points)
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
            brk.var <- seq(from=0, to=x.hist[2]+wd, by=wd)
            xlim.hist <- x.hist
        } else{
            wd <- (set.xlim[2] - mean(set.xlim))/20
            brk <- seq(from=x.hist[1]-wd, to=x.hist[2]+wd, by=wd)
            brk.var <- seq(from=0, to=x.hist[2]+wd, by=wd)
            xlim.hist <- set.xlim
        }

        hists <- list() ## Here there is no ccat because there is no "white region" in the histograms.
        for( w in 1:length(p) ){
            hists[[w]] <- list()
            hist.count <- 1
            for(i in 1:dd){
                for(j in i:dd){
                    if( i == j ){
                        ## Remember that the list 'LL' where made over the lines of the grid.
                        ## This is to calculate the x and y limits of the histogram.
                        hists[[w]][[hist.count]] <- graphics::hist(LL[[w]][[i]][,j], plot=FALSE, breaks=brk.var)
                        y.hist <- max(y.hist, hists[[w]][[hist.count]]$density)
                        hist.count <- hist.count + 1
                    } else{
                        hists[[w]][[hist.count]] <- graphics::hist(LL[[w]][[i]][,j], plot=FALSE, breaks=brk)
                        y.hist <- max(y.hist, hists[[w]][[hist.count]]$density)
                        hist.count <- hist.count + 1
                    }
                }
            }
        }
        ylim.hist <- c(0,y.hist)

        ## Get the xlim and ylim for ellipse plots:
        ell.lim <- lapply( do.call(c, ell.data), function(x) x[[1]] )
        ell.lim <- do.call(rbind, ell.lim)
        ell.lim <- apply(ell.lim, 2, range)
        ## Need to make sure that the y and x axes have the same range. This will make it possible to compare the major axes of
        ##    variances among the traits.
        ell.lim.iso <- c(min(ell.lim[1,]),max(ell.lim[2,]))

        
        
        ## Save the original par:
        old.par <- graphics::par(no.readonly = TRUE)

        ## Set par block:
        graphics::par(mfrow = c(dd[1], dd[1]))
        graphics::par(cex = 0.6)
        graphics::par(mar = c(0, 0, 0, 0), oma = c(2, 2, 2, 2))
        graphics::par(tcl = -0.25)
        graphics::par(mgp = c(2, 0.6, 0))
        ## graphics::par(xpd = NA) ## To make plots outside plotting area.

        ## Create the color for the plot:
        colOff <- grDevices::adjustcolor(col=colors, alpha.f=alphaOff)
        colDiag <- grDevices::adjustcolor(col=colors, alpha.f=alphaDiag)
        colEll <- grDevices::adjustcolor(col=colors, alpha.f=alphaEll)
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
                    graphics::plot(1, xlim=xlim.hist, ylim=ylim.hist, axes=FALSE, type="n", xlab="", ylab="")
                    mid <- mean(xlim.hist)
                    first.quart <- xlim.hist[1] + (mid - xlim.hist[1])/2
                    second.quart <- mid + (xlim.hist[2] - mid)/2
                    graphics::lines(x=c(mid,mid), y=ylim.hist, type="l", lty = 3, col="grey")
                    graphics::lines(x=c(first.quart, first.quart), y=ylim.hist, type="l", lty = 3, col="grey")
                    graphics::lines(x=c(second.quart, second.quart), y=ylim.hist, type="l", lty = 3, col="grey")
                    if(show.zero == TRUE){
                        graphics::lines(x=c(0,0), y=ylim.hist, type="l", lty = 3, col="blue")
                    }
                    graphics::box(col="grey")
                    if(j != i){
                        for( w in 1:length(p) ){
                            graphics::plot(hists[[w]][[hist.plot.count]], add=TRUE, freq=FALSE, col=colOff[w], border="gray")
                        }
                    } else{
                        for( w in 1:length(p) ){
                            graphics::plot(hists[[w]][[hist.plot.count]], add=TRUE, freq=FALSE, col=colDiag[w], border="black")
                        }
                        if(i == dd){ graphics::axis(1, at=round(c(xlim.hist[1],mean(xlim.hist),xlim.hist[2]), digits = 2) ) }
                    }
                    if(i == 1){
                        graphics::mtext(text=set.leg[j], side=3, cex=l.cex)
                    }
                    if(j == 1){
                        graphics::mtext(text=set.leg[i], side=2, cex=l.cex)
                    }
                    if( !is.null(point.matrix) ){
                        for( w in 1:length(p) ){
                            graphics::lines(x=c(point.matrix[[w]][i,j], point.matrix[[w]][i,j]), y=ylim.hist, type="l", col=point.color[w], lwd=point.wd)
                        }
                    }
                    hist.plot.count <- hist.plot.count + 1
                } else{
                    graphics::plot(1, xlim=ell.lim.iso, ylim=ell.lim.iso, axes=FALSE, type="n", xlab="", ylab="")
                    graphics::box(col="grey")
                    for( w in 1:length(p) ){
                    invisible( lapply(ell.data[[w]][[ell.plot.count]][[2]], function(x) graphics::points(x, col = colEll[w]
                                    , type = "l", lwd = ell.wd) ) )
                    }
                    if( !is.null(point.matrix) ){
                        for( w in 1:length(p) ){
                            invisible( graphics::points(ell.point[[w]][[ell.plot.count]], col=point.color[w], type="l", lwd=point.wd) )
                        }
                    }
                    ell.plot.count <- ell.plot.count+1
                    if(j == 1){
                        graphics::mtext(text=set.leg[i], side=2, cex=l.cex)
                    }
                }            
            }
        }
        ## Return the old parameters:
        graphics::par(old.par)
    }
}
