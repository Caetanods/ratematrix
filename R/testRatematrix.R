##' Function uses summary statistics to test for differences between the posterior distribution of parameter estimates for the evolutionary rate matrix regimes.
##'
##' This functions performs a test to check whether the posterior distribution of the fitted matrices are different. It returns the proportion of overlap between regimes. When this proportion is less than 0.05 this means that the posterior distribution of the elements of the evolutionary rate matrices does not overlap more than 5\%. This test statistics is NOT a p value! This is not an estimate of the probability of deviance from a null distribution. It assumes that when the posterior distribution of two or more paramaters do not overlap, then there is a relevant difference between the parameters. \cr
##' \cr
##' The test can be performed using the median overlap of the posterior distribution across all elements of the ratematrix or by contrasting each element separatelly. Checking each element independently provides more information. Using the median overlap will result in a single value returned, but it can be insensitive to important changes in the evolutionary rate matrices between regimes. When a posterior distribution with more than two rate regimes is fitted to the data, the function performs tests for all pairwise combinations.
##' @title Test for difference between evolutionary rate matrix estimates
##' @param chain the posterior distribution of parameter estimates as output of the function 'readMCMC'.
##' @param par the attribute of the rate matrices that are checked by the test. One of 'all', 'correlation', or 'rates' (default is 'all'). Choose 'all' to compute the summary statistics for the overall matrix. Choose 'rates' to check the rates of evolution among the traits. Choose 'correlation' to compute the summary statistics for the evolutionary correlations.
##' @param diff.test whether to return the pairwise difference between the parameters computed from the joint posterior distribution. If set to FALSE (default), then the results will be based on the overlap percentage.
##' @param median.test whether to return a median of the summary statistics across all elements of the evolutionary rate matrices. The default is FALSE.
##' @param regimes a numeric vector or character vector. This is the set of regimes to be compared. If numeric, then the regimes are in the same order as in the 'chain' argument. If character, than names need to match the names of the rate regimes fitted to the phylogenetic tree.
##' @param plot whether to plot the results of the summary statistics test applied to every element of the matrix. Default is FALSE.
##' @return Return a matrix or a list with the value of the test statistics.
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
##' testRatematrix(posterior, par = "all")
##' testRatematrix(posterior, par = "correlation")
##' testRatematrix(posterior, par = "rates")
##' testRatematrix(posterior, par = "correlation", plot = TRUE)
##' }
testRatematrix <- function(chain, par=c("all","correlation","rates"), diff.test=FALSE, median.test=FALSE, regimes=NULL, plot=FALSE){

    par <- match.arg(par)

    if( inherits(chain, what=c("ratematrix_single_chain", "ratematrix_multi_chain")) ){
        if( inherits(chain, what=c("ratematrix_single_chain")) ){
            stop("Cannot perform test with a single regime. \n")
        } else{
            p <- length( chain$matrix )
            r <- ncol( chain$root )
        }
    } else{
        stop("Arguments need to be the posterior distribution of class 'ratematrix_single_chain' or 'ratematrix_multi_chain'. See 'readMCMC'. \n")
    }

    ## Make all the combinations or use the regimes provided.
    if( is.null(regimes) ){
        comb <- utils::combn(1:p, 2)
    } else{
        if( !is.numeric( regimes ) & !is.character( regimes ) ) stop("Argument 'regimes' need to be a numeric or character vector. \n")
        if( is.numeric( regimes ) ){
            comb <- utils::combn(regimes, 2)
        }
        if( is.character( regimes ) ){
            nm <- names( chain$matrix )
            if( as.logical( sum(!regimes %in% nm) ) ) stop("Names of regimes in argument 'regimes' need to match names of regimes of the posterior distribution (see 'names(chain$matrix)' ). \n")
            reg.id <- which( nm %in% regimes )
            comb <- utils::combn(reg.id, 2)
        }
    }
    
    mat.diff <- list()
    if(par == "all"){
        for(i in 1:ncol(comb)){
            mat1 <- t( sapply(chain$matrix[[comb[1,i]]], function(x) c(x) ) )
            mat2 <- t( sapply(chain$matrix[[comb[2,i]]], function(x) c(x) ) )
            mat.diff[[i]] <- mat1 - mat2
        }
    }
    if(par == "rates"){
        for(i in 1:ncol(comb)){
            mat1 <- t( sapply(chain$matrix[[comb[1,i]]], function(x) c( diag(x) ) ) )
            mat2 <- t( sapply(chain$matrix[[comb[2,i]]], function(x) c( diag(x) ) ) )
            mat.diff[[i]] <- mat1 - mat2
        }
    }
    if(par == "correlation"){
        ## This part will break if the rate matrix is 2x2. The median do not need to be calculated.
        if( ncol( chain$matrix[[1]][[1]] ) > 2 ){
            upper <- upper.tri( chain$matrix[[1]][[1]] )
            for(i in 1:ncol(comb)){
                mat1 <- t( sapply(chain$matrix[[comb[1,i]]], function(x) c( stats::cov2cor(x)[upper] ) ) )
                mat2 <- t( sapply(chain$matrix[[comb[2,i]]], function(x) c( stats::cov2cor(x)[upper] ) ) )
                mat.diff[[i]] <- mat1 - mat2
            }
        } else{
            upper <- upper.tri( chain$matrix[[1]][[1]] )
            for(i in 1:ncol(comb)){
                mat1 <- t( sapply(chain$matrix[[comb[1,i]]], function(x) c( stats::cov2cor(x)[upper] ) ) )
                mat2 <- t( sapply(chain$matrix[[comb[2,i]]], function(x) c( stats::cov2cor(x)[upper] ) ) )
                mat.diff[[i]] <- mat1 - mat2
            }
        }
    }

    if(median.test){
        median.diff <- lapply(mat.diff, function(x) apply(x, 1, stats::median))
        cdf.list <- lapply(median.diff, FUN = stats::ecdf)
        qq.list <- lapply(cdf.list, FUN = function(x) x(0) )
        test <- lapply(qq.list, FUN = function(x) 2*apply(cbind(x, 1-x), 1, min) )
        test.dt <- do.call(cbind, test)
        colnames(test.dt) <- paste("mat #",comb[1,], " x #", comb[2,], sep="")
        rownames(test.dt) <- "test value"
        return(test.dt)
    }
    
    if(!median.test){
        overlap <- lapply(mat.diff, FUN=function(x) getOverlap(x, r=r))

        ## The output for the correlations looks strange because we loose the info of the order of the traits.
        ## Here I will re-format the output if the test is made on the correlations.
        if(par == "correlation"){
            overlap.mat <- lapply(1:length(mat.diff), function(x) matrix(ncol=r, nrow=r))
            upper <- upper.tri(x=matrix(ncol=r, nrow=r))
            for( i in 1:length(mat.diff) ) overlap.mat[[i]][upper] <- overlap[[i]]
        } else{
            overlap.mat <- lapply(overlap, function(x) matrix(x, ncol=r, byrow = TRUE) )
        }
        
        if( is.null(names(chain$matrix)) ){
            main <- paste("Regime #",comb[1,], " x #", comb[2,], sep="")
        } else{
            main <- paste("Regime ", names(chain$matrix)[comb[1,]], " x ", names(chain$matrix)[comb[2,]], sep="")
        }
        if(plot){
            n <- length( overlap.mat )
            old.par <- par(no.readonly = TRUE)
            ## Decide the size of the mfrow plot to return. Works up to 12 regimes. Otherwise the function will just return a lot of plots.
            if( n == 1 ) plotHeatmat(mat=overlap.mat[[1]], r=r, par=par, main=main[1])
            if( n > 1 & n < 4 ) par( mfrow = c(1,3) )
            if( n == 4 ) par( mfrow = c(2,2) )
            if( n > 4 & n < 7 ) par( mfrow = c(2,3))
            if( n > 6 & n < 10 ) par( mfrow = c(3,3))
            if( n > 9 & n < 13 ) par( mfrow = c(3,4))
            lapply(1:n, function(x) plotHeatmat(mat=overlap.mat[[x]], r=r, par=par, main=main[x]) )
            par(old.par)
        }
        names( overlap.mat ) <- main
        return( list(overlap.mat) )
    }
}

## Some helping functions for the output:
getOverlap <- function(mt.dif, r){
    if( r == 2 ){
        cdf <- stats::ecdf(mt.dif[1,])
        qq <- cdf(0)
        test <- 2 * min(qq,1-qq)
        return( test ) ## Single number.
    } else{
        nc <- ncol(mt.dif)
        cdf.list <- lapply(1:nc, function(x) stats::ecdf(mt.dif[,x]))
        qq.list <- lapply(cdf.list, FUN = function(x) x(0) )
        test <- lapply(qq.list, FUN = function(x) 2*apply(cbind(x, 1-x), 1, min) )
        test.dt <- do.call(cbind, test)
        return( as.vector(test.dt) )
    }
}

plotHeatmat <- function(mat, r, par, main){
    if( par=="rates" ){
        na.mat <- matrix(NA, nrow=r, ncol=r)
        diag(na.mat) <- unlist( mat )
        plot.mat <- na.mat
    } else{
        plot.mat <- mat
    }
    
    rot.mat <- t(apply(plot.mat, 2, rev))
    ## Need to rotate the matrix 90 degrees because of the behavior of the 'image' function.
    x <- y <- 1:r
    old.par <- par(no.readonly = TRUE)
    par( mar=c(2,2,4,2) )
    ## Set a better color pallette.
    col <- grDevices::heat.colors(15)
    breaks <- c(0, exp(seq(log(0.0001), log(1), length.out = 15)))
    graphics::image(x=x, y=y, z=rot.mat, xaxt= "n", yaxt= "n", xlab="", ylab="", main=main, col=col, breaks=breaks)
    x.text <- rep(1:r, times=r)
    y.text <- rep(1:r, each=r)
    graphics::text(x=x.text, y=y.text, round(rot.mat, digits=4))
    par(old.par)
}
