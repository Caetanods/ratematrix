##' Read the log file produced by the 'ratematrixMCMC' function. Calculates acceptance ration and show the trace plot.
##'
##' 
##' @title Make analysis of the log file of the MCMC.
##' @param out the output object from the 'ratematrixMCMC' function.
##' @param burn the proportion of burn-in. A numeric value between 0 and 1.
##' @param thin the number of generations to skip when reading the posterior distribution.
##' @param show.plots whether to show a trace plot of the log-likelihood and the acceptance ratio.
##' @param print.result whether to print the results of the acceptance ratio to the screen.
##' @param dir the directory where to find the log file. If set to NULL the function will search for the files in the same directory that the MCMC chain was made.
##' @return A named vector with the acceptance ratio for the whole MCMC and each parameter. If a list of phylogenetic trees was provided to the MCMC chain, then the output is a list with the acceptance ration for the parameters and a table showing the frequency in which each of the phylogenies was present when a move step was accepted.
##' @export
##' @author daniel
logAnalizer <- function(out, burn=0.25, thin=1, show.plots=TRUE, print.result=TRUE, dir=NULL){
    ## Need to read the log output. Then make some plots or such.

    ## Check the directory.
    if(is.null(dir)){
        direct <- out$dir
    } else{
        direct <- dir
    }

    ## Read the log and create a matrix.
    post <- seq(round(out$gen * burn)+1, out$gen+1, by=thin) ## First line is the header.
    log.mcmc <- read_lines( file=file.path(direct, paste(out$outname, ".", out$ID, ".log", sep="")) )
    header <- log.mcmc[1]
    header <- as.character( strsplit(x=header, split=";", fixed=TRUE)[[1]] )
    log.mcmc <- log.mcmc[post]
    log.mcmc <- t( sapply(log.mcmc, function(x) as.numeric( strsplit(x=x, split=";", fixed=TRUE)[[1]] )
                        , USE.NAMES=FALSE) )
    colnames( log.mcmc ) <- header

    ## Make analyses:
    accept <- sapply(1:4, function(x) sum(log.mcmc[,x]) / dim(log.mcmc)[1] )
    names(accept) <- c("All", "correlation", "sd", "root")

    mix.phylo <- table( log.mcmc[ which(log.mcmc[,1] == 1), 5] ) / length( log.mcmc[ log.mcmc[,1],5 ] )

    post <- seq(round(2000 * 0.25)+1, 2000+1, by=1) ## First line is the header.
    at.gen <- round( seq(from=1, to=length(post), length.out = 5) )
    labels.gen <- round( seq(from=post[1], to=post[length(post)], length.out = 5) )
    cum.accept <- cumsum(log.mcmc[,1]) / 1:length(log.mcmc[,1])

    if( print.result==TRUE ){
        cat("Acceptance ratio for the MCMC and parameters:\n")
        print( accept )
        if( is.list(out$phy[[1]]) ){
            cat("\n")
            cat("Acceptance ratio for each phylogeny:\n")
            print( mix.phylo )
        }
    }

    if( show.plots==TRUE ){
        ## Log-likelihood and acceptance ration trace plot:
        old.par <- par(no.readonly = TRUE)

        par( mfrow = c(2,1) )
        par(mar = c(1, 0, 0, 0), oma = c(3, 4, 2, 0))
        plot(x=1:dim( log.mcmc )[1], y=log.mcmc[,6], type="l", axes=FALSE, xlab="", ylab="")
        axis(side=2)
        mtext("Log-likelihood", side=2, line=2, cex=1)
        plot(x=1:dim( log.mcmc )[1], y=cum.accept, type="l", axes=FALSE, xlab="", ylab="")
        axis(side=2)
        axis(side=1, at=at.gen, labels=labels.gen)
        mtext("Generations", side=1, line=2, cex=1)
        mtext("Acceptance ratio", side=2, line=2, cex=1)

        par(old.par)
    }

    if( is.list(out$phy[[1]]) ){
        return( list(accept.ratio=accept, mix.phylo=mix.phylo) )
    } else{
        return( accept.ratio=accept )
    }
}
