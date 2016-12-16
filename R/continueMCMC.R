##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Make the thing!
##' @param handle 
##' @param add.gen 
##' @param save.handle 
##' @return The thing!
##' @export
##' @importFrom corpcor decompose.cov
##' @author daniel
continueMCMC <- function(handle, add.gen=NULL, save.handle=TRUE){
    ## Need to use the elements of the handle to continue the mcmc.
    ## The point here is that the file for the MCMC and for the log need to be the same as the one before.
    ## So we need to open the connection to keep appending to it.
    ## Also need to read the last line of the MCMC file and use it as the starting state.
    ## No special step is needed to write to the same file. Just need to take care and not open it again or write
    ##    another header to the same file. This might need another argument for the single and multi MCMC functions.

    if( handle$p == 1 ){

        ## Read the MCMC output in file and grab the last state as the starting state for the run:
        ## 'readMCMC' was not made to do this. So it is not behaving well. Will need a different strategy.

        ## In this version the posterior is in a single file.
        mcmc <- read_lines( file=file.path(handle$dir, paste(handle$outname, ".", handle$ID, ".mcmc", sep="")) )
        total <- length(mcmc)-1
        last <- as.numeric( strsplit(x=mcmc[length(mcmc)], split=";", fixed=TRUE)[[1]] )
        
        ## Define the columns correspondent to the matrix:
        init <- 1
        end <- handle$k^2

        ## Creates the start object:
        start <- list()
        start$root <- as.numeric(last)[ (end+1):(end+handle$k) ]
        start$matrix <- matrix( as.numeric(last)[init:end], nrow=handle$k )
        R <- decompose.cov( start$matrix )
        start$sd <- sqrt(R$v)
        
        if( is.null(add.gen) ){
            continue <- "continue"
            add.gen <- handle$gen - total
            if( add.gen == 1 || add.gen < 0 ) stop("This MCMC run is complete. Use 'add.gen' to add more generations.\n")
            ## Check if 'chunk' is a divisible of the generation number, if not, adjust.
        } else{
            continue <- "add.gen"
            handle$gen <- handle$gen + add.gen ## Update the total gen in the handle.
        }

        ## Checks if the 'chunk' is correct given add.gen.
        if( add.gen < handle$mcmc.par$chunk ){
            handle$mcmc.par$chunk <- add.gen
        }
        if( add.gen %% handle$mcmc.par$chunk > 0 ){ ## Check the remainder of the division
            add.gen <- add.gen - (add.gen %% handle$mcmc.par$chunk)
        }

        singleRegimeMCMC(X=handle$data, phy=handle$phy, start=start, prior=handle$prior, gen=handle$gen, v=handle$mcmc.par$v
                       , w_sd=handle$mcmc.par$w_sd, w_mu=handle$mcmc.par$w_mu, prop=handle$mcmc.par$prop, chunk=handle$mcmc.par$chunk
                       , dir=handle$dir, outname=handle$outname, traits=handle$trait.names, save.handle=save.handle, continue=continue, add.gen=add.gen
                        , ID=handle$ID)

    }
    
    if( handle$p > 1 ){

        ## In this version the posterior is in a single file.
        mcmc <- read_lines( file=file.path(handle$dir, paste(handle$outname, ".", handle$ID, ".mcmc", sep="")) )
        total <- length(mcmc)-1
        last <- as.numeric( strsplit(x=mcmc[length(mcmc)], split=";", fixed=TRUE)[[1]] )
        
        ## Define the columns correspondent to the matrix:
        init <- seq(from=1, to=(handle$k^2)*handle$p, by=handle$k^2)
        end <- rev( seq(from=(handle$k^2)*handle$p, to=1, by=-handle$k^2) )

        ## Creates the start object:
        start <- list()
        
        ## Now find the root values.
        start$root <- as.numeric(last)[ (end[handle$p]+1):(end[handle$p]+handle$k) ]
        
        start$matrix <- list()
        for( i in 1:handle$p ){
            start$matrix[[i]] <- matrix( as.numeric(last)[init[i]:end[i]], nrow=handle$k )
        }
        
        start$sd <- list()
        for( i in 1:handle$p ){
            R <- decompose.cov( start$matrix[[i]] )
            start$sd[[i]] <- sqrt(R$v)
        }

        if( is.null(add.gen) ){
            continue <- "continue"
            add.gen <- handle$gen - total
            if( add.gen == 1 || add.gen < 0 ) stop("This MCMC run is complete. Use 'add.gen' to add more generations.\n")
        } else{
            continue <- "add.gen"
            handle$gen <- handle$gen + add.gen ## Update the total gen in the handle.
        }
        
        ## Checks if the 'chunk' is correct given add.gen.
        if( add.gen < handle$mcmc.par$chunk ){
            handle$mcmc.par$chunk <- add.gen
        }
        if( add.gen %% handle$mcmc.par$chunk > 0 ){ ## Check the remainder of the division
            add.gen <- add.gen - (add.gen %% handle$mcmc.par$chunk)
        }
        
        multRegimeMCMC(X=handle$data, phy=handle$phy, start=start, prior=handle$prior, gen=handle$gen, v=handle$mcmc.par$v
                     , w_sd=handle$mcmc.par$w_sd, w_mu=handle$mcmc.par$w_mu, prop=handle$mcmc.par$prop, chunk=handle$mcmc.par$chunk
                     , dir=handle$dir, outname=handle$outname, regimes=handle$regime.names, traits=handle$trait.names, save.handle=save.handle
                     , continue=continue, add.gen=add.gen, ID=handle$ID)
        
    }

    if( save.handle ) saveRDS(handle, file = file.path(handle$dir, paste(handle$outname,".",handle$ID,".mcmc.handle.rds",sep="")) )    
    return(handle)
    
}
