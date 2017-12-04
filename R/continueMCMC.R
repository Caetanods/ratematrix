##' Function to continue an unfinished MCMC chain or to append new generation to a previously finished MCMC. It works by reading the last state of the chain and the tunning parameters of the previous chain, then continuing it from this step.
##'
##' The function will append the new generations to the files created by the prior run of the 'ratematrixMCMC' function. The function will search for these files in the same directory that the previous run was performed. You should double-check if the files are in the path set in 'handle$dir'. If the files are not in the same location, you can provide a new path (relative or absolute path) to the argument 'dir'. The path provided to 'dir' will override the path pointed by 'handle$dir'. The new 'handle' output from 'continueMCMC' will have a updated total number of generations and will also update the directory path, if required.
##' @title Continue unfinished MCMC chain or add more generations
##' @param handle the output of 'ratematrixMCMC'.
##' @param add.gen number of generations to be added to the chain. This is used after the chain finished the previously set number of generations, but more generations are needed. If 'NULL' (default), the function will only continue unfinished chains.
##' @param save.handle whether to save the updated 'handle' object to the directory.
##' @param dir an optional path to the posterior chain and log files. See 'Details'.
##' @return Function will write the parameter values for each generation and the log to files. The new generations will be appended to the same files created by 'ratematrixMCMC'.
##' @export
##' @importFrom corpcor decompose.cov
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' \donttest{
##' ## Continue unfinished run.
##' data(centrarchidae)
##' ## Run the next line. Then quickly stop R from running to simulate an unexpected termination.
##' set.seed(1234)
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=4000)
##' ## Now continue the same MCMC chain. First we need to read the 'handle' from file, because
##' ##    the MCMC chain did not finished properly and the 'handle' object was not created.
##' handle ## 'object not found'.
##' handle <- readRDS(file="ratematrixMCMC.85636.mcmc.handle.rds")
##' handle.cont <- continueMCMC(handle=handle)
##' ## We can also add some generations to this same chain:
##' handle.add <- continueMCMC(handle=handle.cont, add.gen=1000)
##' }
continueMCMC <- function(handle, add.gen=NULL, save.handle=TRUE, dir=NULL){
    ## Need to use the elements of the handle to continue the mcmc.
    ## The point here is that the file for the MCMC and for the log need to be the same as the one before.
    ## So we need to open the connection to keep appending to it.
    ## Also need to read the last line of the MCMC file and use it as the starting state.
    ## No special step is needed to write to the same file. Just need to take care and not open it again or write
    ##    another header to the same file. This might need another argument for the single and multi MCMC functions.

    if( !is.null(dir) ) handle$dir <- dir

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
        } else{
            continue <- "add.gen"
            handle$gen <- handle$gen + add.gen ## Update the total gen in the handle.
        }

        singleRegimeMCMC(X=handle$data, phy=handle$phy, start=start, prior=handle$prior, gen=handle$gen, v=handle$mcmc.par$v
                       , w_sd=handle$mcmc.par$w_sd, w_mu=handle$mcmc.par$w_mu, prop=handle$mcmc.par$prop, dir=handle$dir
                       , outname=handle$outname, traits=handle$trait.names, save.handle=save.handle
                       , continue=continue, add.gen=add.gen, ID=handle$ID)
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
        
        multRegimeMCMC(X=handle$data, phy=handle$phy, start=start, prior=handle$prior, gen=handle$gen, v=handle$mcmc.par$v
                     , w_sd=handle$mcmc.par$w_sd, w_mu=handle$mcmc.par$w_mu, prop=handle$mcmc.par$prop, dir=handle$dir
                     , outname=handle$outname, regimes=handle$regime.names, traits=handle$trait.names, save.handle=save.handle
                     , continue=continue, add.gen=add.gen, ID=handle$ID)
    }

    if( save.handle ) saveRDS(handle, file = file.path(handle$dir, paste(handle$outname,".",handle$ID,".mcmc.handle.rds",sep="")) )    
    return(handle)
    
}
