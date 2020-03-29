##' Runs a Markov chain Monte Carlo (MCMC) chain to estimate the posterior distribution of a single
##'    evolutionary rate matrix (R) fitted to the phylogeny.
##'
##' MCMC using the inverse-Wishart as a proposal distribution for the covariance matrix and a simple sliding window for the phylogenetic mean in the random walk Metropolis-Hastings algorithm. Function writes to file. Creates an unique identifier for the files to prevent the user to append one MCMC with other working in the same directory. A future version will be able to incorporate measument error originated from population-level data.
##' @title MCMC for a single evolutionary rate matrix.
##' @param X matrix. A matrix with the data. 'rownames(X) == phy$tip.label'.
##' @param phy simmap phylo. A phylogeny of the class "simmap" from the package 'phytools'. Function uses the location information for a number of traits equal to the number of fitted matrices.
##' @param start list. Element [[1]] is the starting value for the phylogenetic mean and element [[2]] is the starting value for the R matrix. Different from 'multi.R.iwish.mcmc.R' this is a single matrix and not a list.
##' @param prior list. Produced by the output of the function 'make.prior.barnard' or 'make.prior.diwish'. First element of the list [[1]] is a prior function for the log density of the phylogenetic mean and the second element [[2]] is a prior function for the evolutionary rate matrix (R).
##' @param gen numeric. Number of generations of the MCMC.
##' @param v numeric. Degrees of freedom parameter for the inverse-Wishart proposal distribution for the evolutionary rate matrix.
##' @param w_sd numeric. Width of the uniform sliding window proposal for the vector of standard deviations.
##' @param w_mu numeric. Width of the uniform sliding window proposal for the vector of phylogenetic means.
##' @param prop vector. The proposal frequencies. Vector with two elements (each between 0 and 1). First is the probability that the phylogenetic mean will be sampled for a proposal step at each generation, second is the probability that the evolutionary rate matrix will be updated instead.
##' @param dir string. Directory to write the files, absolute or relative path. If 'NULL' then output is written to the directory where R is running (see 'getwd()'). If a directory path is given, then function will test if the directory exists and use it. If directiory does not exists the function will try to create one.
##' @param outname string. Name pasted to the files. Name of the output files will start with 'outname'.
##' @param IDlen numeric. Set the length of the unique numeric identifier pasted to the names of all output files. This is set to prevent that multiple runs with the same 'outname' running in the same directory will be lost.Default value of 5 numbers, something between 5 and 10 numbers should be good enough. IDs are generated randomly using the function 'sample'.
##' @param traits 
##' @param save.handle 
##' @param continue 
##' @param add.gen 
##' @param ID 
##' @return Fuction creates files with the MCMC chain. Each run of the MCMC will be identified by a unique identifier to facilitate identification and prevent the function to overwrite results when running more than one MCMC chain in the same directory. See argument 'IDlen'. The files in the directory are: 'outname.ID.loglik': the log likelihood for each generation, 'outname.ID.matrix': the evolutionary rate matrix, one per line, 'outname.ID.root': the root value, one per line. \cr
##' \cr
##' Additionally it returns a list object with information from the analysis to be used by other functions. This list is refered as the 'out' parameter in those functions. The list is composed by: 'acc_ratio' numeric vector with 0 when proposal is rejected and non-zero when proposals are accepted. 1 indicates that root value was accepted, 2 indicates that the evolutionary rate matrix was updated; 'run_time' in seconds; 'k' the number of matrices fitted to the tree. This value will always be 1 for this function, by see 'multi.R.iwish.mcmc'; 'p' the number of traits in the analysis; 'ID' the identifier of the run; 'dir' directory were output files were saved; 'outname' the name of the chain, appended to the names of the files; 'trait.names' A vector of names of the traits in the same order as the rows of the R matrix, can be used as the argument 'leg' for the plotting function 'make.grid.plot'; 'data' the original data for the tips; 'phy' the phylogeny; 'prior' the list of prior functions; 'start' the list of starting parameters for the MCMC run; 'gen' the number of generations of the MCMC.
##' @importFrom geiger treedata
##' @importFrom corpcor decompose.cov
##' @importFrom corpcor rebuild.cov
##' @noRd
singleRegimeMCMC <- function(X, phy, start, prior, gen, v, w_sd, w_mu, prop=c(0.1,0.45,0.45), dir=NULL, outname="single_R_fast", IDlen=5, traits, save.handle, continue=NULL, add.gen, ID=NULL){

    ## Save the 'mcmc.par' list for the mcmc.handle:
    mcmc.par <- list()
    mcmc.par$v <- v
    mcmc.par$w_sd <- w_sd
    mcmc.par$w_mu <- w_mu
    mcmc.par$prop <- prop
    
    ## Creates data and chain cache:
    cache.data <- list()
    cache.chain <- list()
    cache.data$k <- ncol(X) ## Number of traits.
    cache.data$X <- X
    cache.data$prop <- prop ## Necessary to allow for control of sd and corr proposals.

    ## Creates MCMC chain cache:
    ## Here trying to initialize the chain cache with the correct number of elements in the list.
    ## Note that the length is dependent on the 'chunk' and not the total of the 'gen'.
    ## This is a better approach to the loop.
    cache.chain$chain <- start ## Starting value for the chain.
    cache.chain$chain[[4]] <- rebuild.cov(r=stats::cov2cor(start[[2]]), v=start[[3]]^2)
    
    phy_type <- check_phy_list( phy )
    if( phy_type ){ ## The problem here is that a 'phylo' is also a list. So this checks if the first element is a list.
        cache.data$n <- length(phy[[1]]$tip.label) ## Number of tips.
        n.phy <- length(phy)
        init.phylo <- phy[[ sample(1:n.phy, size=1) ]]
        cache.chain$lik <- logLikSingleRegime(data=cache.data, chain=cache.chain, phy=init.phylo
                                               , root=as.vector(cache.chain$chain[[1]])
                                               , R=cache.chain$chain[[4]]) ## Lik start value.
    }

    if( !phy_type ){ ## There is only one phylogeny.
        cache.data$n <- length(phy$tip.label)
        cache.chain$lik <- logLikSingleRegime(data=cache.data, chain=cache.chain, phy=phy
                                               , root=as.vector(cache.chain$chain[[1]])
                                               , R=cache.chain$chain[[4]]) ## Lik start value.
    }
    cat( paste("Starting point log-likelihood: ", cache.chain$lik, "\n", sep="") )
    
    cache.chain$curr.root.prior <- prior[[1]](cache.chain$chain[[1]]) ## Prior log lik starting value.
    cache.chain$curr.r.prior <- prior[[2]](cache.chain$chain[[4]]) ## Prior log lik starting value.

    ## Will need to keep track of the Jacobian for the correlation matrix.
    decom <- decompose.cov( cache.chain$chain[[2]] )
    cache.chain$curr.r.jacobian <- sum( sapply(1:cache.data$k, function(x) log( decom$v[x]) ) ) * log( (cache.data$k-1)/2 )
    
    cache.chain$curr.sd.prior <- prior[[3]](cache.chain$chain[[3]]) ## Prior log lik starting value.

    ## Need to check if this is a continuing MCMC before creating new ID and files:
    if( is.null(continue) ){
        
        ## Generate identifier:
        ID <- paste( sample(x=1:9, size=IDlen, replace=TRUE), collapse="")

        ## Open files to write:
        files <- list( file(file.path(dir, paste(outname,".",ID,".mcmc",sep="")), open="a"),
                      file(file.path(dir, paste(outname,".",ID,".log",sep="")), open="a")
                      )
        
        header <- vector(mode="character")
        for( i in 1:cache.data$k ){
            for( j in 1:cache.data$k ){
                header <- c( header, paste("regime.", i, j, "; ", sep="") )
            }
        }
        
        ## Traitname need also to be changed for the name of the trait.
        header <- c( header, paste("trait.", 1:(cache.data$k-1), "; ", sep=""), paste("trait.", cache.data$k, sep=""), "\n")
        cat(header, sep="", file=files[[1]], append=TRUE) ## Write the header to the file.

        ## Header for the log file:
        cat("accepted; matrix.corr; matrix.sd; root; which.phylo; log.lik \n", sep="", file=files[[2]], append=TRUE)
    } else{
        ## Need to create the files object with the path to the existing files. But do not open again.
        files <- list( file.path(dir, paste(outname,".",ID,".mcmc",sep="") )
                    , file.path(dir, paste(outname,".",ID,".log",sep="") ) )
    }

    ## Build the update.function list:
    ## Now this will have two options. This is the part that the function needs to be updated.
    if( phy_type ){
        cat("MCMC chain using multiple trees/regime configurations.\n")
        update.function <- list( function(...) makePropMeanList(..., n.phy=n.phy), function(...) makePropSingleSigmaList(..., n.phy=n.phy) )
    }
    if( !phy_type ){ ## There is only one phylogeny.
        cat("MCMC chain using a single tree/regime configuration.\n")
        update.function <- list( makePropMean, makePropSingleSigma )
    }

    ## Before running need to exclude the generations already done if continuing.
    ## Also add the option to do additional generations.
    if( is.null(continue) ){
        cat( paste("Start MCMC run ", outname, ".", ID, " with ", gen, " generations.\n", sep="") )
    } else{
        if( continue == "continue" ){
            cat( paste("Continue previous MCMC run ", outname, ".", ID, " for ", add.gen, " generations for a total of ", gen, " generations.\n", sep="") )
            gen <- add.gen
        }
        if( continue == "add.gen" ){
            cat( paste("Adding ", add.gen, " generations to previous MCMC run ", outname, ".", ID, "\n", sep="") )
            gen <- add.gen
        }
    }
    
    ## Start counter for the acceptance ratio and loglik.
    count <- 2

    ## Save the handle object:
    if( save.handle ){
        out <- list(k = cache.data$k, p = 1, ID = ID, dir = dir, outname = outname, trait.names = traits
                  , data = X, phy = phy, prior = prior, start = start, gen = gen
                  , mcmc.par = mcmc.par)
        class( out ) <- "ratematrix_single_mcmc"
        saveRDS(out, file = file.path(dir, paste(outname,".",ID,".mcmc.handle.rds",sep="")) )
    }

    ## Create vector of probabilities for the sample of root vs. matrix:
    if( prop[1] > 1 ) stop("Wrong probabilities set to 'prop' argument. Check help page.")
    mcmc.prop <- c(prop[1], 1-prop[1]) / sum( c(prop[1], 1-prop[1]) )
    ## Loop for the whole MCMC
    for(i in 2:gen ){

        ## Proposals will be sampled given the 'prop' vector of probabilities.

        ## #########################################
        ## Sample which parameter is updated:
        ## 'prop' is a vector of probabilities for 'update.function' 1 or 2.
        ## 1 = phylo root and 2 = R matrix.
        up <- sample(x = c(1,2), size = 1, prob = mcmc.prop)
        ## #########################################

        ## #########################################
        ## Update and accept reject steps:
        cache.chain <- update.function[[up]](cache.data, cache.chain, prior, w_sd, w_mu, v, files, phy)
        ## #########################################

        ## #########################################
        ## Write to file:
        writeToFile(files, cache.chain)
        ## #########################################
        
    }

    ## Close the connections:
    if( is.null(continue) ) lapply(files, close)

    cat( paste("Finished MCMC run ", outname, ".", ID, "\n", sep="") )

    ## Returns 'p = 1' to indentify the results as a single R matrix fitted to the data.
    ## Returns the data, phylogeny, priors and start point to work with other functions.
    out <- list(k = cache.data$k, p = 1, ID = ID, dir = dir, outname = outname
              , trait.names = traits, data = X, phy = phy, prior = prior, start = start, gen = gen
               , mcmc.par=mcmc.par)
    class( out ) <- "ratematrix_single_mcmc"
    return( out )
}
