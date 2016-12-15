##' Function runs a MCMC chain to estimate the posterior distribution of the evolutionary rate matrix and the phylogenetic root value given a phylogeny and trait data. Prior distribution and starting state of the chain can be chosen from one of the options provided. The function also accepts custom priors and starting states. Note that function requires output files in the chosen directory to be created and accessed. See 'Details' and 'Examples' for more information.
##'
##' Talk about function 'x' to collapse two or more regimes of a 'simmap' class tree.
##' Explain the structure of the prior so that people can build their own.
##' Explain how to fix the possible problem of function cannot write to the current directory. Need to be something simple. Maybe bullet proff it somehow.\cr
##' \cr
##' Function will check if the tree has branch lengths and will fail if branch lengths are not present. Function will also rescale the phylogeny so that the distance from every tip to the root is equal to 1. Function will give a warning if the phylogenetic tree is not ultrametric. Although all the comutations can be made with a non-ultrametric tree, the model of evolution assumes that all the species were sampled in the same point in time, so the distance from any given tip to the root of the phylogeny should be the same. This is a common assumption of trait evolution models in comparative methods. One should be careful to interpret the results based on a non-ultrametric tree.\cr
##' \cr
##'Fuction creates files with the MCMC chain. Each run of the MCMC will be identified by a unique identifier to facilitate identification and prevent the function to overwrite results when running more than one MCMC chain in the same directory. See argument 'IDlen'. The files in the directory are: 'outname.ID.loglik': the log likelihood for each generation, 'outname.ID.n.matrix': the evolutionary rate matrix n, one per line. Function will create one file for each R matrix fitted to the tree, 'outname.ID.root': the root value, one per line. \cr
##' \cr
##' Additionally it returns a list object with information from the analysis to be used by other functions. This list is refered as the 'out' parameter in those functions. The list is composed by: 'acc_ratio' numeric vector with 0 when proposal is rejected and non-zero when proposals are accepted. 1 indicates that root value was accepted, 2 and higher indicates that the first or subsequent matrices were updated; 'run_time' in seconds; 'k' the number of matrices fitted to the tree; 'p' the number of traits in the analysis; 'ID' the identifier of the run; 'dir' directory were output files were saved; 'outname' the name of the chain, appended to the names of the files; 'trait.names' A vector of names of the traits in the same order as the rows of the R matrix, can be used as the argument 'leg' for the plotting function 'make.grid.plot'; 'data' the original data for the tips; 'phy' the phylogeny; 'prior' the list of prior functions; 'start' the list of starting parameters for the MCMC run; 'gen' the number of generations of the MCMC.
##' @title Run the MCMC chain for the evolutionary rate matrix model.
##' @param data a matrix with the data. Species names need to be provided as rownames (rownames(data) == phy$tip.label). Each column is a different trait, colnames are used as the names for the traits. If not provided, the function will use default names for the traits.
##' @param phy a phylogeny of the class "simmap" with the mapped regimes. The number of evolutionary rate matrices fitted to the phylogeny is equal to the number of regimes in phy. Regime names will also be used. See 'Details'.
##' @param prior the prior densities for the MCMC. Must be one of (i) "uniform", (ii) "empirical_mean" (the default), (iii) the output of the "makePriorSeparation" function or (iv) a list of functions. See 'Details'.
##' @param start the starting state for the MCMC chain. Must be one of (i) "prior_sample" (the default), (ii) "mle", (iii) a list object. See 'Details'.
##' @param gen number of generations for the chain.
##' @param v value for the degrees of freedom parameter of the inverse-Wishart proposal distribution for the correlation matrix.
##' @param w_sd value for the width of the uniform proposal distribution for the vector of standard deviations.
##' @param w_mu value for the width of the uniform proposal distribution for the vector of root values (phylogenetic mean).
##' @param prop the proposal frequencies for each parameter of the model (default is 'c(0.025,0.975)'). This needs to be a numeric vector of length 2. Each value need to be between 0 and 1 and the sum of the vector equal to be 1. These values are the probability that the phylogenetic mean or the set of evolutionary rate matrices will be updated at each step of the MCMC chain, respectively.
##' @param chunk number of generations that the MCMC chain will be kept in the computer memory before writing to file. Larger values will use more RAM memory.
##' @param dir path of the directory to write the files (default is 'NULL'). If 'NULL' then function will write files to current directory (check 'getwd()'). If directory does not exist, then function will create a new directory. If function does not have permission (or fail) to write files to the pointed directory, then program will stop with an error message.
##' @param outname name for the MCMC chain (default is 'ratematrixMCMC'). Name will be used in all the files alongside a unique ID of numbers with length of 'IDlen'.
##' @param IDlen length of digits of the numeric identifier used to name output files (default is 5).
##' @param singlerate whether the function should fit a single regime and ignore the number of regimes painted to the tree (default is FALSE).
##' @param rescaletree whether the function will rescale the phylogenetic tree so that the depth from the tips to the root is equal to 1. Default is TRUE.
##' @return Function returns a list with the details of the MCMC run and write the MCMC output to the directory in a series of files.
##' @author daniel
##' @export
##' @importFrom mvMORPH mvBM
##' @importFrom corpcor decompose.cov
##' @importFrom ape is.ultrametric
##' @importFrom phytools rescaleSimmap
ratematrixMCMC <- function(data, phy, prior="empirical_mean", start="prior_sample", gen, v=50, w_sd=0.5, w_mu=0.5, prop=c(0.025,0.975), chunk=gen/100, dir=NULL, outname="ratematrixMCMC", IDlen=5, singlerate=FALSE, rescaletree=TRUE){

    ## #######################
    ## Block to check arguments, give warnings and etc.
    if( class(data) == "data.frame" ) data <- as.matrix( data )

    cat("\n")

    ## Inform that default options are being used:
    if( prior == "empirical_mean" ) cat("Using default prior. \n")
    if( start == "prior_sample" ) cat("Using default starting point. \n")
    if( v == 50 && w_sd == 0.5 && w_mu == 0.5 && prop[1] == 0.025 ) cat("Using default proposal settings. \n")

    ## Check if 'phy' is a single phylogeny or a list of phylogenies.
    if( is.list(phy[[1]]) ){ ## Is a list of phylogenies.
        ## Check if the tree is ultrametric, also rescale the tree if needed.
        ultra <- sapply(phy, is.ultrametric)
        if( !sum(ultra)==length(ultra) ) warning("Some (or all) phylogenetic tree are not ultrametric. Continuing analysis. Please check 'details'.")
        if( rescaletree ){
            cat("Rescaling phylogenetic trees so that depth is equal to 1.\n")
            phy <- lapply(phy, function(x) rescaleSimmap(x, model="depth", 1) )
        }
        ## Check if the phylogeny is of 'simmap' class.
        check.simmap <- sapply(phy, function(x) inherits(x, what="simmap") )
        if( !sum(check.simmap)==length(check.simmap) ){
            cat('Some (or all) of the phylogenetic tree are not of class "simmap". Fitting a sigle rate regime to the tree. \n')
            no_phymap <- TRUE
        } else{
            no_phymap <- FALSE
        }
    } else{ ## Is a single phylogeny.
        ## Check if the tree is ultrametric, also rescale the tree if needed.
        if( !is.ultrametric(phy) ) warning("Phylogenetic tree is not ultrametric. Continuing analysis. Please check 'details'.")
        if( rescaletree ){
            cat("Rescaling phylogenetic tree so that depth is equal to 1.\n")
            phy <- rescaleSimmap(phy, model="depth", 1)
        }
        ## Check if the phylogeny is of 'simmap' class.
        if( !inherits(phy, what="simmap") ){
            cat('phy is not of class "simmap". Fitting a sigle rate regime to the tree. \n')
            no_phymap <- TRUE
        } else{
            no_phymap <- FALSE
        }
    }

    ## Make a quick check down the road to see if the prior is working.
    if( !inherits(prior, what="ratematrix_prior_function") ){
        if( inherits(prior, what="character") ){
            if( !prior %in% c("uniform","empirical_mean") ) stop("prior option was not recognized. Check details for the 'prior' argument.")
        } else{
            if( inherits(prior, what="list") ){
                warning("MCMC chain using custom prior as speficied by argument 'prior'.")
            } else{
                stop("Custom prior need to be a list. Check function details.")
            }
        }
    }

    ## Make a quick check down the road to see if the start state is valid.
    if( inherits(start, what="character") ){
        if( !start %in% c("prior_sample","mle") ) stop("start state option was not recognized. Check details for the 'start' argument.")
    } else{
        if( inherits(start, what="list") || inherits(start, what="ratematrix_prior_sample")){
            cat("MCMC chain using custom starting point as speficied by argument 'start'.\n")
        } else{
            stop("Custom start state need to be a list. Check function details.")
        }
    }

    if( !class(gen) == "numeric" ) stop('gen need to be a "numeric" value.')

    ## #######################
    ## Block to create the directory for the output:
    if( is.null(dir) ){
        dir <- "."
        local <- getwd()
        cat( paste("Output files saved to current working directory: ", local, "\n", sep="" ) )
    } else{
        dir.create(file.path(dir), showWarnings = FALSE) ## This line will not modify the previous directory, so great.
        cat( paste("Output files saved to user defined directory: ", dir, "\n", sep="" ) )
    }

    ## #######################
    ## Block to set the regime and trait names:
    if( is.null( colnames(data) ) ){
        trait.names <- paste("trait_", 1:ncol(data), sep="")
    } else{
        trait.names <- colnames(data)
    }
    ## First check if analysis will use regimes.
    if( !no_phymap || singlerate ){
        if( is.list(phy[[1]]) ){ ## Check if phy is a list of phylo.
            if( is.null( colnames(phy[[1]]$mapped.edge) ) ){
                regime.names <- paste("regime_", 1:ncol(phy[[1]]$mapped.edge), sep="")
            } else{
                regime.names <- colnames(phy[[1]]$mapped.edge)
            }
        } else{
            if( is.null( colnames(phy$mapped.edge) ) ){
                regime.names <- paste("regime_", 1:ncol(phy$mapped.edge), sep="")
            } else{
                regime.names <- colnames(phy$mapped.edge)
            }
        }
    }
    

    ## #######################
    ## Block to set the analysis. First division is whether one or more regimes are fitted to the tree.
    if( no_phymap || singlerate ){
        
        r <- ncol( data )

        ## #######################
        ## Block to generate priors.
        prior_run <- prior
        if( inherits(prior, what="character") ){
            if(prior == "uniform"){
                max.data <- apply(data, 2, max)
                min.data <- apply(data, 2, min)
                mean.data <- colMeans(data)
                lower.bound.data <- mean.data - ( (mean.data - min.data) * 10 )
                higher.bound.data <- mean.data + ( (max.data - mean.data) * 10 )
                par.mu <- cbind( lower.bound.data, higher.bound.data )
                par.sd <- c(0,100)
                prior_run <- makePriorSeparation(r=r, p=1, par.mu=par.mu, par.sd=par.sd )
            }
            if(prior == "empirical_mean"){
                mn <- colMeans(data)
                ssd <- apply(data, 2, sd)
                par.mu <- as.matrix( cbind(mn, ssd) )
                par.sd <- c(0,100)
                prior_run <- makePriorSeparation(r=r, p=1, den.mu="norm", par.mu=par.mu, par.sd=par.sd)
            }
        }
        
        ## #######################
        ## Block to generate start point.
        start_run <- start
        if( inherits(start, what="character") ){
            if(start == "prior_sample"){
                start_run <- samplePriorSeparation(n=1, prior=prior_run, sample.sd=FALSE)
            }
            if(start == "mle"){ ## This will break if the phylogeny is a list.
                cat( "Optimizing likelihood for the starting value of the MCMC.\n")
                if( is.list(phy[[1]]) ){ ## Is a list of phylogenies.
                    rr <- sample(1:length(phy), size=1)
                    cat( paste("Using phylogeny number ", rr, " to estimate the MLE.\n", sep="") )
                    phy.sample <- phy[[rr]]
                    mle.fit <- mvBM(tree=phy.sample, data=data, model="BM1", method="pic", echo=FALSE)
                } else{                      
                    mle.fit <- mvBM(tree=phy, data=data, model="BM1", method="pic", echo=FALSE)
                }
                cat("\n")
                decomp.R <- decompose.cov( mle.fit$sigma )
                start_run <- list()
                start_run$mu <- as.numeric( mle.fit$theta )
                start_run$matrix <- unname( as.matrix( decomp.R$r ) )
                start_run$sd <- as.numeric( sqrt(decomp.R$v) )
            }
        }
        
        out_single <- singleRegimeMCMC(X=data, phy=phy, start=start_run, prior=prior_run, gen=gen, v=v, w_sd=w_sd, w_mu=w_mu
                                     , prop=prop, chunk=chunk, dir=dir, outname=outname, IDlen=IDlen, traits=trait.names)
        return( out_single )
        
    } else{

        ## Check if 'phy' is a single phylogeny or a list of phylogenies.
        if( is.list(phy[[1]]) ){ ## Is a list of phylogenies.
            p <- ncol( phy[[1]]$mapped.edge ) ## Multiple regimes.
        } else{ ## Is a single phylogeny.
            p <- ncol( phy$mapped.edge ) ## Multiple regimes.
        }
        
        r <- ncol( data )

        ## #######################
        ## Block to generate priors.
        prior_run <- prior
        if( inherits(prior, what="character") ){
            if(prior == "uniform"){
                max.data <- apply(data, 2, max)
                min.data <- apply(data, 2, min)
                mean.data <- colMeans(data)
                lower.bound.data <- mean.data - ( (mean.data - min.data) * 10 )
                higher.bound.data <- mean.data + ( (max.data - mean.data) * 10 )
                par.mu <- cbind( lower.bound.data, higher.bound.data )
                rep.sd.regime <- rep(c(0,100), times=p)
                par.sd <- matrix(data=rep.sd.regime, nrow=p, ncol=2, byrow=TRUE)
                prior_run <- makePriorSeparation(r=r, p=p, par.mu=par.mu, par.sd=par.sd )
            }
            if(prior == "empirical_mean"){
                mn <- colMeans(data)
                ssd <- apply(data, 2, sd)
                par.mu <- as.matrix( cbind(mn, ssd) )
                rep.sd.regime <- rep(c(0,100), times=p)
                par.sd <- matrix(data=rep.sd.regime, nrow=p, ncol=2, byrow=TRUE)
                prior_run <- makePriorSeparation(r=r, p=p, den.mu="norm", par.mu=par.mu, par.sd=par.sd)
            }
        }
        
        ## #######################
        ## Block to generate start point.
        start_run <- start
        if( inherits(start, what="character") ){
            if(start == "prior_sample"){
                start_run <- samplePriorSeparation(n=1, prior=prior_run, sample.sd=FALSE)
            }
            if(start == "mle"){ ## Need to deal with the list of matrices here.
                cat( "Optimizing likelihood for the starting value of the MCMC.\n")
                if( is.list(phy[[1]]) ){ ## Is a list of phylogenies.
                    rr <- sample(1:length(phy), size=1)
                    cat( paste("Using phylogeny number ", rr, " to estimate the MLE.\n", sep="") )
                    phy.sample <- phy[[r]]
                    mle.fit <- mvBM(tree=phy.sample, data=data, model="BMM", method="rpf", echo=FALSE)
                } else{
                    mle.fit <- mvBM(tree=phy, data=data, model="BMM", method="rpf", echo=FALSE)
                }
                cat( "\n")
                decomp.r <- list()
                decomp.sd <- list()
                for( i in 1:p ){
                    decomp.R <- decompose.cov( mle.fit$sigma[,,i] )
                    decomp.r[[i]] <- unname( as.matrix( decomp.R$r ) )
                    decomp.sd[[i]] <- as.numeric( sqrt(decomp.R$v) )
                }
                start_run <- list()
                start_run$mu <- as.numeric( mle.fit$theta )
                start_run$matrix <- decomp.r
                start_run$sd <- decomp.sd
            }
        }
        
        out_mult <- multRegimeMCMC(X=data, phy=phy, start=start_run, prior=prior_run, gen=gen, v=v, w_sd=w_sd, w_mu=w_mu
                                 , prop=prop, chunk=chunk, dir=dir, outname=outname, IDlen=IDlen, regimes=regime.names, traits=trait.names)
        return( out_mult )
        
    }
    
}
