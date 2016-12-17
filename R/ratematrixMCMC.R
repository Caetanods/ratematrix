##' Function runs a MCMC chain to estimate the posterior distribution of the evolutionary rate matrix (R) and the root value (phylogenetic mean). Prior distribution and starting state for the chain can be chosen among pre-defined options or manually set by the user using accompanying functions. Please note that the function will write files to the directory.
##'
##' The MCMC chain works by proposing values for the evolutionary rate matrices (R) fitted to the tree and the vector of root values (or phylogenetic mean). The proposal for the R matrices works by separating the variance-covariance matrix into a correlation matrix and a vector of standard deviations and making independent proposals for each. This scheme is called the 'separation strategy' and significantly improves the mix of the chain and also provide a intuitive distinction between the evolutionary correlation among the traits (correlation matrix) and the rates of evolution (standard deviation vector). The proposal for the root values are made all in a single step. \cr
##' \cr
##' The function will print a series of messages to the screen. Those provide details of the setup of the chain, the unique identifier for the files and the log-likelihood of the starting value of the chain. Up to now these messages cannot be disabled. \cr
##' \cr
##' SAMPLE OF TREES: The MCMC chain can integrate the phylogenetic uncertainty or the uncertainty in the rate regimes by randomly sampling a phylogenetic tree from a list of trees. To activate this option, provide a list of 'simmap' or 'phylo' trees as the 'phy' argument. The MCMC will randomly sampled a tree each time the likelihood of the proposal is evaluated. Check the 'logAnalizer' function for more information. \cr
##' \cr
##' MCMC DOES NOT START: It is possible that the starting point shows a very low likelihood value, resulting in the collapse of the chain. This might be a result of a random sample from a very unlikely region of the prior. We suggest that another sample of the prior is taken or, if this does not solve the issue, that the starting point be the maximum likelihood estimate or set manually. \cr
##' \cr
##' MCMC DOES NOT CONVERGE OR MIX: If the MCMC is taking too long to converge (something between 500000 and 1500000 is OK) then the parameters of the chain might not be good for your data. First check the 'logAnalyzer' function. The recommended acceptance ratio is ~ 0.24, if it is too high, then the step size of the proposals might be too small, try increasing the step size. In contrast, low acceptance ratio might be due to step sizes too large. Try to decrease the size of the steps. If the effective sample size (ESS) for the chain (see 'checkConvergence' function) is to low for some parameter, then try to increase the proportion of times that the parameter is proposed in the MCMC. Future versions of the package will provide an adaptive MCMC step to automatically tun the sampler. There is no "magic default value" for this, it is normal to adjust the MCMC sampler in function of the characteristics of the data. \cr
##' \cr
##' CANNOT FIND THE POSTERIOR: The function writes the posterior into two files: The '.log' file has the log-likelihood and information about which phylogeny was used, which parameter was proposed and whether the step was accepted or not. The '.mcmc' file has the posterior for the parameters of the model. Those are identified by a name for the chain set by "outname" and an unique series of numbers set by "IDlen". Note that you will need the handle object provided as the output for the function (or saved to the directory if 'save.handle' is TRUE) to be able to load, plot and analyze the posterior distribution. \cr
##' \cr
##' TREE BRANCH LENGTHS: Function will check if the tree has branch lengths and will fail if branch lengths are not present. Function will show a warning if the phylogenetic tree is not ultrametric. Although all computations can be made with a non-ultrametric tree, the model of evolution assumes that all the species were sampled from the same slice of time. Such that the distance from any given tip to the root of the phylogeny should be the same. This is a very common assumption of trait evolution models in comparative methods. Please be careful to interpret parameter estimates in the case of a non-ultrametric tree.
##' @title Estimate the evolutionary rate matrix using Markov-chain Monte Carlo
##' @param data a matrix with the data. Species names need to be provided as rownames (rownames(data) == phy$tip.label). Each column is a different trait. Names for the columns is used as trait labels. If labels are not provided, the function will use default labels.
##' @param phy a phylogeny of the class "simmap" with the mapped regimes for two or more R regimes OR a phylogeny of the class "phylo" for a single regime. The number of evolutionary rate matrices fitted to the phylogeny is equal to the number of regimes in 'phy'. Regime names will also be used. 'phy' can also be a list of phylogenies. See 'Details'.
##' @param prior the prior densities for the MCMC. Must be one of (i) "uniform", (ii) "empirical_mean" (the default), (iii) the output of the "makePrior" function or (iv) a list of functions. See more information on 'makePrior' function.
##' @param start the starting state for the MCMC chain. Must be one of (i) "prior_sample" (the default), (ii) "mle", (iii) a list object. See more information on 'makeStart' function.
##' @param gen number of generations for the chain.
##' @param v value for the degrees of freedom parameter of the inverse-Wishart proposal distribution for the correlation matrix. Smaller values provide larger steps and larger values provide smaller steps. (Yes, it is counterintuitive.)
##' @param w_sd value for the width of the uniform proposal distribution for the vector of standard deviations.
##' @param w_mu value for the width of the uniform proposal distribution for the vector of root values (phylogenetic mean).
##' @param prop the proposal frequencies for each parameter of the model (default is 'c(0.025,0.975)'). This needs to be a numeric vector of length 2. Each value need to be between 0 and 1 and the sum of the vector equal to 1. These values are the probability that the phylogenetic mean or the set of evolutionary rate matrices will be updated at each step of the MCMC chain, respectively.
##' @param chunk number of generations that the MCMC chain will be kept in the computer memory before writing to file. Larger values will use more RAM memory. If 'chunk' is less than 'gen', then 'chunk = gen'. If 'gen' is not divisible by 'chunk', then the remainder of the integer division will be subtracted from 'gen'. (default is 'gen/100')
##' @param dir path of the directory to write the files (default is 'NULL'). If 'NULL', then function will write files to the current working directory (check 'getwd()'). If directory does not exist, then function will create it. The path can be provided both as relative or absolute. It should accept Linux, Mac and Windows path formats.
##' @param outname name for the MCMC chain (default is 'ratematrixMCMC'). Name will be used in all the files alongside a unique ID of numbers with length of 'IDlen'.
##' @param IDlen length of digits of the numeric identifier used to name output files (default is 5).
##' @param singlerate whether the function should fit a single regime regardless of the regimes painted to the tree. (default is FALSE)
##' @param rescaletree whether the function will rescale the phylogenetic tree so that the depth from the tips to the root is equal to 1. (Default is FALSE).
##' @param save.handle whether the handle for the MCMC should be saved to the directory in addition to the output files.
##' @return Function returns the 'handle' object and writes the posterior distribution and log as files in the directory (see 'dir'). The handle is a list with the details of the MCMC chain. It is composed by: *k* the number of traits; *p* the number of R regimes fitted to the tree; *ID* the unique identifier of the run; *dir* the directory where the posterior and log files were saved; *outname* the name for the chain; *trait.names* a vector with the label for the traits; *regime.names* a vector with the label for the rate regimes; *data* the data used in the analysis; *phy* a single phylogeny or the list of phylogenies; *prior* a list with the prior functions; *start* a list with the starting parameters for the chain; *gen* the number of generations for the chain; *mcmc.par* a list with the tunning parameters for the MCMC.
##' @author Daniel S. Caetano and Luke J. Harmon
##' @export
##' @importFrom mvMORPH mvBM
##' @importFrom corpcor decompose.cov
##' @importFrom ape is.ultrametric
##' @importFrom phytools rescaleSimmap
##' @seealso \code{\link{ estimateTimeMCMC }} to estimate the time for the MCMC chain, \code{\link{ readMCMC }} for reading the output files, \code{\link{ plotPrior }} for plotting the prior, \code{\link{ plotRatematrix }} and \code{\link{ plotRootValue }} for plotting the posterior,  \code{\link{ checkConvergence }} to check convergence, \code{\link{ testRatematrix }} to perform tests, and \code{\link{ logAnalizer }} to read and analyze the log file.
##' @examples
##' \donttest{
##' ## Load data
##' data(centrarchidae)
##' ## Run MCMC. This is just a very short chain.
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=1000)
##' ## Load posterior distribution, make plots and check the log.
##' posterior <- readMCMC(handle, burn=0.25, thin=1)
##' plotRatematrix(posterior)
##' plotRootValue(posterior)
##' plotPrior(handle)
##' logAnalizer(handle)
##' }
ratematrixMCMC <- function(data, phy, prior="empirical_mean", start="prior_sample", gen, v=50, w_sd=0.5, w_mu=0.5, prop=c(0.025,0.975), chunk=gen/100, dir=NULL, outname="ratematrixMCMC", IDlen=5, singlerate=FALSE, rescaletree=FALSE, save.handle=TRUE){

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

    ## Check if 'chunk' is a divisible of the generation number, if not, adjust.
    ## Default is: gen/100
    if( chunk > gen ){
        chunk <- gen
    }        
    if( gen %% chunk > 0 ){ ## Check the remainder of the division
        gen <- gen - (gen %% chunk)
    }   

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
                prior_run <- makePrior(r=r, p=1, par.mu=par.mu, par.sd=par.sd )
            }
            if(prior == "empirical_mean"){
                mn <- colMeans(data)
                ssd <- apply(data, 2, sd)
                par.mu <- as.matrix( cbind(mn, ssd) )
                par.sd <- c(0,100)
                prior_run <- makePrior(r=r, p=1, den.mu="norm", par.mu=par.mu, par.sd=par.sd)
            }
        }
        
        ## #######################
        ## Block to generate start point.
        start_run <- start
        if( inherits(start, what="character") ){
            if(start == "prior_sample"){
                start_run <- samplePrior(n=1, prior=prior_run, sample.sd=FALSE)
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

        cat(paste0("chunk", chunk, "\n"))
        cat(paste0("gen", gen, "\n"))
        out_single <- singleRegimeMCMC(X=data, phy=phy, start=start_run, prior=prior_run, gen=gen, v=v, w_sd=w_sd, w_mu=w_mu
                                     , prop=prop, chunk=chunk, dir=dir, outname=outname, IDlen=IDlen, traits=trait.names
                                      , save.handle=save.handle)
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
                prior_run <- makePrior(r=r, p=p, par.mu=par.mu, par.sd=par.sd )
            }
            if(prior == "empirical_mean"){
                mn <- colMeans(data)
                ssd <- apply(data, 2, sd)
                par.mu <- as.matrix( cbind(mn, ssd) )
                rep.sd.regime <- rep(c(0,100), times=p)
                par.sd <- matrix(data=rep.sd.regime, nrow=p, ncol=2, byrow=TRUE)
                prior_run <- makePrior(r=r, p=p, den.mu="norm", par.mu=par.mu, par.sd=par.sd)
            }
        }
        
        ## #######################
        ## Block to generate start point.
        start_run <- start
        if( inherits(start, what="character") ){
            if(start == "prior_sample"){
                start_run <- samplePrior(n=1, prior=prior_run, sample.sd=FALSE)
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
                                 , prop=prop, chunk=chunk, dir=dir, outname=outname, IDlen=IDlen, regimes=regime.names, traits=trait.names
                                   , save.handle=save.handle)
        return( out_mult )
        
    }
    
}
