##' Function runs a MCMC chain to estimate the posterior distribution of the evolutionary rate matrix (R) and the root value (phylogenetic mean). Prior distribution and starting state for the chain can be chosen among pre-defined options or manually set by the user using accompanying functions (see function 'makePrior' for more information). User NEED to provide a directory to write the files (See 'Details'). Use dir="." option to write files to the current directory or provide a name of a folder to be created.
##'
##' The MCMC chain works by proposing values for the evolutionary rate matrices (R) fitted to the tree and the vector of root values (or phylogenetic mean). The proposal for the R matrices works by separating the variance-covariance matrix into a correlation matrix and a vector of standard deviations and making independent proposals for each. This scheme is called the 'separation strategy' and significantly improves the mix of the chain and also provide a intuitive distinction between the evolutionary correlation among the traits (correlation matrix) and the rates of evolution (standard deviation vector). The proposal for the root values are made all in a single step. \cr
##' \cr
##' The function will print a series of messages to the screen. Those provide details of the setup of the chain, the unique identifier for the files and the log-likelihood of the starting value of the chain. Up to now these messages cannot be disabled. \cr
##' \cr
##' DIRECTORY TO WRITE FILES: User need to specify the directory to write the files. Use "." to write to the current directoy. Or, for example, use "MCMC_files" to create the folder "MCMC_files" in the current directory. RCran policy prohibits the package to automaticaly write to the current directory. \cr
##' \cr
##' DEFAULT PRIOR: The default prior distribution ('uniform_scaled') is composed by a uniform distribution on the root values with range equal to the range observed at the tip data. The size of the window used at each proposal step for the root values is equal to the width of the prior divided by 10 units. For the evolutionary rate matrix, this prior sets a uniform distribution on the correlations (spanning all possible correlation structures) and also a uniform distribution on the vector of standard deviations. The limits of the prior on the standard deviation is computed by first doing a quick Maximum Likelihood estimate of each trait under a single rate BM model and using the results to inform the magnitude of the rates. This default prior distribution might not be the best for your dataset. Keep in mind that the default behavior of the MCMC is to draw a starting point from the prior distribution. Please check the 'makePrior' function for more information on priors and how to make a custom prior distribution.\cr
##' \cr
##' SAMPLE OF TREES: The MCMC chain can integrate the phylogenetic uncertainty or the uncertainty in the rate regimes by randomly sampling a phylogenetic tree from a list of trees. To activate this option, provide a list of 'simmap' or 'phylo' trees as the 'phy' argument. The MCMC will randomly sample a tree each proposal step. Check the 'logAnalyzer' function for more information. \cr
##' \cr
##' MCMC DOES NOT START: It is possible that the starting point shows a very low likelihood value, resulting in the collapse of the chain. This might be a result of a random sample from a very unlikely region of the prior. We suggest that another sample of the prior is taken or that the user make a more suitable prior using the function 'makePrior'. \cr
##' \cr
##' MCMC DOES NOT CONVERGE OR MIX: If the MCMC is taking too long to converge then the parameters of the chain might not be good for your data. First check the 'logAnalyzer' function as well as the 'computeESS'. The recommended acceptance ratio is ~ 0.24, if it is too high, then the step size of the proposals might be too small, try increasing the step size. In contrast, low acceptance ratio might be due to step sizes too large. Try to decrease the size of the steps. If the effective sample size (ESS) for the chain (see 'checkConvergence' and 'computeESS' functions) is to low for some parameter, then try to increase the proportion of times that the parameter is proposed in the MCMC. \cr
##' \cr
##' CANNOT FIND THE POSTERIOR: The function writes the posterior into two files: The '.log' file has the log-likelihood and information about which phylogeny was used, which parameter was proposed and whether the step was accepted or not. The '.mcmc' file has the posterior for the parameters of the model. Those are identified by a name for the chain set by "outname" and an unique series of numbers set by "IDlen". Note that you will need the handle object provided as the output for the function (or saved to the directory if 'save.handle' is TRUE) to be able to load, plot and analyze the posterior distribution.
##' @title Estimate the evolutionary rate matrix using Markov-chain Monte Carlo
##' @param data a matrix with the data. Species names need to be provided as rownames (rownames(data) == phy$tip.label). Each column is a different trait. Names for the columns is used as trait labels. If labels are not provided, the function will use default labels.
##' @param phy a phylogeny of the class "simmap" with the mapped regimes for two or more R regimes OR a phylogeny of the class "phylo" for a single regime. The number of evolutionary rate matrices fitted to the phylogeny is equal to the number of regimes in 'phy'. Regime names will also be used. 'phy' can also be a list of phylogenies. See 'Details'.
##' @param prior the prior densities for the MCMC. Must be one of "uniform", "uniform_scaled" (the default, see 'Details'), "empirical_mean", or the output of the "makePrior" function. See more information on 'makePrior' and in the examples below.
##' @param start the starting state for the MCMC chain. Must be one of "prior_sample" (the default), "mle", or a sample from the prior generated with the "samplePrior" functions.
##' @param gen number of generations for the chain. The final number of posterior samples is equal to "(gen*burn)/thin".
##' @param burn the fraction of the chain for the burnin (not written to file). A numeric value between 0 and 1 (i.e., 0 means no burnin). The final number of posterior samples is equal to "(gen*burn)/thin".
##' @param thin the number of generations to be skipped between each sample of the posterior. The final number of posterior samples is equal to "(gen*burn)/thin".
##' @param v value for the degrees of freedom parameter of the inverse-Wishart proposal distribution for the correlation matrix. Smaller values provide larger steps and larger values provide smaller steps. (Yes, it is counterintuitive.) This needs to be a single value applied to all regimes or a vector with the same length as the number of regimes.
##' @param w_sd the multiplying factor for the multiplier proposal on the vector of standard deviations. This can be a single value to be used for the sd of all traits for all regimes or a matrix with number of columns equal to the number of regimes and number of rows equal to the number of traits. If a matrix, then each element will be used to control the correspondent width of the standard deviation.
##' @param w_mu value for the width of the sliding window proposal for the vector of root values (phylogenetic mean). This can be a single value to be used for the root value of all traits or a vector of length equal to the number of traits. If a vector, then each element will be used as the width of the proposal distribution for each trait in the same order as the columns in 'data'. When 'prior="uniform_scaled"' (the default) this parameter is computed from the data.
##' @param prop a numeric vector of length 3 with the proposal frequencies for each parameter of the model. The vector need to sum to 1. These values are the probability that the phylogenetic mean (prop[1]), the vector of standard deviations (prop[2]), and the correlation matrix (prop[3]) will be updated at each step of the MCMC chain, respectively. Default value is 'c(0.05, 0.475, 0.475)'.
##' @param dir path of the directory to write the files. Has no default value (due to RCran policy). The path can be provided both as relative or absolute. It should accept Linux, Mac and Windows path formats.
##' @param outname name for the MCMC chain (default is 'ratematrixMCMC'). Name will be used in all the files alongside a unique ID of numbers with length of 'IDlen'.
##' @param IDlen length of digits of the numeric identifier used to name output files (default is 5).
##' @param save.handle whether the handle for the MCMC should be saved to the directory in addition to the output files.
##' @return Function returns the 'handle' object and writes the posterior distribution and log as files in the directory (see 'dir'). The handle is a list with the details of the MCMC chain. It is composed by: *k* the number of traits; *p* the number of R regimes fitted to the tree; *ID* the unique identifier of the run; *dir* the directory where the posterior and log files were saved; *outname* the name for the chain; *trait.names* a vector with the label for the traits; *regime.names* a vector with the label for the rate regimes; *data* the data used in the analysis; *phy* a single phylogeny or the list of phylogenies; *prior* a list with the prior functions; *start* a list with the starting parameters for the chain; *gen* the number of generations for the chain; *mcmc.par* a list with the tunning parameters for the MCMC.
##' @author Daniel S. Caetano and Luke J. Harmon
##' @references Revell, L. J., and L. J. Harmon. 2008. Testing quantitative genetic hypotheses about the evolutionary rate matrix for continuous characters. Evolutionary Ecology Research 10:311.
##' @references Revell, L. J., and D. C. Collar. 2009. Phylogenetic Analysis of the Evolutionary Correlation Using Likelihood. Evolution 63:1090–1100.
##' @references Caetano, D. S., and L. J. Harmon. 2017. ratematrix: An R package for studying evolutionary integration among several traits on phylogenetic trees. Methods in Ecology and Evolution 8:1920–1927.
##' @references Caetano, D. S., and L. J. Harmon. 2018. Estimating Correlated Rates of Trait Evolution with Uncertainty. Systematic Biology, doi: 10.1093/sysbio/syy067.
##' @export
##' @importFrom mvMORPH mvBM
##' @importFrom corpcor decompose.cov
##' @importFrom ape is.ultrametric
##' @importFrom ape Ntip
##' @importFrom geiger fitContinuous
##' @importFrom stats coef
##' @examples
##' \donttest{
##' data( centrarchidae )
##' ## Set the limits of the uniform prior on the root based on the observed traits
##' data.range <- t( apply( centrarchidae$data, 2, range ) )
##' ## The step size for the root value can be set given the range we need to sample from:
##' w_mu <- ( data.range[,2] - data.range[,1] ) / 10
##' ## Set a reasonable value for the uniform prior distribution for the standard deviation.
##' ## Here the minimum rate for the traits is 0 and the maximum is 10 ( using 'sqrt(10)' to 
##' ##      transform to standard deviation).
##' par.sd <- cbind(c(0,0), sqrt( c(10,10) ))
##' ## The proposal step on the standard deviation is a multiplier. So 0.2 is good enough 
##' ##       for most cases.
##' w_sd <- matrix(0.2, ncol = 2, nrow = 2)
##' prior <- makePrior(r = 2, p = 2, den.mu = "unif", par.mu = data.range, den.sd = "unif"
##'                    , par.sd = par.sd)
##' ## Run multiple MCMC chains.
##' handle.list <- lapply(1:4, function(x) ratematrixMCMC(data=centrarchidae$data
##'                       , phy=centrarchidae$phy.map, prior=prior, gen=10000
##'                       , w_mu=w_mu, w_sd=w_sd, dir=tempdir()) )
##' ## Read all to a list
##' posterior.list <- lapply(handle.list, readMCMC)
##' ## Check for convergence (it might not converge with only 10000 steps)
##' checkConvergence(posterior.list)
##' ## Merge all posteriors in the list.
##' merged.posterior <- mergePosterior(posterior.list)
##' ## PLot results:
##' plotRatematrix(merged.posterior)
##' plotRootValue(merged.posterior)
##' }
ratematrixMCMC <- function(data, phy, prior="uniform_scaled", start="prior_sample", gen = 1000000, burn = 0.25, thin = 100, v=50, w_sd=2, w_mu=0.5, prop=c(0.05, 0.475, 0.475), dir=NULL, outname="ratematrixMCMC", IDlen=5, save.handle=TRUE){

    ## #######################
    ## Block to check arguments, give warnings and etc.

    ## Check burn and thin and create the vector of generations.
    ## Note here that the first generation is gen 0.
    post_seq <- seq(from = gen * burn, to = (gen-1), by = thin)
    post_seq <- post_seq + 1 ## Bounce the generation vector forward to start from 1.
    
    ## Quickly check if a directory was provide. If not return an error.
    if( is.null(dir) ) stop('Need to provide a path to write MCMC samples. Use dir="." to write files to current directory.')
    if( !inherits(dir, what="character") ) stop("Value for argument 'dir' need to be a character. See help page.")

    ## Corrects the data if necessary.
    if( class(data) == "data.frame" ) data <- as.matrix( data )

    cat("\n")

    ## Inform that default options are being used:
    if( inherits(prior, what="character") && prior == "empirical_mean" ) cat("Using old default prior. \n")
    if( inherits(prior, what="character") && prior == "uniform_scaled" ) cat("Using new default prior. \n")
    if( inherits(start, what="character") && start == "prior_sample" ) cat("Using default starting point. \n")
    if( v == 50 && w_sd == 0.2 && w_mu == 0.5 && prop[1] == 0.05 && prop[2] == 0.475 ) cat("Using default proposal settings. \n")

    ## Check if provided prior has the correct dimension:
    if( inherits(prior, what="ratematrix_prior_function") ){
        if( !prior$pars$r == ncol(data) ) stop("Wrong number of traits specified on the prior.")
    }
    
    ## Check formats for 'w_mu'. Need to delay check for 'w_sd' until we get the number of regimes.
    if( length( w_mu ) > 1 ){
        if( !length(w_mu) == ncol(data) ) stop("Length of 'w_mu' need to be 1 or equal to number of traits.")
    } else{
        w_mu <- rep(w_mu, times=ncol(data))
    }
    
    ## Check if 'phy' is a single phylogeny or a list of phylogenies.
    if( is.list(phy[[1]]) ){ ## Is a list of phylogenies.
        ## Check if the tree is ultrametric, also rescale the tree if needed.
        ultra <- sapply(phy, is.ultrametric)
        if( !sum(ultra)==length(ultra) ) warning("Some (or all) phylogenetic tree are not ultrametric. Continuing analysis. Please check 'details'.")
        ## Check if the phylogeny is of 'simmap' class.
        check.simmap <- sapply(phy, function(x) inherits(x, what="simmap") )
        if( !sum(check.simmap)==length(check.simmap) ){
            cat('Some (or all) of the phylogenetic tree are not of class "simmap". Fitting a sigle rate regime to the tree. \n')
            no_phymap <- TRUE
        } else{
            no_phymap <- FALSE
        }
        ## Check if data match the tree. If not, break and return error message.
        equaln <- sapply(phy, function(x) Ntip( x ) == nrow( data ) )
        if( !sum(equaln, na.rm=TRUE) == length(phy) ) stop("Number of species in data not equal to tree.")
        same.spp <- sapply(phy, function(x) sum(x$tip.label %in% rownames( data ), na.rm=TRUE) == nrow(data) )
        if( !sum(same.spp) == length(phy) ) stop("Rownames of data does not match the tip labels of the tree.")
        ## Reorder the data to be the same order as the tip labels.
        mm <- match(phy[[1]]$tip.label, rownames(data))
        data <- data[mm,]
        
    } else{ ## Is a single phylogeny.
        ## Check if the tree is ultrametric, also rescale the tree if needed.
        if( !is.ultrametric(phy) ) warning("Phylogenetic tree is not ultrametric. Continuing analysis. Please check 'details'.")
        ## Check if the phylogeny is of 'simmap' class.
        if( !inherits(phy, what="simmap") ){
            cat('phy is not of class "simmap". Fitting a sigle rate regime to the tree. \n')
            no_phymap <- TRUE
        } else{
            no_phymap <- FALSE
        }
        ## Check if data match the tree. If not, break and return error message.
        equaln <- Ntip( phy ) == nrow( data )
        if( !equaln ) stop("Number of species in data not equal to tree.")
        same.spp <- sum(phy$tip.label %in% rownames( data ), na.rm=TRUE) == nrow(data)
        if( !same.spp ) stop("Rownames of data does not match the tip labels of the tree.")
        ## Reorder the data to be the same order as the tip labels.
        mm <- match(phy$tip.label, rownames(data))
        data <- data[mm,]
    }

    ## Check the 'w_sd' parameter. Need to know if phylo is a list and if it is simmap.
    ## Get the number of regimes and check the prior, if provided.
    if( is.list(phy[[1]]) ){ ## Is a list of phylogenies.
        if( is.null( phy[[1]]$mapped.edge ) ){
            n_regimes <- 1
            if( inherits(prior, what="ratematrix_prior_function") ){
                if( !prior$pars$p == n_regimes ) stop("Number of regimes specified on prior does not match the phylogeny.")
            }
        } else{
            n_regimes <- ncol(phy[[1]]$mapped.edge)
            if( inherits(prior, what="ratematrix_prior_function") ){
                if( !prior$pars$p == n_regimes ) stop("Number of regimes specified on prior does not match the phylogeny.")
            }
        }
    } else{
        if( is.null( phy$mapped.edge ) ){
            n_regimes <- 1
            if( inherits(prior, what="ratematrix_prior_function") ){
                if( !prior$pars$p == n_regimes ) stop("Number of regimes specified on prior does not match the phylogeny.")
            }
        } else{
            n_regimes <- ncol(phy$mapped.edge)
            if( inherits(prior, what="ratematrix_prior_function") ){
                if( !prior$pars$p == n_regimes ) stop("Number of regimes specified on prior does not match the phylogeny.")
            }
        }
    }
    
    if( is.matrix( w_sd ) ){
        if( !ncol(w_sd) == n_regimes ) stop("ncol(w_sd) need to be equal to the number of regimes.")
        if( !nrow(w_sd) == ncol(data) ) stop("nrow(w_sd) need to be equal to the number of traits.")
    } else{
        if( n_regimes > 1 & !length( w_sd ) == 1 ) stop(" 'w_sd' need to be a single numeric value or a matrix with ncol equal to the number of regimes and nrow equal to the number of traits.")
        if( length( w_sd ) == 1 ){
            w_sd <- rep(w_sd, times=ncol(data))
        }
        if( !length( w_sd ) == ncol(data) ) stop("length of 'w_sd' need to be equal to the number of traits.")
        if( n_regimes > 1 ){ ## w_sd needs to be a matrix!
            temp_mat <- matrix(nrow=length(w_sd), ncol=n_regimes)
            for(i in 1:n_regimes) temp_mat[,i] <- w_sd
            w_sd <- temp_mat
        }
    }

    ## Make a quick check down the road to see if the prior is working.
    if( !inherits(prior, what="ratematrix_prior_function") ){
        if( inherits(prior, what="character") ){
            if( !prior %in% c("uniform","empirical_mean", "uniform_scaled") ) stop("prior option was not recognized. Check details for the 'prior' argument.")
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
        ## dir <- "."
        ## local <- getwd()
        ## cat( paste("Output files saved to current working directory: ", local, "\n", sep="" ) )
        stop('Need to provide a path to write MCMC samples. Use dir="." to write files to current directory.')
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
    if( !no_phymap ){
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
    if( no_phymap ){
        
        r <- ncol( data )

        ## Check if 'w_sd' is a matrix, then make it a vector:
        if( is.matrix( w_sd ) ){
            if( ncol(w_sd) == 1 & nrow == r ){
                w_sd <- as.numeric( w_sd[,1] )
            } else{
                stop( "Parameter w_sd needs to be a single value or a matrix with number of columns equal to the number of regimes (i.e., 1) and number of rows equal to the number of traits. See help page." )
            }
        } else{
            if( length( w_sd > 1 ) ){
                if( !length( w_sd == r ) ) stop( "Length of w_sd vector needs to be equal to the number of traits.")
            }
        }

        ## #######################
        ## Block to generate priors.
        prior_run <- prior
        if( inherits(prior, what="character") ){
            if(prior == "uniform"){
                max.data <- apply(data, 2, max)
                min.data <- apply(data, 2, min)
                par.mu <- cbind( min.data, max.data )
                par.sd <- c(0, sqrt(10)) ## Prior for the standard deviation.
                prior_run <- makePrior(r=r, p=1, par.mu=par.mu, par.sd=par.sd )
            }
            if(prior == "empirical_mean"){
                mn <- colMeans(data)
                ssd <- apply(data, 2, stats::sd) * 2
                par.mu <- as.matrix( cbind(mn, ssd) )
                par.sd <- c(0, sqrt(10)) ## Prior for the standard deviation.
                prior_run <- makePrior(r=r, p=1, den.mu="norm", par.mu=par.mu, par.sd=par.sd)
            }
            if(prior == "uniform_scaled"){
                data.range <- t( apply(data, 2, range) )
                ## Step on the root is 1/10 of the space.
                ## This overrides the user parameters!
                cat("Computing step size for root value proposal from the data. \n")
                w_mu <- ( data.range[,2] - data.range[,1] ) / 10
                cat("Guessing magnitude of rates from the data. \n")
                if( is.list(phy[[1]]) ){ ## List of phylo.
                    ## Using a single core to compute BM model.
                    fit <- lapply(1:ncol(data), function(x) fitContinuous(phy = phy[[1]], dat=data[,x], model = "BM", ncores = 1) )
                } else{ ## Only one phylo.
                    ## Using a single core to compute BM model.
                    fit <- lapply(1:ncol(data), function(x) fitContinuous(phy = phy, dat=data[,x], model = "BM", ncores = 1) )
                }
                guess.rates <- sapply(fit, function(x) coef(x)[1])
                top.sd <- sqrt( ceiling( max(guess.rates) ) * 10 )
                bottom.sd <- 0
                par.sd <- cbind(bottom.sd, top.sd)
                prior_run <- makePrior(r=r, p=1, den.mu="norm", par.mu=data.range, par.sd=par.sd)
            }
        }
        
        ## #######################
        ## Block to generate start point.
        start_run <- start
        if( inherits(start, what="character") ){
            if(start == "prior_sample"){
                cat( "Taking sample from prior as starting point. \n" )
                start_run <- samplePrior(n=1, prior=prior_run)
            }
            if(start == "mle"){ ## This will break if the phylogeny is a list.
                cat( "Optimizing likelihood for the starting value of the MCMC.\n")
                if( is.list(phy[[1]]) ){ ## Is a list of phylogenies.
                    rr <- sample(1:length(phy), size=1)
                    cat( paste("Using phylogeny number ", rr, " to estimate the MLE.\n", sep="") )
                    phy.sample <- phy[[rr]]
                    mle.fit <- mvBM(tree=phy.sample, data=data, model="BM1", method="pic", echo=FALSE, diagnostic=FALSE)
                } else{                      
                    mle.fit <- mvBM(tree=phy, data=data, model="BM1", method="pic", echo=FALSE, diagnostic=FALSE)
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
                                     , prop=prop, dir=dir, outname=outname, IDlen=IDlen, traits=trait.names
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
                par.mu <- cbind( min.data, max.data )
                rep.sd.regime <- rep(c(0,sqrt(10)), times=p)
                par.sd <- matrix(data=rep.sd.regime, nrow=p, ncol=2, byrow=TRUE)
                prior_run <- makePrior(r=r, p=p, par.mu=par.mu, par.sd=par.sd )
            }
            if(prior == "empirical_mean"){
                mn <- colMeans(data)
                ssd <- apply(data, 2, stats::sd)
                par.mu <- as.matrix( cbind(mn, ssd) )
                rep.sd.regime <- rep(c(0,sqrt(10)), times=p)
                par.sd <- matrix(data=rep.sd.regime, nrow=p, ncol=2, byrow=TRUE)
                prior_run <- makePrior(r=r, p=p, den.mu="norm", par.mu=par.mu, par.sd=par.sd)
            }
            if(prior == "uniform_scaled"){
                data.range <- t( apply(data, 2, range) )
                ## Step on the root is 1/10 of the space.
                ## This overrides the user parameters!
                cat("Computing step size for root value proposal from the data. \n")
                w_mu <- ( data.range[,2] - data.range[,1] ) / 10
                cat("Guessing magnitude of rates from the data. \n")
                if( is.list(phy[[1]]) ){ ## List of phylo.
                    ## Using a single core to compute BM model.
                    fit <- lapply(1:ncol(data), function(x) fitContinuous(phy = phy[[1]], dat=data[,x], model = "BM", ncores = 1) )
                } else{ ## Only one phylo.
                    ## Using a single core to compute BM model.
                    fit <- lapply(1:ncol(data), function(x) fitContinuous(phy = phy, dat=data[,x], model = "BM", ncores = 1) )
                }
                guess.rates <- sapply(fit, function(x) coef(x)[1])
                top.sd <- sqrt( ceiling( max(guess.rates) ) * 10 )
                rep.sd.regime <- rep(c(0,top.sd), times=p)
                par.sd <- matrix(data=rep.sd.regime, nrow=p, ncol=2, byrow=TRUE)
                prior_run <- makePrior(r=r, p=p, den.mu="norm", par.mu=data.range, par.sd=par.sd)
            }
        }
        
        ## #######################
        ## Block to generate start point.
        start_run <- start
        if( inherits(start, what="character") ){
            if(start == "prior_sample"){
                cat( "Taking sample from prior as starting point. \n" )
                start_run <- samplePrior(n=1, prior=prior_run)
            }
            if(start == "mle"){ ## Need to deal with the list of matrices here.
                cat( "Optimizing likelihood for the starting value of the MCMC.\n")
                if( is.list(phy[[1]]) ){ ## Is a list of phylogenies.
                    rr <- sample(1:length(phy), size=1)
                    cat( paste("Using phylogeny number ", rr, " to estimate the MLE.\n", sep="") )
                    phy.sample <- phy[[r]]
                    mle.fit <- mvBM(tree=phy.sample, data=data, model="BMM", method="rpf", echo=FALSE, diagnostic=FALSE)
                } else{
                    mle.fit <- mvBM(tree=phy, data=data, model="BMM", method="rpf", echo=FALSE, diagnostic=FALSE)
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
        
        out_mult <- multRegimeMCMC(X=data, phy=phy, start=start_run, prior=prior_run, gen=gen, post_seq=post_seq, v=v, w_sd=w_sd, w_mu=w_mu
                                 , prop=prop, dir=dir, outname=outname, IDlen=IDlen, regimes=regime.names, traits=trait.names
                                   , save.handle=save.handle, add.gen=NULL)
        return( out_mult )
        
    }
    
}
