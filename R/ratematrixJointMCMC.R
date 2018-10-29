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
##' @title Estimate the evolutionary rate matrix together with the regimes using Markov-chain Monte Carlo
##' @param data_BM a matrix with the data. Species names need to be provided as rownames (rownames(data) == phy$tip.label). Each column is a different trait. Names for the columns are used as trait labels. If labels are not provided, the function will use default labels.
##' @param data_Mk a named vector with the discrete data for the tips. Species names need to be provided as names for the vector. States are used to estimate the Markov model for the predictor regimes and perform stochastic mapping simulations.
##' @param phy a phylogeny of the class "phylo". Here the stochastic maps will be simulated together with the MCMC estimation. The regimes will follow the data provided as 'data_Mk'.
##' @param prior_BM the prior densities for the multivariate Brownian motion model. Must be one of "uniform", "uniform_scaled" (the default, see 'Details'), "empirical_mean", or the output of the "makePrior" function. See more information on 'makePrior' and in the examples below.
##' @param prior_Mk the prior density for the Markov model for the predictor regimes. Must be one of "uniform" or "exponential" (default is "uniform").
##' @param par_prior_Mk either a numeric vector with length 2 with the min and max for the uniform prior when 'prior_Mk = "uniform"' or a single value for the rate of the exponential distribution when 'prior_Mk = "exponential"'.
##' @param Mk_model the Markov model fitted to the predictor regimes and to make the stochastic map simulations. One of "SYM" for symmetrical rates (default value), "ARD" for all rates different, or "ER" for equal reates.
##' @param root_Mk the root type for the Mk model for the predictor regimes. Can be one of "madfitz" (default) or "equal".
##' @param smap_limit the maximum number of times that a stochastic map for a given branch can be attempted. If the simulation is not finished by this number of trials then the stochastic map is rejected. If set to '0' then there is no limit. The default value is 1e6.
##' @param start the starting state for the MCMC chain. Must be one of "prior_sample" (the default), "mle", or a sample from the prior generated with the "samplePrior" functions.
##' @param start_Q A matrix with the starting state for the Q matrix. Default is 'NULL' and the Q matrix is sampled from its prior distribution.
##' @param gen number of generations for the chain.
##' @param v value for the degrees of freedom parameter of the inverse-Wishart proposal distribution for the correlation matrix. Smaller values provide larger steps and larger values provide smaller steps. (Yes, it is counterintuitive.) This needs to be a single value applied to all regimes or a vector with the same length as the number of regimes.
##' @param w_sd the multiplying factor for the multiplier proposal on the vector of standard deviations. This can be a single value to be used for the sd of all traits for all regimes or a matrix with number of columns equal to the number of regimes and number of rows equal to the number of traits. If a matrix, then each element will be used to control the correspondent width of the standard deviation.
##' @param w_q the multiplying factor for the multiplier proposal on the transition matrix for the Markov model fitted to the predictor traits. Need to be a single value.
##' @param w_mu value for the width of the sliding window proposal for the vector of root values (phylogenetic mean). This can be a single value to be used for the root value of all traits or a vector of length equal to the number of traits. If a vector, then each element will be used as the width of the proposal distribution for each trait in the same order as the columns in 'data'. When 'prior="uniform_scaled"' (the default) this parameter is computed from the data.
##' @param prop a numeric vector of length 5 with the proposal frequencies for each parameter of the model. The vector need to sum to 1. Values are i this order: phylogenetic mean (prop[1]), standard deviations (prop[2]), correlation matrix (prop[3]), transition matrix (Mk) (prop[4]), and stochastic maps (prop[5]). Default value is 'c(0.050, 0.300, 0.300, 0.175, 0.175)'.
##' @param dir path of the directory to write the files. Has no default value (due to RCran policy). The path can be provided both as relative or absolute. It should accept Linux, Mac and Windows path formats.
##' @param outname name for the MCMC chain (default is 'ratematrixMCMC'). Name will be used in all the files alongside a unique ID of numbers with length of 'IDlen'.
##' @param IDlen length of digits of the numeric identifier used to name output files (default is 5).
##' @param save.handle whether the handle for the MCMC should be saved to the directory in addition to the output files.
##' @return Function returns the 'handle' object and writes the posterior distribution and log as files in the directory (see 'dir'). The handle is a list with the details of the MCMC chain. It is composed by: *k* the number of traits; *p* the number of R regimes fitted to the tree; *ID* the unique identifier of the run; *dir* the directory where the posterior and log files were saved; *outname* the name for the chain; *trait.names* a vector with the label for the traits; *regime.names* a vector with the label for the rate regimes; *data* the data used in the analysis; *phy* a single phylogeny or the list of phylogenies; *prior* a list with the prior functions; *start* a list with the starting parameters for the chain; *gen* the number of generations for the chain; *mcmc.par* a list with the tunning parameters for the MCMC.
##' @author Daniel S. Caetano and Luke J. Harmon
##' @export
##' @importFrom mvMORPH mvBM
##' @importFrom corpcor decompose.cov
##' @importFrom ape is.ultrametric
##' @importFrom ape Ntip
##' @importFrom geiger fitContinuous
##' @importFrom stats coef rexp
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
ratematrixJointMCMC <- function(data_BM, data_Mk, phy, prior_BM="uniform_scaled", prior_Mk="uniform", par_prior_Mk=c(0, 100), Mk_model = "SYM", root_Mk = "madfitz", smap_limit=1e6, start="prior_sample", start_Q = NULL, gen, v=50, w_sd=0.2, w_q=0.2, w_mu=0.5, prop=c(0.05, 0.3, 0.3, 0.175, 0.175), dir=NULL, outname="ratematrixJointMCMC", IDlen=5, save.handle=TRUE){

    ## #######################
    ## Block to check arguments, give warnings and etc.

    ## Quickly check if a directory was provide. If not return an error.
    if( is.null(dir) ) stop('Need to provide a path to write MCMC samples. Use dir="." to write files to current directory.')
    if( !inherits(dir, what="character") ) stop("Value for argument 'dir' need to be a character. See help page.")

    ## Corrects the data if necessary.
    if( class(data_BM) == "data.frame" ) data_BM <- as.matrix( data_BM )

    cat("\n")

    ## Inform that default options are being used:
    if( inherits(prior_BM, what="character") && prior_BM == "empirical_mean" ) cat("Using old default prior. \n")
    if( inherits(prior_BM, what="character") && prior_BM == "uniform_scaled" ) cat("Using new default prior. \n")
    if( inherits(start, what="character") && start == "prior_sample" ) cat("Using default starting point. \n")
    if( v == 50 && w_sd == 0.2 && w_q == 0.2 && w_mu == 0.5 && prop[1] == 0.05 && prop[2] == 0.3 && prop[3] == 0.3 && prop[4] == 0.175) cat("Using default proposal settings. \n")

    ## Check if the proposal frequency vector is good:
    if( length(prop) != 5 ) stop(" Vector of proposal frequencies 'prop' need to have length 5. ")
    prop <- prop / sum(prop)
    
    ## Check if provided prior has the correct dimension:
    if( inherits(prior_BM, what="ratematrix_prior_function") ){
        if( !prior_BM$pars$r == ncol(data_BM) ) stop("Wrong number of traits specified on the prior.")
    }
    
    ## Check formats for 'w_mu'. Need to delay check for 'w_sd' until we get the number of regimes.
    if( length( w_mu ) > 1 ){
        if( !length(w_mu) == ncol(data_BM) ) stop("Length of 'w_mu' need to be 1 or equal to number of traits.")
    } else{
        w_mu <- rep(w_mu, times=ncol(data_BM))
    }

    ## Check format for 'w_q':
    if( length( w_q ) > 1 ) stop(" Scaling factor for Q need to be a single numeric value. ")

    ## Check if the phylogeny is of type 'phylo' and if it is a single phylo.
    if( is.list(phy[[1]]) ) stop("The joint MCMC for the ratematrix and regime only works with a single phylogeny.")  
    ## Check if the tree is ultrametric, also rescale the tree if needed.
    if( !is.ultrametric(phy) ) warning("Phylogenetic tree is not ultrametric. Continuing analysis. Please check 'details'.")
    
    ## Check the data for the Mk model:
    equaln <- Ntip( phy ) == length( data_Mk )
    if( !equaln ) stop("Number of species in 'data_Mk' data not equal to tree.")
    same.spp <- sum(phy$tip.label %in% names( data_Mk), na.rm=TRUE) == length( data_Mk )
    if( !same.spp ) stop("Names of 'data_Mk' do not match the tip labels of the tree.")
    mm1 <- match(phy$tip.label, names(data_Mk))
    data_Mk <- data_Mk[mm1]

    ## Check if data for the mvBM model.
    equaln <- Ntip( phy ) == nrow( data_BM )
    if( !equaln ) stop("Number of species in data_BM not equal to tree.")
    same.spp <- sum(phy$tip.label %in% rownames( data_BM ), na.rm=TRUE) == nrow(data_BM)
    if( !same.spp ) stop("Rownames of data_BM does not match the tip labels of the tree.")
    ## Reorder the data_BM to be the same order as the tip labels.
    mm2 <- match(phy$tip.label, rownames(data_BM))
    data_BM <- data_BM[mm2,]

    ## Check the number of regimes in the Mk data with the prior:
    if( inherits(prior_BM, what="ratematrix_prior_function") ){
        if( !prior_BM$pars$p == length(unique(data_Mk)) ) stop("Number of regimes specified on prior does not match the states in the 'data_Mk' . See command: 'table(data_Mk)' ")        
    }

    n_regimes <- length(unique(data_Mk))
    if( is.matrix( w_sd ) ){
        if( !ncol(w_sd) == n_regimes ) stop("ncol(w_sd) need to be equal to the number of regimes.")
        if( !nrow(w_sd) == ncol(data_BM) ) stop("nrow(w_sd) need to be equal to the number of traits.")
    } else{
        if( n_regimes > 1 & !length( w_sd ) == 1 ) stop(" 'w_sd' need to be a single numeric value or a matrix with ncol equal to the number of regimes and nrow equal to the number of traits.")
        if( length( w_sd ) == 1 ){
            w_sd <- rep(w_sd, times=ncol(data_BM))
        }
        if( !length( w_sd ) == ncol(data_BM) ) stop("length of 'w_sd' need to be equal to the number of traits.")
        if( n_regimes > 1 ){ ## w_sd needs to be a matrix!
            temp_mat <- matrix(nrow=length(w_sd), ncol=n_regimes)
            for(i in 1:n_regimes) temp_mat[,i] <- w_sd
            w_sd <- temp_mat
        }
    }

    ## Make a quick check down the road to see if the prior is working.
    if( !inherits(prior_BM, what="ratematrix_prior_function") ){
        if( inherits(prior_BM, what="character") ){
            if( !prior_BM %in% c("uniform","empirical_mean", "uniform_scaled") ) stop("prior option was not recognized. Check details for the 'prior' argument.")
        } else{
            if( inherits(prior_BM, what="list") ){
                warning("MCMC chain using custom prior as speficied by argument 'prior'.")
            } else{
                stop("Custom prior need to be a list. Check function details.")
            }
        }
    }

    ## Check if the choice for the prior on the Mk model is correct:
    prior_Mk <- match.arg(prior_Mk, choices=c("uniform","exponential"), several.ok=FALSE)
    if( prior_Mk == "uniform" && length(par_prior_Mk) != 2) stop("Prior on the Q matrix is set to uniform. Parameter 'par_prior_Mk' needs to be a numeric vector with min and max for the uniform distribution.")
    if( prior_Mk == "exponential" && length(par_prior_Mk) != 1) stop("Prior on the Q matrix is set to uniform. Parameter 'par_prior_Mk' needs to be a single numeric value with the exponential rate.")
    ## Check the model type for the Mk model:
    Mk_model <- match.arg(Mk_model, choices=c("ER","ARD", "SYM"), several.ok=FALSE)
    root_Mk <- match.arg(root_Mk, choices=c("madfitz","equal"), several.ok=FALSE)

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
    if( is.null( colnames(data_BM) ) ){
        trait.names <- paste("trait_", 1:ncol(data_BM), sep="")
    } else{
        trait.names <- colnames(data_BM)
    }

    ## #######################
    ## Block to set the analysis.
    p <- length( unique( data_Mk ) )
    r <- ncol( data_BM)

    ## #######################
    ## Block to generate priors the multivariate BM model.
    ## This step does not apply for the Mk model.
    prior_run <- prior_BM
    if( inherits(prior_BM, what="character") ){
        if(prior_BM == "uniform"){
            max.data <- apply(data_BM, 2, max)
            min.data <- apply(data_BM, 2, min)
            par.mu <- cbind( min.data, max.data )
            rep.sd.regime <- rep(c(0,sqrt(10)), times=p)
            par.sd <- matrix(data=rep.sd.regime, nrow=p, ncol=2, byrow=TRUE)
            prior_run <- makePrior(r=r, p=p, par.mu=par.mu, par.sd=par.sd )
        }
        if(prior_BM == "empirical_mean"){
            mn <- colMeans(data_BM)
            ssd <- apply(data_BM, 2, stats::sd)
            par.mu <- as.matrix( cbind(mn, ssd) )
            rep.sd.regime <- rep(c(0,sqrt(10)), times=p)
            par.sd <- matrix(data=rep.sd.regime, nrow=p, ncol=2, byrow=TRUE)
            prior_run <- makePrior(r=r, p=p, den.mu="norm", par.mu=par.mu, par.sd=par.sd)
        }
        if(prior_BM == "uniform_scaled"){
            data.range <- t( apply(data_BM, 2, range) )
            ## Step on the root is 1/10 of the space.
            ## This overrides the user parameters!
            cat("Computing step size for root value proposal from the data. \n")
            w_mu <- ( data.range[,2] - data.range[,1] ) / 10
            cat("Guessing magnitude of rates from the data. \n")
            fit <- lapply(1:ncol(data_BM), function(x) fitContinuous(phy = phy, dat=data_BM[,x], model = "BM") )
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
            mle.fit <- mvBM(tree=phy, data=data_BM, model="BMM", method="rpf", echo=FALSE, diagnostic=FALSE)
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

    ## Gereate starting point and regime names for the Mk model.
    ## Keep the order of the regime states on the Q matrix. We will need it later.
    regime.names <- unique( data_Mk )

    if( is.null(start_Q) ){ ## Sample from the prior.
    if( Mk_model == "ER" ){
        if( prior_Mk == "uniform" ){
            par_Q <- runif(n=1, min = par_prior_Mk[1], max = par_prior_Mk[2])
            Q <- matrix(data=par_Q, nrow=p, ncol=p)
            diag(Q) <- -1 * ( rowSums(Q) - par_Q )
        } else{
            par_Q <- rexp(n=1, rate = par_prior_Mk)
            Q <- matrix(data=par_Q, nrow=p, ncol=p)
            diag(Q) <- -1 * ( rowSums(Q) - par_Q )
        }
    }
    if( Mk_model == "ARD" ){
        if( prior_Mk == "uniform" ){
            Qsize <- (p * p) - p
            par_Q <- runif(n=Qsize, min = par_prior_Mk[1], max = par_prior_Mk[2])
            Q <- matrix(data=0, nrow=p, ncol=p)
            count <- 1
            for( i in 1:p){
                for(j in 1:p){
                    if( i == j ) next
                    Q[i,j] <- par_Q[count]
                    count <- count + 1
                }
            }
            diag(Q) <- -1 * rowSums(Q)
        } else{
            Qsize <- (p * p) - p
            par_Q <- rexp(n=Qsize, rate = par_prior_Mk)
            Q <- matrix(data=0, nrow=p, ncol=p)
            count <- 1
            for( i in 1:p){
                for(j in 1:p){
                    if( i == j ) next
                    Q[i,j] <- par_Q[count]
                    count <- count + 1
                }
            }
            diag(Q) <- -1 * rowSums(Q)
        }
    }
    if( Mk_model == "SYM" ){
        if( prior_Mk == "uniform" ){
            Qsize <- ((p * p) - p) / 2
            par_Q <- runif(n=Qsize, min = par_prior_Mk[1], max = par_prior_Mk[2])
            Q <- matrix(data=0, nrow=p, ncol=p)
            count <- 1
            for( i in 1:p){
                for(j in 1:p){
                    if( i >= j ) next
                    Q[i,j] <- par_Q[count]
                    Q[j,i] <- par_Q[count]
                    count <- count + 1
                }
            }
            diag(Q) <- -1 * rowSums(Q)
        } else{
            Qsize <- ((p * p) - p) / 2
            par_Q <- rexp(n=Qsize, rate = par_prior_Mk)
            Q <- matrix(data=0, nrow=p, ncol=p)
            count <- 1
            for( i in 1:p){
                for(j in 1:p){
                    if( i >= j ) next
                    Q[i,j] <- par_Q[count]
                    Q[j,i] <- par_Q[count]
                    count <- count + 1
                }
            }
            diag(Q) <- -1 * rowSums(Q)
        }
    }
    } else { ## Use the provided Q matrix.
        Q <- start_Q
        if( is.null( colnames(Q) ) ) stop( "Provided Q matrix need to have colnames equal to the states." )
        if( all(colnames(Q) %in% data_Mk) )
        regime.names <- colnames(Q)
        if( !is.matrix(Q) ) stop( "Provided Q for the starting state need to be a matrix" )
        if( ncol(Q) != nrow(Q) | ncol(Q) != p ) stop( "Wrong number of dimensions for Q matrix." )
    }

    ## Make a stochastic map draw for the starting state of the mapped_edge matrix
    ## Need to reformat the 'data_Mk' object first to be able to work with this.
    ## The format is a matrix with present and absent data.
    matrix_Mk <- makeDataTips(X = data_Mk, states = regime.names)
    
    prun.phy <- reorder.phylo(x = phy, order = "postorder")
    edge_mat <- prun.phy$edge
    root_type <- as.numeric( switch(root_Mk, "madfitz" = 1, "equal" = 0) )

    ## This is the starting stochastic map.
    ## We need at least one to start. Let's try 100 times. If we cannot get a stochastic map after trying 100 times, then return and error.
    ntimes <- 1
    while( TRUE ){
    mapped.edge <- makeSimmapMappedEdge(n_nodes=Nnode(prun.phy), n_tips=Ntip(prun.phy), n_states=p
                                      , edge_len=prun.phy$edge.length, sims_limit=smap_limit
                                      , edge_mat=prun.phy$edge, parents=unique( prun.phy$edge[,1] )
                                      , X=matrix_Mk, Q=Q, root_node=(Ntip(prun.phy)+1), root_type=root_type)
    if( sum( mapped.edge ) > prun.phy$edge.length[1] ){
        break
    }
    ntimes <- ntimes + 1
    if( ntimes >= 101 ){
        stop("Starting Q matrix failed to generate a stochastic maps within the 'smap_limit'. Please read details and consider increasing value of 'smap_limit' or changing the starting Q matrix.")
        }
    }
        
    out_mult <- multRegimeJointMCMC(X_BM=data_BM, X_Mk=matrix_Mk, phy=prun.phy
                                  , start=start_run, smap_limit=smap_limit
                                  , prior=prior_run, start_Q = Q, root_Mk = root_type
                                  , start_mapped.edge = mapped.edge, prior_Mk = prior_Mk
                                  , par_prior_Mk = par_prior_Mk, Mk_model = Mk_model
                                  , gen=gen, v=v, w_sd=w_sd, w_mu=w_mu, w_q=w_q
                                  , prop=prop, dir=dir, outname=outname, IDlen=IDlen
                                  , regimes=regime.names, traits=trait.names
                                  , save.handle=save.handle, add.gen=NULL)
    return( out_mult )
}
