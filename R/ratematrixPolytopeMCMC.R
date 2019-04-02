##' Function runs a MCMC chain to estimate the multivariate phenotypic space occupied by the ancestrals given the data observed at the tips of the phylogeny and the model of trait macroevolution. The trait evolution model incorporates separated rates of evolution for each trait as well as the evolutionary correlation among the traits. This method is an extension of the model described in Caetano and Harmon (2018). The multivariate phenotypic space is treated as a n-dimensional polytope. We use data augmentation to estimate the ancestral states following previous work by Lemey et al. (2010) and Bouckaert et al. (2012).
##'
##' The input data is a named list with a matrix of observations for each of the species on the tree. For each matrix, the rows are individual observations and the columns are each of the traits. The number of rows between matrices can vary. At least two observations are necessary for each of the species. The number of columns need to be the same and in the same order among matrices. The names of the list elements need to match the names of the species in the tree. Names for the columns of each matrix are optional.
##'
##' This function does not allow for a list of phylogenies as the input. Since the aim is to estimate the phenotypic space occupied by ancestors, it is unclear how to summarize information across distinct topologies.
##' 
##' Please check 'ratematrixMCMC' for more details. This function follows the same algorithm in 'ratematrixMCMC' except for the data augmentation approach used to incorportate the phenotypic space for the traits at the tips and nodes of the tree. Also note that this function does not estimate the root value for the tree. Here we use the restricted likelihood for the multivariate Brownian motion model.
##' 
##' @title Joint estimate the multidimentional phenotypic space for ancestrals and the evolutionary rate matrix
##' @param data a named list with multiple trait values per species. Names of list elements need to match tip labels. list of matrices with observations as rows and trait values as columns. See "Details" for more information.
##' @param phy a single phylogeny of the class "simmap" with the mapped regimes for two or more R regimes OR a phylogeny of the class "phylo" for a single regime. The number of evolutionary rate matrices fitted to the phylogeny is equal to the number of regimes in 'phy'. Regime names will also be used.
##' @param ancestral_priors the prior for the ancestral nodes. This can be "descendants" or "global". If "descendants" then the prior at each internal node comprise the morphospace occupied by all descendants of that node. If "global" then all internal nodes share the same prior distribution with a morphospace spanning the phenotypic space observed by the pool of all species in the tree.
##' @param prior the prior densities for the MCMC. Must be one of "uniform", "uniform_scaled" (the default, see 'Details'), or the output of the "makePrior" function. See more information on 'makePrior' and in the examples below.
##' @param start the starting state for the MCMC chain. Must be one of "prior_sample" (the default) or a sample from the prior generated with the "samplePrior" functions.
##' @param gen number of generations for the chain.
##' @param burn percentage to discard as burn in. Set 'burn' in to 0.0 to start writing from the first generation.
##' @param thin number of generations to be skipped as thinning. Set 'thin' to 0 to write every generation to the file.
##' @param v value for the degrees of freedom parameter of the inverse-Wishart proposal distribution for the correlation matrix. Smaller values provide larger steps and larger values provide smaller steps. (Yes, it is counterintuitive.) This needs to be a single value applied to all regimes or a vector with the same length as the number of regimes.
##' @param w_sd the multiplying factor for the multiplier proposal on the vector of standard deviations. This can be a single value to be used for the sd of all traits for all regimes or a matrix with number of columns equal to the number of regimes and number of rows equal to the number of traits. If a matrix, then each element will be used to control the correspondent width of the standard deviation.
##' @param prop a numeric vector of length 2 with the proposal frequencies for the vector of standard deviations (prop[1]) and the correlation matrix (prop[2]). Default value is 'c(0.5, 0.5)'.
##' @param n_nodes_prop the number of nodes (among tip nodes and internal nodes) that are going to be updated during each step of the MCMC. Min of 1 and max should be no more than the number of tips in the tree.
##' @param dir path of the directory to write the files. Has no default value (due to RCran policy). The path can be provided both as relative or absolute. It should accept Linux, Mac and Windows path formats.
##' @param outname name for the MCMC chain (default is 'ratematrixMCMC'). Name will be used in all the files alongside a unique ID of numbers with length of 'IDlen'.
##' @param IDlen length of digits of the numeric identifier used to name output files (default is 5).
##' @param save.handle whether the handle for the MCMC should be saved to the directory in addition to the output files.
##' @return Function returns the 'handle' object and writes the posterior distribution and log as files in the directory (see 'dir'). The handle is a list with the details of the MCMC chain. It is composed by: *k* the number of traits; *p* the number of R regimes fitted to the tree; *ID* the unique identifier of the run; *dir* the directory where the posterior and log files were saved; *outname* the name for the chain; *trait.names* a vector with the label for the traits; *regime.names* a vector with the label for the rate regimes; *data* the data used in the analysis; *phy* a single phylogeny or the list of phylogenies; *prior* a list with the prior functions; *start* a list with the starting parameters for the chain; *gen* the number of generations for the chain; *mcmc.par* a list with the tunning parameters for the MCMC.
##' @author Daniel S. Caetano and Luke J. Harmon
##' @references Revell, L. J., and L. J. Harmon. 2008. Testing quantitative genetic hypotheses about the evolutionary rate matrix for continuous characters. Evolutionary Ecology Research 10:311.
##' @references Revell, L. J., and D. C. Collar. 2009. Phylogenetic Analysis of the Evolutionary Correlation Using Likelihood. Evolution 63:1090–1100.
##' @references Lemey, P., A. Rambaut, J. J. Welch, and M. A. Suchard. 2010. Phylogeography Takes a Relaxed Random Walk in Continuous Space and Time. Mol Biol Evol 27:1877–1885.
##' @references Bouckaert, R., P. Lemey, M. Dunn, S. J. Greenhill, A. V. Alekseyenko, A. J. Drummond, R. D. Gray, M. A. Suchard, and Q. D. Atkinson. 2012. Mapping the Origins and Expansion of the Indo-European Language Family. Science 337:957–960.
##' @references Caetano, D. S., and L. J. Harmon. 2017. ratematrix: An R package for studying evolutionary integration among several traits on phylogenetic trees. Methods in Ecology and Evolution 8:1920–1927.
##' @references Caetano, D. S., and L. J. Harmon. 2018. Estimating Correlated Rates of Trait Evolution with Uncertainty. Systematic Biology, doi: 10.1093/sysbio/syy067.
##' @export
##' @importFrom mvMORPH mvBM
##' @importFrom corpcor decompose.cov
##' @importFrom ape is.ultrametric
##' @importFrom ape Ntip
##' @importFrom geiger fitContinuous
##' @importFrom geiger tips
##' @importFrom stats coef
##' @examples
##' \donttest{
##' ## Need to add examples.
##' }
ratematrixPolytopeMCMC <- function(data, phy, ancestral_priors="descendants", prior="uniform_scaled", start="prior_sample", gen, burn = 0.25, thin = 100, v=50, w_sd=0.2, prop=c(0.5, 0.5), n_nodes_prop=1, dir=NULL, outname="ratematrixPolyMCMC", IDlen=5, save.handle=TRUE){

    ## #######################
    ## Block to check arguments, give warnings and etc.

    ## Check burn and thin and create the vector of generations.
    ## Note here that the first generation is gen 0.
    post_seq <- seq(from = gen * burn, to = (gen-1), by = thin)
    post_seq <- post_seq + 1 ## Bounce the generation vector forward to start from 1.
    
    ## Quickly check if a directory was provided. If not return an error.
    if( is.null(dir) ) stop('Need to provide a path to write MCMC samples. Use dir="." to write files to current directory.')
    if( !inherits(dir, what="character") ) stop("Value for argument 'dir' need to be a character. See help page.")

    ## Data here needs to be a named list.
    if( !inherits(data, what="list") ) stop("Data needs to be a named list. See help page.")
    ## Functions only allow for a single phylogeny:
    if( is.list(phy[[1]]) ) stop("Only a single phylogeny is allowed.")
    
    if( length(data) != Ntip(phy) ) stop("Data needs to be a named list with the same length of the tips of the phylogeny. See help page.")
    name.list <- names(data)
    if( is.null(name.list) ) stop("Data needs to be a named list. See help page.")
    if( !all( name.list %in% phy$tip.label ) ) stop("List names need to match the names of the tip labels of the phylogeny.")
    
    ## Check if some of the species have a single observation.
    ## This can be either a matrix with one row or a vector.
    data.vec <- !sapply(data, is.matrix)
    if( any(data.vec) ){
        ## At least one of the elements of the data is a vector (not a matrix).
        for( i in which(data.vec) ){
            data[[i]] <- matrix(data[[i]], ncol=length( data[[i]] ), nrow=1)
        }
    }
    single.obs <- sapply(data, function(x) nrow(x) == 1)
    if( any(single.obs) ){
        cat("Some species have a single observation. Using fixed value. \n")
        ## Single we compute the range to set up the matrix for the MCMC, need to replicate the single obs for the range function to work.
        for( i in which(single.obs) ){
            data[[i]] <- rbind( data[[i]], data[[i]] )
        }
    }

    ## Check if all matrices have the same number of columns.
    ncol.list <- sapply(data, ncol)
    if( length( unique( ncol.list ) ) > 1 ) stop("All matrices with trait observations need to have the same number of traits.")
    ## Check if all matrices have at least two rows. Unclear what to do with a single observation at the moment.
    nrow.list <- sapply(data, nrow)
    if( any( nrow.list < 2 ) ) stop("All matrices with trait observations need to have at least 2 rows.")
    
    ## Check which of the matrices are data.frame, transform those to matrices.
    ## Crude, but does the job.
    which.dataframe <- sapply(data, function(x) class(x) == "data.frame")
    for( i in 1:length(data) ){
        if( which.dataframe[[i]] ){
            data[[i]] <- as.matrix( data[[i]] )
        }
    }

    ## Make sure the list is in the correct order.
    order.ID <- match(phy$tip.label, name.list)
    ordered.list <- list()
    for( i in 1:length(data) ){
        ordered.list[[i]] <- data[[ order.ID[i] ]]
    }
    names( ordered.list ) <- name.list[ order.ID ]
    
    cat("\n")

    ## Check the characteristics of the prop vector:
    if( length( prop ) != 2 ) stop("prop vector needs length of 2.")
    prop <- prop / sum( prop ) ## Make sure it sums to 1.

    ## Check the prior options for the ancestral nodes:
    ancestral_priors <- match.arg(ancestral_priors, choices=c("descendants", "global"), several.ok=FALSE)
    cat(paste0("Using ", ancestral_priors, " option for prior at internal nodes. \n"))
    
    ## Inform that default options are being used:
    if( inherits(prior, what="character") && prior == "uniform_scaled" ) cat("Using new default prior. \n")
    if( inherits(start, what="character") && start == "prior_sample" ) cat("Using default starting point. \n")
    if( v == 50 && w_sd == 0.2 && prop[1] == 0.5 && prop[2] == 0.5 ) cat("Using default proposal settings. \n")

    ## Check if provided prior has the correct dimension:
    if( inherits(prior, what="ratematrix_prior_function") ){
        if( !prior$pars$r == ncol(data[[1]]) ) stop("Wrong number of traits specified on the prior.")
    }

    ## Check attributes of the phylogeny.
    ## Check if the tree is ultrametric, also rescale the tree if needed.
    if( !is.ultrametric(phy) ) warning("Phylogenetic tree is not ultrametric. Continuing analysis. Please check 'details'.")
    ## Check if the phylogeny is of 'simmap' class.
    if( !inherits(phy, what="simmap") ){
        cat('phy is not of class "simmap". Fitting a sigle rate regime to the tree. \n')
        no_phymap <- TRUE
    } else{
        no_phymap <- FALSE
    }

    ## Check the 'w_sd' parameter.
    ## Get the number of regimes and check the prior, if provided.
    if( no_phymap ){
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
    
    if( is.matrix( w_sd ) ){
        if( !ncol(w_sd) == n_regimes ) stop("ncol(w_sd) need to be equal to the number of regimes.")
        if( !nrow(w_sd) == ncol(data[[1]]) ) stop("nrow(w_sd) need to be equal to the number of traits.")
    } else{
        if( n_regimes > 1 & !length( w_sd ) == 1 ) stop(" 'w_sd' need to be a single numeric value or a matrix with ncol equal to the number of regimes and nrow equal to the number of traits.")
        if( length( w_sd ) == 1 ){
            w_sd <- rep(w_sd, times=ncol(data[[1]]))
        }
        if( !length( w_sd ) == ncol(data[[1]]) ) stop("length of 'w_sd' need to be equal to the number of traits.")
        if( n_regimes > 1 ){ ## w_sd needs to be a matrix!
            temp_mat <- matrix(nrow=length(w_sd), ncol=n_regimes)
            for(i in 1:n_regimes) temp_mat[,i] <- w_sd
            w_sd <- temp_mat
        }
    }

    ## Make a quick check down the road to see if the prior is working.
    if( !inherits(prior, what="ratematrix_prior_function") ){
        if( inherits(prior, what="character") ){
            if( !prior %in% c("uniform", "uniform_scaled") ) stop("prior option was not recognized. Check details for the 'prior' argument.")
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
        if( start != "prior_sample" ) stop("start state option was not recognized. Check details for the 'start' argument.")
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
    if( is.null( colnames(data[[1]]) ) ){
        trait.names <- paste("trait_", 1:ncol(data[[1]]), sep="")
    } else{
        trait.names <- colnames( data[[1]] )
    }
    
    ## First check if analysis will use regimes.
    if( !no_phymap ){
        if( is.null( colnames(phy$mapped.edge) ) ){
            regime.names <- paste("regime_", 1:ncol(phy$mapped.edge), sep="")
        } else{
            ## In the case we do not have names.
            regime.names <- colnames(phy$mapped.edge)
        }
    } else{
        ## In the case of a single regime analysis:
        regime.names <- NA
    }

    ## #######################
    ## Run the analyses with multiple regimes.
    ## This block is also making the estimates for the case of a single rate model.
    ## So pay attention to the differences between these cases.
    ## #######################

    ## Number of regimes.
    if( no_phymap ){
        p <- 1
    } else{
        p <- ncol( phy$mapped.edge ) ## Multiple regimes.
    }
    k <- ncol( data[[1]] ) ## Number of traits.

    ## #######################
    ## Block to generate priors.
    prior_run <- prior
    if( inherits(prior, what="character") ){
        if(prior == "uniform"){
            data.range <- t( apply(data[[1]], 2, range) )
            rep.sd.regime <- rep(c(0,sqrt(10)), times=p)
            par.sd <- matrix(data=rep.sd.regime, nrow=p, ncol=2, byrow=TRUE)
            prior_run <- makePrior(r=k, p=p, den.mu = "unif", par.mu=data.range, par.sd=par.sd )
        }
        if(prior == "uniform_scaled"){
            data.range <- t( apply(data[[1]], 2, range) )
            cat("Guessing magnitude of rates from the data. \n")
            ## Using a single random observation for each of the species.
            simple.data <- matrix(ncol = ncol(data[[1]]), nrow = Ntip(phy) )
            for( i in 1:Ntip(phy) ){
                get.mat <- data[[i]]
                get.mat.row <- sample(1:nrow(get.mat), size = 1)
                simple.data[i,] <- get.mat[get.mat.row,]
            }
            rownames( simple.data ) <- names( ordered.list )
            fit <- lapply(1:ncol(simple.data), function(x) fitContinuous(phy = phy, dat=simple.data[,x], model = "BM", ncores = 1) )
            guess.rates <- sapply(fit, function(x) coef(x)[1])
            top.sd <- sqrt( ceiling( max(guess.rates) ) * 10 )
            rep.sd.regime <- rep(c(0,top.sd), times=p)
            par.sd <- matrix(data=rep.sd.regime, nrow=p, ncol=2, byrow=TRUE)
            prior_run <- makePrior(r=k, p=p, den.mu="unif", par.mu=data.range, par.sd=par.sd)
        }
    }
    
    ## #######################
    ## Block to generate start point.
    ## #######################
    
    start_run <- start
    if( inherits(start, what="character") ){
        if(start == "prior_sample"){
            cat( "Taking sample from prior as starting point. \n" )
            start_run <- samplePrior(n=1, prior=prior_run)
        } else{
            stop("Only prior sample type allowed is 'prior_sample'.")
        }
    }

    ## #######################
    ## Block to make the estimate of the model.
    ## #######################

    ## The steps below are similar to the steps of the function 'multRegimeMCMC'.
    ## Decided not to call yet another function to do this.

    ## Check the value of v.
    if( p > 1 ){
        if( length( v ) > 1 ){
            if( !length( v ) == p ) stop( "Length of v need to be 1 or equal to the number of regimes." )
        }
        if( length( v ) == 1 ){
            v <- rep(v, times=p)
        }
    } else{
        ## Using a single regime. Keep only the first element of the vector.
        if( length( v ) > 1 ){
            v <- v[1]
        }
    }    

    ## Take care about the format of the data matrix.
    ## Need to transpose the data matrix.
    ## X <- t(X)

    ## The proposal for the standard deviation.
    prob_sample_sd <- prop[1]

    ## Save the list with the MCMC parameters.
    mcmc.par <- list()
    mcmc.par$v <- v
    mcmc.par$w_sd <- w_sd
    mcmc.par$w_mu <- NA ## Not using this parameter here.
    mcmc.par$prob_sample_root <- NA
    mcmc.par$prob_sample_sd <- prob_sample_sd
    mcmc.par$prop <- prop

    ## I dropped the option to continue the MCMC. So here we always start a new one.
    ## if( !is.null(continue) ){
    ##     ## Use the provided ID number for the run.
    ##     mcmc_file_name <- file.path(dir, paste(outname,".", ID, ".mcmc",sep=""))
    ##     log_file_name <- file.path(dir, paste(outname,".", ID, ".log",sep=""))
    ##     poly_file_name <- file.path(dir, paste(outname,".", ID, ".traitspace",sep=""))
    ##     write_header <- 0
    ##     gen <- add.gen
    ## } else{
    ## Generate identifier and name for the files:
    new.ID <- paste( sample(x=1:9, size=IDlen, replace=TRUE), collapse="")
    mcmc_file_name <- file.path(dir, paste(outname, ".", new.ID, ".mcmc", sep=""))
    log_file_name <- file.path(dir, paste(outname, ".", new.ID, ".log", sep=""))
    poly_file_name <- file.path(dir, paste(outname, ".", new.ID, ".traitspace", sep=""))
    write_header <- 1
    ## }

    ## Compute important info from the phylogeny:
    ## Order for traversal.
    ord.id <- reorder.phylo(phy, order="postorder", index.only = TRUE)
    if( no_phymap ){
        ## If the model is a single rate model, then we need to pass the branch lengths of the phylogeny as the mapped.edge.
        ## We are doing this because I merged the functions that deal with multiple regimes or a single regime.
        ## Note that, although this is a matrix, we will only use a single regime (first column) of it.
        mapped.edge <- cbind( phy$edge.length, phy$edge.length )
        mapped.edge <- mapped.edge[ord.id,] ## Need to reorganize the edges.
    } else{
        mapped.edge <- phy$mapped.edge[ord.id,] ## The regimes.
    }
    ## Need to take care how to match the regimes and the R matrices.
    anc <- phy$edge[ord.id,1] ## Ancestral edges.
    des <- phy$edge[ord.id,2] ## Descendent edges.
    nodes <- unique(anc) ## The internal nodes we will traverse.

    ## Set the types for each of the nodes that are going to be visited.
    node.to.tip <- which( tabulate( anc[which(des <= length(phy$tip.label))] ) == 2 )
    node.to.node <- which( tabulate( anc[which(des > length(phy$tip.label))] ) == 2 )
    node.to.tip.node <- unique( anc )[!unique( anc ) %in% c(node.to.node, node.to.tip)]
    ## 1) nodes to tips: nodes that lead only to tips, 2) nodes to nodes: nodes that lead only to nodes, 3) nodes to tips and nodes: nodes that lead to both nodes and tips.
    names(anc) <- rep(1, times=length(anc))
    names(anc)[which(anc %in% node.to.node)] <- 2
    names(anc)[which(anc %in% node.to.tip.node)] <- 3
    names_anc <- as.numeric( names(anc) )

    ## Before running need to exclude the generations already done if continuing.
    ## Also add the option to do additional generations.
    ## We dropped the option to continue the MCMC.
    ## if( is.null(continue) ){
    cat( paste("Start MCMC run ", outname, ".", new.ID, " with ", gen, " generations.\n", sep="") )
    ## } else{
    ##     if( continue == "continue" ){
    ##         cat( paste("Continue previous MCMC run ", outname, ".", ID, " for ", add.gen, " generations for a total of ", gen, " generations.\n", sep="") )
    ##         gen <- add.gen
    ##         new.ID <- ID
    ##     }
    ##     if( continue == "add.gen" ){
    ##         cat( paste("Adding ", add.gen, " generations to previous MCMC run ", outname, ".", ID, "\n", sep="") )
    ##         gen <- add.gen
    ##         new.ID <- ID
    ##     }
    ## }

    ## Save the handle object:
    if( save.handle ){
        out <- list(k = k, p = p, ID = new.ID, dir = dir, outname = outname
                  , trait.names = trait.names, regime.names = regime.names
                  , data = data, phy = phy, prior = prior_run, start = start_run
                  , gen = gen, mcmc.par = mcmc.par)
        class( out ) <- "ratematrix_poly_mcmc"
        saveRDS(out, file = file.path(dir, paste(outname,".",new.ID,".mcmc.handle.rds",sep="")) )
    }

    ## Set the objects holding the initial state for the chain.
    ## These are array and matrix classes. 'start_run' here is produced by the 'samplePrior' function.
    startR <- array(dim=c(k, k, p))
    startCorr <- array(dim=c(k, k, p))
    startvar <- matrix(nrow=k, ncol=p)
    if( no_phymap ){
        ## In the case of a single regime model the sample from the prior distribution will be different.
        startR.list <- rebuild.cov(r = cov2cor(start_run$matrix), v = start_run$sd^2)
        startR[,,1] <- startR.list
        startCorr[,,1] <- start_run$matrix
        startvar[,1] <- start_run$sd^2
    } else{
        startR.list <- lapply(1:p, function(x) rebuild.cov(r=cov2cor(start_run$matrix[[x]]), v=start_run$sd[[x]]^2) )
        for( i in 1:p ){
            startR[,,i] <- startR.list[[i]]
            startCorr[,,i] <- start_run$matrix[[i]]
            startvar[,i] <- start_run$sd[[i]]^2
        }
    }

    ## Get info from the prior object.
    den_sd <- prior_run$pars$den.sd
    par_sd <- prior_run$pars$par.sd

    if( prior_run$pars$unif.corr ){
        sigma.mat <- diag(nrow=k)
        sigma_array <- array(dim=c(k, k, p))
        for( i in 1:p){
            sigma_array[,,i] <- sigma.mat
        }
        nu <- rep(k+1, times=p)
    } else{
        if( length(prior_run$pars$Sigma) != p ) stop( "Length of Sigma need to be equal to number of regimes." )
        sigma_array <- sapply(prior_run$pars$Sigma, identity, simplify="array")
        nu <- prior_run$pars$nu
        if( length(nu) != p ) stop( "Length of nu need to be equal to number of regimes." )
    }

    ## Prepare the specific objects for the Polytope version of the MCMC:
    ## For the data at the tips we will use the min and maximum.
    X_poly <- matrix(nrow = Ntip(phy), ncol = k * 2)
    for( i in 1:Ntip(phy) ){
        range_poly <- vector( mode = "numeric" )
        for( j in 1:k ){
            range_poly <- c(range_poly, range( data[[i]][,j] ) )
        }
        X_poly[i,] <- range_poly
    }
    rownames( X_poly ) <- names( data )

    ## Set up the boundaries for the internal nodes.
    ## This case will depend on the type of prior used.
    if( ancestral_priors == "descendants" ){
        ## Range of the descendant lineages
        ## This will depend on visiting each node, selecting the descendant tips and computing the range.
        ## The anc_poly matrix is in the same order as the nodes.
        anc_poly <- matrix(nrow = length(nodes)-1, ncol = k*2)
        for( i in 1:nrow(anc_poly) ){
            tmp_des_mat <- X_poly[tips(phy = phy, node = nodes[i]), ]
            tmp_vec <- vector()
            for( j in 1:(k*2) ){
                if((j %% 2) < 1){ ## Test for zero can create problems.
                    ## Take the max
                    tmp_vec <- c(tmp_vec, max(tmp_des_mat[,j]))
                }else{
                    ## Take the min
                    tmp_vec <- c(tmp_vec, min(tmp_des_mat[,j]))
                }
            }
            anc_poly[i,] <- tmp_vec
        }

    } else{
        ## Range of all the tips pooled.
        anc_global <- vector()
        for( i in 1:(k*2) ){
            if((i %% 2) < 1){ ## Test for zero can create problems.
                ## Take the max
                anc_global <- c(anc_global, max(X_poly[,i]))
            }else{
                ## Take the min
                anc_global <- c(anc_global, min(X_poly[,i]))
            }
        }
        ## All internal nodes share the same bounds.
        ## The matrix has nrows equal to the number of nodes-1.
        anc_poly <- matrix(anc_global, nrow = length(nodes)-1, ncol = k*2, byrow = TRUE)
    }

    ## In the case of single regime analyses, some of the objects might be passed as a vector and need to be a matrix:
    ## sd need to be a matrix
    ## W_sd need to be a matrix
    ## par_prior_sd need to be a matrix
    ## Most simple way out is to duplicate the vectors to create a matrix and ignore one of the element insite the MCMC function.
    sd_mat_par <- sqrt(startvar)
    if( p == 1 ){
        ## Create a matrix out of some of the parameters to be able to run the model.
        sd_mat_par <- cbind( sd_mat_par, sd_mat_par)
        w_sd <- cbind(w_sd, w_sd)
        par_sd <- rbind(par_sd, par_sd)
    }    
    
    ## Pass the parameters and run the C++ MCMC function:
    runRatematrixPolytopeMCMC(X_poly=X_poly, anc_poly=anc_poly, n_input_move=n_nodes_prop
                            , k=k, p=p, nodes=nodes, des=des, anc=anc
                            , names_anc=names_anc, mapped_edge=mapped.edge, R=startR
                            , sd=sd_mat_par, Rcorr=startCorr, w_sd=w_sd
                            , par_prior_sd=par_sd, den_sd=den_sd, nu=nu
                            , sigma=sigma_array, v=v, log_file=log_file_name
                            , mcmc_file=mcmc_file_name, poly_file=poly_file_name
                            , prob_sample_sd=prob_sample_sd, gen=gen, post_seq=post_seq
                            , write_header=write_header)

    cat( paste("Finished MCMC run ", outname, ".", new.ID, "\n", sep="") )

    out <- list(k = k, p = p, ID = new.ID, dir = dir, outname = outname
              , trait.names = trait.names, regime.names = regime.names, data = data
              , phy = phy, prior = prior_run, start = start_run, gen = gen
              , mcmc.par=mcmc.par)
    class( out ) <- "ratematrix_poly_mcmc"
    return( out )
    
}
