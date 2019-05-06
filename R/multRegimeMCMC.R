##' Runs a Markov chain Monte Carlo (MCMC) chain to estimate the posterior distribution of two or more
##'    evolutionary rate matrices (R) fitted to the phylogeny.
##'
##' MCMC using the inverse-Wishart as a proposal distribution for the covariance matrix and a simple sliding
##'    window for the phylogenetic mean in the random walk Metropolis-Hastings algorithm. At the moment the
##'    function only applies the 'rpf' method for calculation of the log likelihood implemented in the
##'    package 'mvMORPH'. Future versions should offer different log likelihood methods for the user.
##' This version is using the pruning algoritm to c
##' @title MCMC for two or more evolutionary rate matrices.
##' @param X matrix. A matrix with the data. 'rownames(X) == phy$tip.label'.
##' @param phy simmap phylo. A phylogeny of the class "simmap" from the package 'phytools'. Function uses the location information for a number of traits equal to the number of fitted matrices.
##' @param start list. Element [[1]] is the starting value for the phylogenetic mean and element [[2]] is the starting value for the R matrices. Element [[2]] is also a list with length equal to the number of matrices to be fitted to the data.
##' @param prior list. Produced by the output of the function 'make.prior.barnard' or 'make.prior.diwish'. First element of the list [[1]] is a prior function for the log density of the phylogenetic mean and the second element [[2]] is a prior function for the evolutionary rate matrix (R). The prior can be shared among the rate matrices or be set a different prior for each matrix. At the moment the function only produces a shared prior among the fitted matrices. Future versions will implement independent priors for each of the fitted matrices.
##' @param gen numeric. Number of generations of the MCMC.
##' @param post_seq parameter
##' @param v numeric. Degrees of freedom parameter for the inverse-Wishart proposal distribution for the evolutionary rate matrix.
##' @param w_sd numeric. Width of the uniform sliding window proposal for the vector of standard deviations.
##' @param w_mu numeric. Width of the uniform sliding window proposal for the vector of phylogenetic means.
##' @param prop vector. The proposal frequencies. Vector with two elements (each between 0 and 1). First is the probability that the phylogenetic mean will be sampled for a proposal step at each genetarion, second is the probability that the evolutionary rate matrix will be updated instead. First the function sample whether the root value or the matrix should be updated. If the matrix is selected for an update, then one of the matrices fitted to the phylogeny is selected to be updated at random with the same probability.
##' @param dir string. Directory to write the files, absolute or relative path. If 'NULL' then output is written to the directory where R is running (see 'getwd()'). If a directory path is given, then function will test if the directory exists and use it. If directiory does not exists the function will try to create one.
##' @param outname string. Name pasted to the files. Name of the output files will start with 'outname'.
##' @param IDlen numeric. Set the length of the unique numeric identifier pasted to the names of all output files. This is set to prevent that multiple runs with the same 'outname' running in the same directory will be lost.Default value of 5 numbers, something between 5 and 10 numbers should be good enough. IDs are generated randomly using the function 'sample'.
##' @param regimes parameter
##' @param traits parameter
##' @param save.handle parameter
##' @param continue parameter
##' @param add.gen parameter
##' @param ID parameter
##' @return Fuction creates files with the MCMC chain. Each run of the MCMC will be identified by a unique identifier to facilitate identification and prevent the function to overwrite results when running more than one MCMC chain in the same directory. See argument 'IDlen'. The files in the directory are: 'outname.ID.loglik': the log likelihood for each generation, 'outname.ID.n.matrix': the evolutionary rate matrix n, one per line. Function will create one file for each R matrix fitted to the tree, 'outname.ID.root': the root value, one per line. \cr
##' \cr
##' Additionally it returns a list object with information from the analysis to be used by other functions. This list is refered as the 'out' parameter in those functions. The list is composed by: 'acc_ratio' numeric vector with 0 when proposal is rejected and non-zero when proposals are accepted. 1 indicates that root value was accepted, 2 and higher indicates that the first or subsequent matrices were updated; 'run_time' in seconds; 'k' the number of matrices fitted to the tree; 'p' the number of traits in the analysis; 'ID' the identifier of the run; 'dir' directory were output files were saved; 'outname' the name of the chain, appended to the names of the files; 'trait.names' A vector of names of the traits in the same order as the rows of the R matrix, can be used as the argument 'leg' for the plotting function 'make.grid.plot'; 'data' the original data for the tips; 'phy' the phylogeny; 'prior' the list of prior functions; 'start' the list of starting parameters for the MCMC run; 'gen' the number of generations of the MCMC.
##' @importFrom geiger treedata
##' @importFrom ape reorder.phylo
##' @importFrom corpcor rebuild.cov
##' @importFrom stats cov2cor
##' @noRd
multRegimeMCMC <- function(X, phy, start, prior, gen, post_seq, v, w_sd, w_mu, prop, dir, outname, IDlen, regimes, traits, save.handle, continue=NULL, add.gen=NULL, ID=NULL){

    ## Get the number of regimes.
    p <- length(start[[2]])
    
    ## Check the value of v.
    if( length( v ) > 1 ){
        if( !length( v ) == p ) stop( "Length of v need to be 1 or equal to the number of regimes." )
    }
    if( length( v ) == 1){
        v <- rep(v, times=p)
    }
    
    ## Need to transpose the data matrix.
    X <- t(X)
    k <- nrow(X) ## Number of traits.

    ## Translate the proportion to the objects:
    if( !length(prop) == 3 ) stop("Length of prop need to be equal to 3. Check help page.")
    if( !sum(prop) == 1 ) stop("Values of prop need to sum to 1.")
    prob_sample_root <- prop[1]
    ## Need to get the chance to sample sd rescaled.
    prob_sample_sd <- prop[2]/sum(prop[2:3])

    ## Save the list with the MCMC parameters.
    mcmc.par <- list()
    mcmc.par$v <- v
    mcmc.par$w_sd <- w_sd
    mcmc.par$w_mu <- w_mu
    mcmc.par$prob_sample_root <- prob_sample_root
    mcmc.par$prob_sample_sd <- prob_sample_sd
    mcmc.par$prop <- prop
    
    if( !is.null(continue) ){
        ## Use the provided ID number for the run.
        mcmc_file_name <- file.path(dir, paste(outname,".", ID, ".mcmc",sep=""))
        log_file_name <- file.path(dir, paste(outname,".", ID, ".log",sep=""))
        write_header <- 0
        gen <- add.gen
    } else{
        ## Generate identifier and name for the files:
        new.ID <- paste( sample(x=1:9, size=IDlen, replace=TRUE), collapse="")
        mcmc_file_name <- file.path(dir, paste(outname,".",new.ID,".mcmc",sep=""))
        log_file_name <- file.path(dir, paste(outname,".",new.ID,".log",sep=""))
        write_header <- 1
    }
    
    
    if( !is.list(phy[[1]]) ){ ## There is only one phylogeny.
        cat("MCMC chain using a single tree/regime configuration.\n")
        ord.id <- reorder.phylo(phy, order="postorder", index.only = TRUE) ## Order for traversal.
        mapped.edge <- phy$mapped.edge[ord.id,] ## The regimes.
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
    }

    if( is.list(phy[[1]]) ){ ## Phy is a list of phylogenies.
        cat("MCMC chain using multiple trees/regime configurations.\n")
        ord.id <- lapply(phy, function(x) reorder.phylo(x, order="postorder", index.only = TRUE))
        ntrees <- length(ord.id) ## Number of trees in the data.
        mapped.edge <- sapply(1:ntrees, function(x) phy[[x]]$mapped.edge[ord.id[[x]],], simplify="array") ## The regimes.
        ## Need to take care how to match the regimes and the R matrices.
        anc <- sapply(1:ntrees, function(x) phy[[x]]$edge[ord.id[[x]],1]) ## Ancestral edges.
        des <- sapply(1:ntrees, function(x) phy[[x]]$edge[ord.id[[x]],2]) ## Descendent edges.
        nodes <- sapply(1:ntrees, function(x) unique(anc[,x])) ## The internal nodes we will traverse.

        ## Set the types for each of the nodes that are going to be visited.
        node.to.tip <- sapply(1:ntrees, function(x) which( tabulate( anc[which(des[,x] <= length(phy[[x]]$tip.label)), x] ) == 2 ) )
        node.to.node <- sapply(1:ntrees, function(x) which( tabulate( anc[which(des[,x] > length(phy[[x]]$tip.label)), x] ) == 2 ) )
        node.to.tip.node <- sapply(1:ntrees, function(x) unique( anc[,x] )[!unique( anc[,x] ) %in% c(node.to.node[,x], node.to.tip[,x])] )
        ## 1) nodes to tips: nodes that lead only to tips, 2) nodes to nodes: nodes that lead only to nodes, 3) nodes to tips and nodes: nodes that lead to both nodes and tips.
        names_anc <- matrix(data=1, nrow=nrow(anc), ncol=ncol(anc))
        for( i in 1:ncol(anc) ){
            names_anc[which(anc[,i] %in% node.to.node[,i]), i] <- 2
            names_anc[which(anc[,i] %in% node.to.tip.node[,i]), i] <- 3            
        }
    }
    
    ## Before running need to exclude the generations already done if continuing.
    ## Also add the option to do additional generations.
    if( is.null(continue) ){
        cat( paste("Start MCMC run ", outname, ".", new.ID, " with ", gen, " generations.\n", sep="") )
    } else{
        if( continue == "continue" ){
            cat( paste("Continue previous MCMC run ", outname, ".", ID, " for ", add.gen, " generations for a total of ", gen, " generations.\n", sep="") )
            gen <- add.gen
            new.ID <- ID
        }
        if( continue == "add.gen" ){
            cat( paste("Adding ", add.gen, " generations to previous MCMC run ", outname, ".", ID, "\n", sep="") )
            gen <- add.gen
            new.ID <- ID
        }
    }

    ## Save the handle object:
    if( save.handle ){
        out <- list(k = k, p = p, ID = new.ID, dir = dir, outname = outname, trait.names = traits
                  , regime.names = regimes, data = t(X), phy = phy, prior = prior, start = start, gen = gen
                  , mcmc.par = mcmc.par, post_seq = post_seq)
        class( out ) <- "ratematrix_multi_mcmc"
        saveRDS(out, file = file.path(dir, paste(outname,".",new.ID,".mcmc.handle.rds",sep="")) )
    }

    ## Set the objects holding the initial state for the chain.
    ## These are array and matrix classes.
    ## start here is produced by the 'samplePrior' function.
    startR.list <- lapply(1:p, function(x) rebuild.cov(r=cov2cor(start$matrix[[x]]), v=start$sd[[x]]^2) )
    startR <- array(dim=c(k, k, p))
    startCorr <- array(dim=c(k, k, p))
    startvar <- matrix(nrow=k, ncol=p)
    for( i in 1:p){
        startR[,,i] <- startR.list[[i]]
        startCorr[,,i] <- start$matrix[[i]]
        startvar[,i] <- start$sd[[i]]^2
    }

    ## Get info from the prior object.
    den_mu <- prior$pars$den.mu
    par_mu <- prior$pars$par.mu
    den_sd <- prior$pars$den.sd
    par_sd <- prior$pars$par.sd

    if( prior$pars$unif.corr ){
        sigma.mat <- diag(nrow=k)
        sigma_array <- array(dim=c(k, k, p))
        for( i in 1:p){
            sigma_array[,,i] <- sigma.mat
        }
        nu <- rep(k+1, times=p)
    } else{
        if( length(prior$pars$Sigma) !=p ) stop( "Length of Sigma need to be equal to number of regimes." )
        sigma_array <- sapply(prior$pars$Sigma, identity, simplify="array")
        nu <- prior$pars$nu
        if( length(nu) !=p ) stop( "Length of nu need to be equal to number of regimes." )
    }

    ## Pass the arguments and start the MCMC.
    if( !is.list(phy[[1]]) ){ ## There is only one phylogeny.
        runRatematrixMCMC_C(X=X, k=k, p=p, nodes=nodes, des=des, anc=anc, names_anc=names_anc
                          , mapped_edge=mapped.edge, R=startR, mu=start$root, sd=sqrt(startvar), Rcorr=startCorr, w_mu=w_mu
                          , par_prior_mu=par_mu, den_mu=den_mu, w_sd=w_sd, par_prior_sd=par_sd, den_sd=den_sd
                          , nu=nu, sigma=sigma_array, v=v, log_file=log_file_name, mcmc_file=mcmc_file_name
                          , prob_sample_root = prob_sample_root, prob_sample_sd = prob_sample_sd, gen = gen
                          , post_seq = post_seq, write_header = write_header)
    }
    if( is.list(phy[[1]]) ){ ## Phy is a list of phylogenies.
        runRatematrixMultiMCMC_C(X=X, k=k, p=p, nodes=nodes, des=des, anc=anc, names_anc=names_anc
                               , mapped_edge=mapped.edge, R=startR, mu=start$root, sd=sqrt(startvar), Rcorr=startCorr, w_mu=w_mu
                               , par_prior_mu=par_mu, den_mu=den_mu, w_sd=w_sd, par_prior_sd=par_sd, den_sd=den_sd
                               , nu=nu, sigma=sigma_array, v=v, log_file=log_file_name, mcmc_file=mcmc_file_name
                               , prob_sample_root = prob_sample_root, prob_sample_sd = prob_sample_sd
                               , gen = gen, post_seq = post_seq, write_header = write_header)

    }

    cat( paste("Finished MCMC run ", outname, ".", new.ID, "\n", sep="") )

    ## Returns 'p = 1' to indentify the results as a single R matrix fitted to the data.
    ## Returns the data, phylogeny, priors and start point to work with other functions.
    out <- list(k = k, p = p, ID = new.ID, dir = dir, outname = outname, trait.names = traits
              , regime.names = regimes, data = t(X), phy = phy, prior = prior, start = start, gen = gen
               , mcmc.par=mcmc.par, post_seq = post_seq)
    class( out ) <- "ratematrix_multi_mcmc"
    return( out )
}
