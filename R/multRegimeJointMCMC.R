##' Runs a Markov chain Monte Carlo (MCMC) chain to estimate the posterior distribution of two or more
##'    evolutionary rate matrices (R) fitted to the phylogeny.
##'
##' MCMC using the inverse-Wishart as a proposal distribution for the covariance matrix and a simple sliding
##'    window for the phylogenetic mean in the random walk Metropolis-Hastings algorithm. At the moment the
##'    function only applies the 'rpf' method for calculation of the log likelihood implemented in the
##'    package 'mvMORPH'. Future versions should offer different log likelihood methods for the user.
##' This version is using the pruning algoritm to c
##' @title MCMC for two or more evolutionary rate matrices.
##' @param X_BM parameter
##' @param X_Mk parameter
##' @param phy simmap phylo. A phylogeny of the class "simmap" from the package 'phytools'. Function uses the location information for a number of traits equal to the number of fitted matrices.
##' @param start list. Element [[1]] is the starting value for the phylogenetic mean and element [[2]] is the starting value for the R matrices. Element [[2]] is also a list with length equal to the number of matrices to be fitted to the data.
##' @param prior list. Produced by the output of the function 'make.prior.barnard' or 'make.prior.diwish'. First element of the list [[1]] is a prior function for the log density of the phylogenetic mean and the second element [[2]] is a prior function for the evolutionary rate matrix (R). The prior can be shared among the rate matrices or be set a different prior for each matrix. At the moment the function only produces a shared prior among the fitted matrices. Future versions will implement independent priors for each of the fitted matrices.
##' @param start_Q parameter
##' @param start_mapped.edge parameter 
##' @param prior_Mk parameter
##' @param par_prior_Mk parameter
##' @param Mk_model parameter
##' @param root_Mk parameter
##' @param smap_limit parameter
##' @param gen numeric. Number of generations of the MCMC.
##' @param v numeric. Degrees of freedom parameter for the inverse-Wishart proposal distribution for the evolutionary rate matrix.
##' @param w_sd numeric. Width of the uniform sliding window proposal for the vector of standard deviations.
##' @param w_mu numeric. Width of the uniform sliding window proposal for the vector of phylogenetic means.
##' @param w_q parameter
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
##' @param post_seq parameter
##' @return Fuction creates files with the MCMC chain. Each run of the MCMC will be identified by a unique identifier to facilitate identification and prevent the function to overwrite results when running more than one MCMC chain in the same directory. See argument 'IDlen'. The files in the directory are: 'outname.ID.loglik': the log likelihood for each generation, 'outname.ID.n.matrix': the evolutionary rate matrix n, one per line. Function will create one file for each R matrix fitted to the tree, 'outname.ID.root': the root value, one per line. \cr
##' \cr
##' Additionally it returns a list object with information from the analysis to be used by other functions. This list is refered as the 'out' parameter in those functions. The list is composed by: 'acc_ratio' numeric vector with 0 when proposal is rejected and non-zero when proposals are accepted. 1 indicates that root value was accepted, 2 and higher indicates that the first or subsequent matrices were updated; 'run_time' in seconds; 'k' the number of matrices fitted to the tree; 'p' the number of traits in the analysis; 'ID' the identifier of the run; 'dir' directory were output files were saved; 'outname' the name of the chain, appended to the names of the files; 'trait.names' A vector of names of the traits in the same order as the rows of the R matrix, can be used as the argument 'leg' for the plotting function 'make.grid.plot'; 'data' the original data for the tips; 'phy' the phylogeny; 'prior' the list of prior functions; 'start' the list of starting parameters for the MCMC run; 'gen' the number of generations of the MCMC.
##' @importFrom geiger treedata
##' @importFrom ape reorder.phylo Nnode Ntip
##' @importFrom corpcor rebuild.cov
##' @importFrom stats cov2cor
##' @noRd
multRegimeJointMCMC <- function(X_BM, X_Mk, phy, start, prior, start_Q, start_mapped.edge, prior_Mk, par_prior_Mk, Mk_model, root_Mk, smap_limit, gen, v, w_sd, w_mu, w_q, prop, dir, outname, IDlen, regimes, traits, save.handle, continue=NULL, add.gen=NULL, ID=NULL, post_seq){

    ## This is the call for the C++ function:
    ## std::string runRatematrixMCMC_jointMk_C(arma::mat X, arma::vec datMk, int k, int p, arma::vec nodes, arma::uvec des, arma::uvec anc, arma::uvec names_anc, arma::mat mapped_edge, arma::mat edge_mat, int n_nodes, arma::mat Q, double w_Q, std::string model_Q, int root_type, std::string den_Q, arma::vec par_prior_Q, arma::cube R, arma::vec mu, arma::mat sd, arma::cube Rcorr, arma::vec w_mu, arma::mat par_prior_mu, std::string den_mu, arma::mat w_sd, arma::mat par_prior_sd, std::string den_sd, arma::vec nu, arma::cube sigma, arma::vec v, std::string log_file, std::string mcmc_file, std::string Q_mcmc_file, arma::vec par_prob, int gen, int write_header){
    
    ## Get the number of regimes.
    p <- length( start[[2]] )
    
    ## Check the value of v.
    if( length( v ) > 1 ){
        if( !length( v ) == p ) stop( "Length of v need to be 1 or equal to the number of regimes." )
    }
    if( length( v ) == 1){
        v <- rep(v, times=p)
    }
    
    ## Need to transpose the data matrix.
    X_BM <- t(X_BM)
    k <- nrow(X_BM) ## Number of traits.

    ## Save the list with the MCMC parameters.
    mcmc.par <- list()
    mcmc.par$v <- v
    mcmc.par$w_sd <- w_sd
    mcmc.par$w_mu <- w_mu
    mcmc.par$w_q <- w_q
    mcmc.par$prop <- prop
    
    if( !is.null(continue) ){
        ## Use the provided ID number for the run.
        mcmc_file_name <- file.path(dir, paste(outname,".", ID, ".mcmc",sep=""))
        Q_mcmc_file_name <- file.path(dir, paste(outname,".",ID,".Q_mcmc",sep=""))
        log_file_name <- file.path(dir, paste(outname,".", ID, ".log",sep=""))
        write_header <- 0
        gen <- add.gen
    } else{
        ## Generate identifier and name for the files:
        new.ID <- paste( sample(x=1:9, size=IDlen, replace=TRUE), collapse="")
        mcmc_file_name <- file.path(dir, paste(outname,".",new.ID,".mcmc",sep=""))
        Q_mcmc_file_name <- file.path(dir, paste(outname,".",new.ID,".Q_mcmc",sep=""))
        log_file_name <- file.path(dir, paste(outname,".",new.ID,".log",sep=""))
        write_header <- 1
    }
    
    ## ord.id <- reorder.phylo(phy, order="postorder", index.only = TRUE) ## Order for traversal.
    ## mapped.edge <- start_mapped.edge ## The regimes.
    ## ## Need to take care how to match the regimes and the R matrices.
    ## anc <- phy$edge[ord.id,1] ## Ancestral edges.
    ## des <- phy$edge[ord.id,2] ## Descendent edges.
    ## edge_mat <- phy$edge[ord.id,]
    ## nodes <- unique(anc) ## The internal nodes we will traverse.

    mapped.edge <- start_mapped.edge ## The regimes.
    ## Need to take care how to match the regimes and the R matrices.
    anc <- phy$edge[,1] ## Ancestral edges.
    des <- phy$edge[,2] ## Descendent edges.
    edge_mat <- phy$edge
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
                  , regime.names = regimes, data = t(X_BM), data_Mk= t(X_Mk), phy = phy, prior = prior
                  , start = start, start_Q = start_Q, prior_Mk = prior_Mk, par_prior_Mk = par_prior_Mk
                  , gen = gen, mcmc.par = mcmc.par)
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
    runRatematrixMCMC_jointMk_C(X=X_BM, datMk=X_Mk, k=k, p=p, nodes=nodes, des=des
                              , anc=anc, names_anc=names_anc, n_tips=Ntip(phy)
                              , mapped_edge=start_mapped.edge, edge_mat=edge_mat
                              , n_nodes=Nnode(phy), Q=start_Q, w_Q=w_q
                              , model_Q=Mk_model, root_type=root_Mk, den_Q = prior_Mk
                              , par_prior_Q=par_prior_Mk, R=startR, mu=start$root
                              , sd=sqrt(startvar), Rcorr=startCorr, w_mu=w_mu
                              , par_prior_mu=par_mu, den_mu=den_mu, w_sd=w_sd
                              , par_prior_sd=par_sd, den_sd=den_sd
                              , nu=nu, sigma=sigma_array, v=v, log_file=log_file_name
                              , mcmc_file=mcmc_file_name, Q_mcmc_file=Q_mcmc_file_name
                              , par_prob = prop, gen = gen, post_seq = post_seq
                              , write_header = write_header, sims_limit = smap_limit)

    cat( paste("Finished MCMC run ", outname, ".", new.ID, "\n", sep="") )

    ## Returns the data, phylogeny, priors and start point to work with other functions.
    out <- list(k = k, p = p, ID = new.ID, dir = dir, outname = outname, trait.names = traits
              , regime.names = regimes, data = t(X_BM), data_Mk= t(X_Mk), phy = phy, prior = prior
              , start = start, start_Q = start_Q, prior_Mk = prior_Mk, par_prior_Mk = par_prior_Mk
              , gen = gen, mcmc.par = mcmc.par)
    class( out ) <- "ratematrix_multi_mcmc"
    return( out )
}
