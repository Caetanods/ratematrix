##' Estimate time minimum time needed to run the MCMC.
##'
##' Function will estimate the time based on the computation of the log-likelihood, prior density, and the Jacobian of the proposal step. The time estimated is a minimum bound based on the processing power of the current computer. Running the MCMC in different computers might change the time. Other factors, such as writing the posterior samples to large files, can influence the time to run the MCMC.
##' @title Time estimate to complete a MCMC chain
##' @param data a matrix with the data. Each column is a different trait and species names need to be provided as rownames (rownames(data) == phy$tip.label).
##' @param phy a phylogeny of the class "simmap" with the mapped regimes for two or more R regimes OR a phylogeny of the class "phylo" for a single regime. The number of evolutionary rate matrices fitted to the phylogeny is equal to the number of regimes in 'phy'. Regime names will also be used.
##' @param gen number of generations of the complete MCMC chain. This is used to create the time estimate for the analysis.
##' @param eval.times number of replicates to compute the likelihood (default is 5). A time average across replicates will be used in order to account for the uncertainty associated with computing times.
##' @param singlerate whether the function should fit a single regime and ignore the number of regimes painted to the tree (default is FALSE).
##' @return Function returns a numeric value with the time estimate in hours and prints a message to the screen with the result.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' \donttest{
##' data(centrarchidae)
##' estimateTimeMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=10000)
##' }
estimateTimeMCMC <- function(data, phy, gen, eval.times=5, singlerate=FALSE){
    ## This version of the function does not take into account writing the MCMC to files and making the proposal steps.
    ## A quick check with an overestimation of the process (by running the whole 'ratematrixMCMC' tree pre-processing) increased time estimates by 30 min only. So, the function the way it is seems fine.
    
    ## Define some objects for the loglik to work.
    v <- 50
    w_sd <- 0.5
    w_mu <- 0.5
    prop <- c(0.1,0.9)
    chunk <- gen/100

    ## Change the data to matrix:
    if( inherits(x = data, what = "data.frame") ) data <- as.matrix( data )

    ## Check if single or multi regime:
    if( !inherits(phy, what="simmap") ){
        cat('phy is not of class "simmap". Fitting a sigle rate regime to the tree. \n')
        no_phymap <- TRUE
    } else{
        no_phymap <- FALSE
    }

    ## Use a random phylogeny from the sample, if necessary:
    if( is.list(phy[[1]]) ){
        ## Multiple phylogenies.
        binary_tree <- sapply(phy, ape::is.binary)
        if( !all(binary_tree) ) stop("Phylogeny need to be fully resolved. Try using 'multi2di' function.")
        phy <- phy[[ sample(1:length(phy), size=1) ]]
    } else{
        binary_tree <- ape::is.binary(phy)
        if( !binary_tree ) stop("Phylogeny need to be fully resolved. Try using 'multi2di' function.")
    }

    ## Separate the analysis for the single rate or for the multiple rate regime:
    if( no_phymap || singlerate ){
        ## Generate the prior distribution.
        r <- ncol( data )
        mn <- colMeans(data)
        ssd <- apply(data, 2, stats::sd)
        par.mu <- as.matrix( cbind(mn, ssd) )
        par.sd <- c(0,100)
        prior <- makePrior(r=r, p=1, den.mu="norm", par.mu=par.mu, par.sd=par.sd)
        ## Generate the state to compute the log lik.
        start <- samplePrior(n=1, prior=prior, sample.sd=FALSE)

        ## Creates data cache:
        cache.data <- list()
        cache.data$n <- length(phy$tip.label) ## Number of tips.
        cache.data$k <- ncol(data) ## Number of traits.
        cache.data$X <- data
        cache.data$traits <- colnames(data) ## Get names for the traits.
        cache.chain <- list()
        cache.chain$chain <- vector(mode="list", length=chunk+1) ## Chain list.
        cache.chain$chain[[1]] <- start ## Starting value for the chain.
        cache.chain$chain[[1]][[4]] <- rebuild.cov(r=stats::cov2cor(start[[2]]), v=start[[3]]^2)

        evaluateStepMCMC <- function(){
            logLikSingleRegime(data=cache.data, chain=cache.chain, phy=phy
                             , root=as.vector(cache.chain$chain[[1]][[1]])
                             , R=cache.chain$chain[[1]][[4]]) ## Lik start value.
            prior[[1]](cache.chain$chain[[1]][[1]]) ## Prior log lik starting value.
            prior[[2]](cache.chain$chain[[1]][[4]]) ## Prior log lik starting value.
            decom <- decompose.cov( cache.chain$chain[[1]][[2]] )
            sum( sapply(1:cache.data$k, function(x) log( decom$v[x]) ) ) * log( (cache.data$k-1)/2 )
            prior[[3]](cache.chain$chain[[1]][[3]]) ## Prior log lik starting value.
        }
        
    } else{
        ## Generate prior distribution.
        p <- ncol( phy$mapped.edge ) ## Multiple regimes.
        r <- ncol( data )
        mn <- colMeans(data)
        ssd <- apply(data, 2, stats::sd)
        par.mu <- as.matrix( cbind(mn, ssd) )
        rep.sd.regime <- rep(c(0,100), times=p)
        par.sd <- matrix(data=rep.sd.regime, nrow=p, ncol=2, byrow=TRUE)
        prior <- makePrior(r=r, p=p, den.mu="norm", par.mu=par.mu, par.sd=par.sd)
        ## Generate the state to compute the log lik.
        start <- samplePrior(n=1, prior=prior, sample.sd=FALSE)

        ## Prepare the likelihood to run:

        ## Cache for the data and for the chain:
        cache.data <- list()
        cache.chain <- list()    
        cache.data$X <- data
        cache.data$data_cor <- stats::cov2cor( stats::var( data ) ) ## This is to use the correlation of the data to draw proposals for the root.
        cache.data$k <- r ## Number of traits.
        cache.data$p <- p ## Number of regimes fitted to the tree.
        ord.id <- reorder.phylo(phy, order="postorder", index.only = TRUE) ## Order for traversal.
        cache.data$mapped.edge <- phy$mapped.edge[ord.id,] ## The regimes.
        anc <- phy$edge[ord.id,1] ## Ancestral edges.
        cache.data$des <- phy$edge[ord.id,2] ## Descendent edges.
        cache.data$nodes <- unique(anc) ## The internal nodes we will traverse.
        node.to.tip <- which( tabulate( anc[which(cache.data$des <= length(phy$tip.label))] ) == 2 )
        node.to.node <- which( tabulate( anc[which(cache.data$des > length(phy$tip.label))] ) == 2 )
        node.to.tip.node <- unique( anc )[!unique( anc ) %in% c(node.to.node, node.to.tip)]
        names(anc) <- rep(1, times=length(anc))
        names(anc)[which(anc %in% node.to.node)] <- 2
        names(anc)[which(anc %in% node.to.tip.node)] <- 3
        cache.data$anc <- anc ## This need to come with the names.
        cache.chain$chain <- vector(mode="list", length=chunk+1) ## Chain list.
        cache.chain$chain[[1]] <- start ## Starting value for the chain.
        cache.chain$chain[[1]][[4]] <- list()
        for(i in 1:cache.data$p) cache.chain$chain[[1]][[4]][[i]] <- rebuild.cov(r=stats::cov2cor(start[[2]][[i]]), v=start[[3]][[i]]^2)

        ## Function that evaluates the loglik, prior and Jacobian.
        evaluateStepMCMC <- function(){
            logLikPrunningMCMC(cache.data$X, cache.data$k, cache.data$p, cache.data$nodes, cache.data$des, cache.data$anc, cache.data$mapped.edge
                             , R=cache.chain$chain[[1]][[2]], mu=as.vector(cache.chain$chain[[1]][[1]]) )    
            prior[[1]](cache.chain$chain[[1]][[1]])
            prior[[2]](cache.chain$chain[[1]][[4]])
            decom <- lapply(1:cache.data$p, function(x) decompose.cov( cache.chain$chain[[1]][[2]][[x]] ) )
            lapply(1:cache.data$p, function(y) sum( sapply(1:cache.data$k, function(x) log( decom[[y]]$v[x]) ) ) * log( (cache.data$k-1)/2 ) )
            prior[[3]](cache.chain$chain[[1]][[3]])
        }
        
    }

    cat("\n")
    cat("Computing time...\n")
    tpm <- proc.time()
    void <- sapply(1:eval.times, function(x) evaluateStepMCMC() )
    tpm.diff <- proc.time() - tpm
    eval.time <- as.numeric( tpm.diff[3] ) / eval.times
    estimate <- ( eval.time * gen ) / 3600
    cat(paste( "Time estimated with current processing power is at least ", round(estimate, 1), " hours.\n", sep=""))
    return( estimate )
}
