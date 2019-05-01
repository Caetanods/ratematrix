##' This function uses Gibbs sampling to sample ancestral states conditioned on the posterior distribution for the other parameters of the model.
##'
##' Note that instead of randomly sampling from the posterior distribution, ancestral states are sampled for each of the samples in the posterior distribution in order. For each of the sample 'n' Gibbs samples are taken.
##' @title Sample ancestral states using Gibbs sampling
##' @param phy a phylogeny in 'phylo' format. Same used to estimate the posterior distribution for the ratematrix model.
##' @param mcmc a posterior distribution as estimated with the 'ratematrixPolytopeMCMC' function.
##' @param trait a data matrix with the traits for the tips. Same format as for the 'ratematrixMCMC' function. If `trait=NULL`, then the argument `mcmc` needs to be of type `ratematrixPolytopeMCMC`.
##' @param n number of Gibbs samples to be taken for each of the posterior samples.
##' @return an array with the samples at the nodes.
##' @author Daniel Caetano
##' @export
##' @importFrom ape reorder.phylo Ntip
sampleAncestralStates <- function(phy, mcmc, trait = NULL, n = 1){
    ## Function to perform Gibbs sampling conditioned on the samples from the posterior distribution.
    ## At the moment we are assuming a single regime model. But this can be extended to work with more regimes.

    ## Note that to take the samples for the Gibbs sampler we DO NOT need to reorder the phylogeny!

    if( !inherits(phy, what="simmap") ){
        no_phymap <- TRUE
    } else{
        no_phymap <- FALSE
    }
    
    if( no_phymap ){
        mapped.edge <- cbind( phy$edge.length, phy$edge.length )
    } else{
        mapped.edge <- phy$mapped.edge ## The regimes.
    }

    ## Need to take care how to match the regimes and the R matrices.
    anc <- phy$edge[,1] ## Ancestral edges.
    des <- phy$edge[,2] ## Descendent edges.
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

    ## Get the posterior samples from the mcmc object.
    ## k is the number of traits.
    if( is.null(trait) ){
        k <- ncol(mcmc$tip_samples) / Ntip(phy)
    } else{
        k <- ncol(trait)
    }
    
    ## R is a sample for the R matrix. This code is assuming a single regime for the moment.
    ## R needs to be a cube type object.
    R <- array( dim = c(nrow(mcmc$matrix[[1]]), ncol(mcmc$matrix[[1]]), 1) )
    R[,,1] <- mcmc$matrix[[1]]
    
    ## poly_tips is a matrix with nrow equal to the number of traits and ncols equal to the number of species.
    if( is.null(trait) ){
        poly_tips <- matrix(mcmc$tip_samples[1,], nrow = k, ncol = Ntip(phy), byrow = FALSE)
    } else{
        ## This is the matrix based on the provided data.
        poly_tips <- t( trait ) ## Need to transpose the matrix here.
    }
    sigma_vec <- diag( R[,,1] )
    ## Here we will use maximum likelihood to find create an starting state for the ancestral values.
    ## The starting state will refer to the MLE for the model with 0 covariance.
    ## We might want to run some draws before start sampling for real.
    anc_start <- t( sapply(1:k, function(x) get.ML.anc(tree = phy, x = poly_tips[x,], rate = sigma_vec[x]) ) )

    ## Start to take the samples:
    n_samples <- length( mcmc$matrix )
    poly_nodes_sample <- array( dim = c(nrow(anc_start), ncol(anc_start), n_samples * n )
                             , dimnames = list(NULL, colnames( anc_start ), NULL) )
    
    ## Make the first sample (conditioned on the starting state):
    poly_nodes_sample[,,1] <- Gibbs_sample_nodes(poly_tips = poly_tips, poly_nodes = anc_start, p = 1
                                               , nodes = nodes, des = des, anc = anc, names_anc = names_anc
                                               , mapped_edge = mapped.edge, R = R)

    if( n > 1 ){
        rates_to_visit <- rep(2:n_samples, each = n)
        rates_to_visit <- c(rep(1, times = n-1), rates_to_visit)
    } else{
        rates_to_visit <- 2:n_samples
    }

    count <- 2
    for( i in rates_to_visit ){
        ## For each sample need to get the rate and the tip states.
        R[,,1] <- mcmc$matrix[[i]]
        ## Update the tip traits.
        if( is.null(trait) ){
            ## Need to update the sample for the tips if we have them.
            poly_tips <- matrix(mcmc$tip_samples[i,], nrow = k, ncol = Ntip(phy), byrow = FALSE)
        }
        poly_nodes_sample[,,count] <- Gibbs_sample_nodes(poly_tips = poly_tips
                                                       , poly_nodes = poly_nodes_sample[,,count-1]
                                                       , p = 1, nodes = nodes, des = des, anc = anc
                                                       , names_anc = names_anc, mapped_edge = mapped.edge
                                                       , R = R)
        count <- count + 1
    }

    return( poly_nodes_sample )
}

