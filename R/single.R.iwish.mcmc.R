##' Runs a Markov chain Monte Carlo (MCMC) chain to estimate the posterior distribution of a single
##'    evolutionary rate matrix (R) fitted to the phylogeny.
##'
##' MCMC using the inverse-Wishart as a proposal distribution for the covariance matrix and a simple sliding
##'    window for the phylogenetic mean in the random walk Metropolis-Hastings algorithm. At the moment the
##'    function only applies the 'rpf' method for calculation of the log likelihood implemented in the
##'    package 'mvMORPH'. Future versions should offer different log likelihood methods for the user.
##' @title MCMC for a single evolutionary rate matrix.
##' @param X matrix. A matrix with the data. 'rownames(X) == phy$tip.label'.
##' @param phy simmap phylo. A phylogeny of the class "simmap" from the package 'phytools'. Function uses the location information for a number of traits equal to the number of fitted matrices.
##' @param start list. Element [[1]] is the starting value for the phylogenetic mean and element [[2]] is the starting value for the R matrix. Different from 'multi.R.iwish.mcmc.R' this is a single matrix and not a list.
##' @param prior list. Produced by the output of the function 'make.prior.barnard' or 'make.prior.diwish'. First element of the list [[1]] is a prior function for the log density of the phylogenetic mean and the second element [[2]] is a prior function for the evolutionary rate matrix (R).
##' @param gen numeric. Number of generations of the MCMC.
##' @param v numeric. Degrees of freedom parameter for the inverse-Wishart proposal distribution for the evolutionary rate matrix.
##' @param w numeric. Width of the sliding window proposal step for the phylogenetic mean.
##' @param prop vector. The proposal frequencies. Vector with two elements (each between 0 and 1). First is the probability that the phylogenetic mean will be sampled for a proposal step at each generation, second is the probability that the evolutionary rate matrix will be updated instead.
##' @param chunk numeric. Number of generations that the MCMC chain will be stored in memory before writing to file. At each 'chunk' generations the function will write the block stored in memory to a file and erase all but the last generation, which is used to continue the MCMC chain.
##' @param dir string. Directory to write the files, absolute or relative path. If 'NULL' then output is written to the directory where R is running (see 'getwd()'). If a directory path is given, then function will test if the directory exists and use it. If directiory does not exists the function will try to create one.
##' @param outname string. Name pasted to the files. Name of the output files will start with 'outname'.
##' @param IDlen numeric. Set the length of the unique numeric identifier pasted to the names of all output files. This is set to prevent that multiple runs with the same 'outname' running in the same directory will be lost.Default value of 5 numbers, something between 5 and 10 numbers should be good enough. IDs are generated randomly using the function 'sample'.
##' @return Fuction creates files with the MCMC chain. Additionally it returns a list object with information from the analysis to be used by other functions. This list is refered as the 'out' parameter in those functions. The list is composed by: 'acc_ratio' numeric vector with 0 when proposal is rejected and non-zero when proposals are accepted. 1 indicates that root value was accepted, 2 indicates that the evolutionary rate matrix was updated; 'run_time' in seconds; 'k' the number of matrices fitted to the tree. This value will always be 1 for this function, by see 'multi.R.iwish.mcmc'; 'p' the number of traits in the analysis; 'ID' the identifier of the run; 'dir' directory were output files were saved; 'outname' the name of the chain, appended to the names of the files; 'trait.names' A vector of names of the traits in the same order as the rows of the R matrix, can be used as the argument 'leg' for the plotting function 'make.grid.plot'; 'data' the original data for the tips; 'phy' the phylogeny; 'prior' the list of prior functions; 'start' the list of starting parameters for the MCMC run; 'gen' the number of generations of the MCMC.
##' @export
single.R.iwish.mcmc <- function(X, phy, start, prior, gen, v, w, prop=c(0.3,0.7), chunk, dir=NULL, outname="single_R_fast", IDlen=5){
    ## MCMC using the inverse-Wishart as a proposal for the random walk Metropolis-Hastings.
    ## Estimate the single evolutionary matrix model by Revell and Harmon 2008.
    ## The MCMC also make proposals for the phylogenetic mean.    
    ## Function writes to file. Creates an unique identifier for the files to prevent the user to
    ##      append one MCMC with other working in the same directory.
    ## An future version will be able to incorporate measument error originated from population-level data.
    ## X = A matrix with the data. rownames == phy$tip.label.
    ## phy = The phylogeny object.
    ## start = List with [[1]] for the phylogenetic mean and [[2]] for the R matrix.
    ## prior = List with [[1]] a function for the phylogenetic mean and [[2]] a function for the R matrix.
    ## gen = number of generations.
    ## v = tunning parameter for the inverse-Wishart proposal density.
    ## w = width parameter for the uniform sliding window proposal density.
    ## prop = vector with two floating points numbers. First is the probability that the phylogenetic
    ##        mean will be sampled for a proposal step, second is the evolutionary matrix (R) chance.
    ## chunk = number of generations that chain will be in memory before writing to file.
    ##       'chunk' need to be a fraction of the number of generations 'gen'. If 'gen/chunk != integer'
    ##       then the code will stop and return an error message.
    ## dir = directory to write the files. Optional argument. If empty then output is written to the
    ##       directory where R is running. If a directory path is given, then test if the directory
    ##       phylo.mean.step.fastexists and use it. If directiory does not exists the function creates it.
    ## outname = idetifier name of the run. files are goint to start with this string.
    ## IDlen = lenght of the unique identifier added to the name of the files. This prevent that two runs
    ##       of the MCMC will append results to the same files. Something between 5 and 10 is good enough.

    ## Verify the directory:
    if( is.null(dir) ){
        dir <- "."
    } else{
        dir.create(file.path(dir), showWarnings = FALSE)
    }

    ## Creates data cache:
    cache.data <- list()
    cache.data$n <- length(phy$tip.label) ## Number of tips.
    cache.data$k <- dim(start[[2]])[2] ## Number of traits.
    cache.data$C <- vcv(phy) ## The phylogenetic covariance matrix.
    ## This version will need to design matrix (D).
    cache.data$D <- matrix(0, nrow = cache.data$n*cache.data$k, ncol = cache.data$k)
    for(i in 1:cache.data$k) cache.data$D[((cache.data$n*(i-1))+1):(cache.data$n*i),i] <- 1
    cache.data$X <- X[rownames(cache.data$C),] ## Matching rownames of X and C.
    cache.data$traits <- colnames(X) ## Get names for the traits.
    #cache.data$y <- matrix( c( as.matrix(cache.data$X) ) ) ## Column vector format.

    ## Creates MCMC chain cache:
    ## Here trying to initialize the chain cache with the correct number of elements in the list.
    ## Note that the length is dependent on the 'chunk' and not the total of the 'gen'.
    ## This is a better approach to the loop.
    cache.chain <- list()
    cache.chain$chain <- vector(mode="list", length=chunk+1) ## Chain list.
    cache.chain$chain[[1]] <- start ## Starting value for the chain.
    cache.chain$acc <- vector(mode="integer", length=gen) ## Vector for acceptance ratio.
    ## Acceptance ratio is not recycled at each 'chunk'.
    cache.chain$acc[1] <- 1 ## Represents the starting value.
    ## The loglik function from mvMORPH does not need the b vector, but the root value.
    #cache.chain$root.curr <- as.vector(cache.chain$chain[[1]][[1]])
    cache.chain$lik <- vector(mode="numeric", length=chunk+1) ## Lik vector.

    cache.chain$lik[1] <- singleR.loglik(X=cache.data$X, root=as.vector(cache.chain$chain[[1]][[1]])
                                       , R=cache.chain$chain[[1]][[2]], C=cache.data$C, D=cache.data$D
                                       , n=cache.data$n, r=cache.data$k, method="rpf") ## Lik start value.
    ## cache.chain$lik[1] <- singleR.loglik(R=cache.chain$chain[[1]][[2]], C=cache.data$C, y=cache.data$y
    ##                                    , b=cache.chain$b.curr, n=cache.data$n, r=cache.data$k) ## Lik start value.
    cache.chain$curr.root.prior <- prior[[1]](cache.chain$chain[[1]][[1]]) ## Prior log lik starting value.
    cache.chain$curr.vcv.prior <- prior[[2]](cache.chain$chain[[1]][[2]]) ## Prior log lik starting value.

    ## Generate identifier:
    ID <- paste( sample(x=1:9, size=IDlen, replace=TRUE), collapse="")

    ## Open files to write:
    files <- list(log.lik=file(file.path(dir, paste(outname,".",ID,".loglik",sep="")), open="a"),
                  root=file(file.path(dir, paste(outname,".",ID,".root",sep="")), open="a"),
                  matrix=file(file.path(dir, paste(outname,".",ID,".matrix",sep="")),open="a")
                  )

    ## Build the update.function list:
    update.function <- list(phylo.mean.step.fast, R.matrix.step.fast)

    ## Starts the clock for the MCMC loop:
    ptm <- proc.time()

    ## Calculate chunks and create write point.
    block <- gen/chunk
    ##if(is.integer(block) == FALSE) stop("The division 'gen/chunk' needs to return an INTEGER value.")

    ## Start counter for the acceptance ratio and loglik.
    count <- 2

    ## Make loop equal to the number of blocks:
    for(jj in 1:block){

        ## Loop over the generations in each chunk:
        for(i in 2:(chunk+1) ){

            ## Proposals will be sampled given the 'prop' vector of probabilities.

            ###########################################
            ## Sample which parameter is updated:
            ## 'prop' is a vector of probabilities for 'update.function' 1 or 2.
            ## 1 = phylo root and 2 = R matrix.
            up <- sample(x = c(1,2), size = 1, prob = prop)
            ###########################################

            ###########################################
            ## Update and accept reject steps:
            cache.chain <- update.function[[up]](cache.data, cache.chain, prior, w, v, iter=i, count)
            ## Update counter.
            count <- count+1
            ###########################################
            
        }
        
        ###########################################
        ## Write to file:
        cache.chain <- write.to.file(files, cache.chain, chunk)
        ###########################################
        
    }

    ## Close the connections:
    lapply(files, close)
    
    ## Stops the clock.
    time <- proc.time() - ptm

    ## Create table of acceptance ratio.
    acc.mat <- matrix(table(cache.chain$acc), nrow=1)
    colnames(acc.mat) <- c("reject","root","R")

    ## Returns 'p = 1' to indentify the results as a single R matrix fitted to the data.
    ## Returns the data, phylogeny, priors and start point to work with other functions.
    return( list(acc_ratio = acc.mat, run_time = time[1]/3200, k = cache.data$k, p = 1
               , ID = ID, dir = dir, outname = outname, trait.names = cache.data$traits, data = X
               , phy = phy, prior = prior, start = start, gen = gen) )
}
