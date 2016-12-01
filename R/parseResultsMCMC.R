##' Function to read and parse the results from the simulations.
##'
##' Function to read and parse the results from the simulations. At this point the function will only work
##'   when there is two chains. This need to be extended and fixed for it to run with one or two chains.
##' @title Parse results of simulations.
##' @param out.one output of the MCMC
##' @param out.two output of the MCMC
##' @param chain.one The MCMC chain one.
##' @param chain.two The MCMC chain two.
##' @param flag Name for the run.
##' @return Plots and etc.
##' @export
parseResultsMCMC <- function(out.one, out.two, chain.one, chain.two, flag=NULL){
    ## Function to read and parse the results from the simulations.
    ## This also make plots. The function creates a directory to store all the plots and outputs.
    ## The chain is not stored. Since it is already in file.
	## This version is different from 'parseResults' because this takes the mcmc chain and not the
	##      out object to read the MCMC.
    
    ## library(readr)
    ## library(coda)
    ## library(geiger)
    ## library(phytools)
    
    if( is.null(flag) ) stop("flag need to be a string for naming the output.")

    out.one <- readRDS(out.one)
    out.two <- readRDS(out.two)
    
    mcmc.one <- readRDS(chain.one)
    mcmc.two <- readRDS(chain.two)

    ## Check convergence. Here using only a single chain.
    ## Need to separate the root and the matrices. Matrices are transformed using 'c'.
    root.one.mcmc <- coda::mcmc( mcmc.one[[1]] )
    root.two.mcmc <- coda::mcmc( mcmc.two[[1]] )
    matrix.one.mcmc <- coda::mcmc( do.call(rbind, lapply(mcmc.one[[2]], c) ) )
    matrix.two.1.mcmc <- coda::mcmc( do.call(rbind, lapply(mcmc.two[[2]][[1]], c) ) )
    matrix.two.2.mcmc <- coda::mcmc( do.call(rbind, lapply(mcmc.two[[2]][[2]], c) ) )

    ## Using the Heidelberger and Welch diagnostic for convergence.
    ## Note that better would be to use Gelman's R. But we only have one chain.
    hei.root.one <- checkHeidelTest(root.one.mcmc)
    hei.matrix.one <- checkHeidelTest(matrix.one.mcmc)
    hei.root.two <- checkHeidelTest(root.two.mcmc)
    hei.matrix.two.1 <- checkHeidelTest(matrix.two.1.mcmc)
    hei.matrix.two.2 <- checkHeidelTest(matrix.two.2.mcmc)
    passed.hei.diag <- data.frame(hei.root.one, hei.matrix.one, hei.root.two
                                , hei.matrix.two.1, hei.matrix.two.2)

    ## Calculate the Effective sample size for each of the parameters:
    eff.root.one <- range( effectiveSize(root.one.mcmc) )
    eff.matrix.one <- range( effectiveSize(matrix.one.mcmc) )
    eff.root.two <- range( effectiveSize(root.two.mcmc) )
    eff.matrix.two.1 <- range( effectiveSize(matrix.two.1.mcmc) )
    eff.matrix.two.2 <- range( effectiveSize(matrix.two.2.mcmc) )
    effRange <- data.frame(eff.root.one, eff.matrix.one, eff.root.two, eff.matrix.two.1, eff.matrix.two.2)
    rownames(effRange) <- c("min","max")

    ## Create directory to store the results
    dir.create(path = "./output_analyses", showWarnings = FALSE)

    ## Write results to file:
    converge <- rbind( as.matrix(effRange), as.matrix(passed.hei.diag) )
    colnames(converge) <- c("root.one","matrix.one","root.two","matrix.two.1","matrix.two.2")
    write.csv(converge, file=paste("./output_analyses/", flag,".converge.res.csv",sep="") )

    ## Calculate the DIC for the model test:    
    dic.one.mat <- dicMCMC(out.one, mcmc.one)
    dic.two.mat <- dicMCMC(out.two, mcmc.two)
    cat( paste("Analysis: ", flag, " | DIC1 = ", dic.one.mat, " | DIC2 = ", dic.two.mat,
               " | diff = ", dic.one.mat - dic.two.mat, "\n", sep="") )

    ## Compare with the MLE estimate. The MLE estimate function from phytools already perform
    ##      the comparison with the single matrix model.
    mle <- evol.vcv(tree=out.two$phy, X=out.two$data, maxit=4000)
    saveRDS(mle, file=paste("./output_analyses/", flag,".mle.rds",sep="") )
    cat( "MLE estimate \n" )
    cat( paste("Log lik for the simple model: ", mle$logL1, "\n", sep="") )
    cat( paste("Log lik for the two rate matrices model: ", mle$logL.multiple, "\n", sep="") )
    cat( paste("P-value for the LRT: ", mle$P.chisq, "\n", sep="") )

    ## Save DIC and LRT in a csv file:
    model.test <- matrix(data=c(dic.one.mat, dic.two.mat, dic.one.mat - dic.two.mat,
                                mle$logL1, mle$logL.multiple, mle$P.chisq),
                         nrow=3, ncol=2, byrow = FALSE)
    colnames(model.test) <- c("DIC results", "MLE results")
    rownames(model.test) <- c("One matrix", "Two matrix", "Test")
    write.csv(model.test, file=paste("./output_analyses/", flag,".test.csv",sep="") )

    ## Make plots:
    pdf( paste("./output_analyses/", flag,".gridplot.onemat.pdf",sep="") )
    plotRatematrix(mat1=mcmc.one[[2]], mat2=NULL, mle1=mle$R.single)
    dev.off()

    pdf( paste("./output_analyses/", flag,".gridplot.twomat.pdf",sep="") )
    plotRatematrix(mat1=mcmc.two[[2]][[1]], mat2=mcmc.two[[2]][[2]], mle1=mle$R.multiple[,,1]
                 , mle2=mle$R.multiple[,,2])
    dev.off()
}
