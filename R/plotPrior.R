##' Function uses the output object from the MCMC chain performed by 'ratematrixMCMC' to plot samples from the prior distribution.
##'
##' 
##' @title Plot the prior distribution used in the MCMC analysis.
##' @param out the output object from the 'ratematrixMCMC' function.
##' @param n number of samples from the prior to be plotted (default is 1000).
##' @param color color for the plot.
##' @param ... parameters to be passed to the function 'plotRatematrix'. See help for 'plotRatematrix' for list of possible parameters.
##' @return A plot similar to 'plotRatematrix'.
##' @export
##' @author daniel
plotPrior <- function(out, n=1000, color="black", ...){
    ## Right now it assumes that the prior for the p regimes fitted to the tree are the same. This does not need to be the case. Need to extend the function so that one can plot the prior when the prior is different between regimes.
    ## Take samples from the prior and use plotRatematrix to plot it.

    ## samplePriorSeparationLimits need to be implemented and used here.
    
    samples <- samplePrior(n=n, prior=out$prior, sample.sd=TRUE)
    plotRatematrix(chain=samples, p=1, colors=color, ...)
}
