##' Function plots the prior distribution used in the MCMC analysis.
##'
##' Function will make a plot of the prior for the evolutionary rate matrix by default. One can plot the prior for the root value instead by setting 'root' to TRUE. \cr
##' \cr
##' The prior distribution often has a range of parameter values order of magnitude larger than the posterior distribution. In this case, it is important to observe the scale of the x axis when comparing the prior and the posterior distribution. One can use 'set.xlim' parameter to restrict the x axis for plotting the prior to be similar to the posterior distribution. However, often the region of parameter space of the posterior distribution has a low likelihood under the prior. This results in problems to take samples from that region to make the plot. This problem can be identified when the 'set.xlim' argument is changed and the plot shows only a few samples. \cr
##' \cr
##' A solution for this issue will be implemented on future versions of the package. By now, pay attention to the scale of the axis when comparing plots!
##' @title Plot the prior distribution used in the MCMC analysis
##' @param handle the output object from the 'ratematrixMCMC' function.
##' @param n number of samples from the prior to be plotted (default is 1000).
##' @param root whether to plot the prior for the root value instead of the evolutionary rate matrix (default is FALSE).
##' @param color color for the plot (default is "black").
##' @param ... other parameters to be passed to the function 'plotRatematrix' or 'plotRootValue'. See help page for list of possible parameters.
##' @return A plot similar to 'plotRatematrix'.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' \donttest{
##' ## Load data
##' data(centrarchidae)
##' ## Run MCMC. This is just a very short chain.
##' handle <- ratematrixMCMC(data=centrarchidae$data, phy=centrarchidae$phy.map, gen=1000)
##' ## Load posterior distribution, make plots and check the log.
##' posterior <- readMCMC(handle, burn=0.25, thin=1)
##' plotRatematrix(posterior)
##' plotRootValue(posterior)
##' plotPrior(handle)
##' plotPrior(handle, root=TRUE)
##' logAnalyzer(handle)
##' }
plotPrior <- function(handle, n=1000, root=FALSE, color="black", ...){
    
    samples <- samplePrior(n=n, prior=handle$prior, sample.sd=TRUE)
    
    if(root){
        plotRootValue(chain=samples, color=color, ...)
    } else{
        plotRatematrix(chain=samples, p=1, colors=color, ...)
    }
}
