##' Function plots the prior distribution used in the MCMC analysis.
##'
##' Function will make a plot of the prior for the evolutionary rate matrix by default. One can plot the prior for the root value instead by setting 'root' to TRUE. \cr
##' \cr
##' The prior distribution often has a different range of parameter values when compared to the posterior distribution. Depending on the prior configuration the range of the prior can be orders of magnitude larger than the posterior distribution. In this case, it is important to observe the scale of the x axis when comparing the prior and the posterior distribution. One can use the 'set.xlim' parameter to restrict the x axis for plotting the prior to be similar to the posterior distribution. However, often the region of parameter space of the posterior distribution has a low likelihood under the prior. This results in problems to take samples from that region to make the plot. This problem can be identified when the 'set.xlim' argument is changed and the plot shows only a few samples. \cr
##' @title Plot the prior distribution used in the MCMC analysis
##' @param handle the output object from the 'ratematrixMCMC' function.
##' @param n number of samples from the prior to be plotted (default is 1000).
##' @param root whether to plot the prior for the root value instead of the evolutionary rate matrix (default is FALSE).
##' @param color color for the plot (default is "black").
##' @param ... other parameters to be passed to the function 'plotRatematrix' or 'plotRootValue'. See help page for list of possible parameters.
##' @return A plot similar to 'plotRatematrix'.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
plotPrior <- function(handle, n=1000, root=FALSE, color="black", ...){

    cat( "Plotting the prior distribution. \n" )
    cat( "IMPORTANT NOTE: Pay attention to the scale of the rates when comparing the posterior and prior distribution plots!!. \n" )
    
    samples <- samplePrior(n=n, prior=handle$prior, sample.sd=TRUE)
    
    if(root){
        plotRootValue(chain=samples, color=color, ...)
    } else{
        plotRatematrix(chain=samples, p=1, colors=color, ...)
    }
}
