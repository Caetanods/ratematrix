##' Plot histogram of rates.
##'
##' Make a custom histogram to plot the posterior distribution of the rates.
##' @title Histogram of rates.
##' @param x list. Posterior distribution of the evolutionary covariance matrix (R).
##' @param main string. Title for the plot.
##' @return Plot.
##' @author Daniel S. Caetano, \email{caetanods1@@gmail.com}
hist.rate <- function(x, main = NULL){
    ## Make a custom histogram to plot the posterior distribution of the rates.
    ##  x = posterior distribution for the rate.

    h0 <- hist(x, plot = FALSE)

    ## Get probability density:
    h0$counts <- h0$counts/sum(h0$counts)
    
    xlim <- range(x)
    ylim <- c(0, ceiling(max(h0$counts)*10)/10 + 0.1 )
    plot(h0, main = "", xlab = "", ylab = "", xlim = xlim, ylim = ylim, axes = FALSE)

    ## The 95% HPD lines:
    div1 <- quantile(ecdf(x), probs = c(0.025, 0.975))
    segments(x0 = div1[1], y0 = -0.01, x1 = div1[2], y1 = -0.01, lwd = 5, lend = 1)

    ## The median lines:
    md <- median(x)
    segments(x0 = md, y0 = 0, x1 = md, y1 = (ylim[2]-0.05), lwd = 2.5, lend = 1, lty = 2)
    text(x = md, y = (ylim[2]-0.025), labels = as.character(round(md, digits = 3)))
    
    ## Make the axis and legends.
    title(main = main)
    axis(side = 2); axis(side = 1, pos = -0.01)
    mtext(expression(sigma^2), 1, line=2.5, cex = 1.5)
    mtext("Probability density", 2, line=2.5, cex = 1.5)
}
