##' Only works for a 2x2 matrix.
##'
##' Unclear
##' @title Unclear
##' @param mat list of covariance matrices.
##' @param mu matrix with the centroid.
##' @param burn proportion of burnin
##' @param thin thinning.
##' @param sample.line number of lines to be plotted.
##' @param n.points number of points used to calculate the ellipses.
##' @param start The start value used for the MCMC chain
##' @param true.vcv The point estimate or true value (in case of a simulation). This is going to be plotted to indicate the glocal maximum.
##' @param true.mu The point estimate or the true value (in case of a simulation).
##' @return Plot the spread of the proposals. The matrices should start close to the start value and then gravitate toward the point estimate (or true value).
vcv.proposal.plate <- function(mat, mu, burn=0, thin, sample.line, n.points = 200, start = NULL, true.vcv = NULL, true.mu = NULL){
    ## Make a plot of the proposals for the vcv matrix.
    ## Visualise the distribution of the proposals for the vcv matrix.
    ## NOTE: This will only work well when the matrix is 2x2 only.
    ## mat = list of vcv matrices.
    ## mu = a large matrix with the mu vectors.
    ## burn = burnin for the posterior
    ## thin = thinning for some plots.
    ## sample.line = number of lines to plot in the ellipse plot.
    ## n.points = the number of points to approximate the ellipse plot.
    ## start = the start parameter values for the MCMC run.
    ## true.vcv = the true vcv matrix.
    ## true.mu = the true vector of means.

    par(mfrow = c(3,3))
    ## Effective variance:
    e.var <- effective.var(mat, burn = burn, thin = thin)
    hist(e.var, main = "Effective variance")
    if( is.matrix(true.vcv) ){
        k <- nrow(true.vcv)
        true.e.var <- ( det(true.vcv) )^(1/k)
        abline(v=true.e.var, col = "red", lwd = 1.5)
    }
    ## Effective dependence:
    e.dep <- effective.dep(mat, burn = burn, thin = thin)
    hist(e.dep, main = "Effective dependence")
    if( is.matrix(true.vcv) ){
        k <- nrow(true.vcv)
        true.mcor <- cov2cor(true.vcv)
        true.e.dep <- 1 - ( det(true.mcor) )^(1/k)
        abline(v=true.e.dep, col = "red", lwd = 1.5)
    }
    ## Standard deviation plots:
    out <- sd.and.corr(mat, burn = burn, thin = thin)
    plot(x = log(out$sd[,2]), y = log(out$sd[,1]), main = "Standard deviation"
       , xlab = expression(log(sigma[2])), ylab = expression(log(sigma[1])) )
    if( is.matrix(true.vcv) ){
        true.sd <- sqrt( diag(true.vcv) )
        points(x = log(true.sd[2]), y = log(true.sd[1]), col = "red", pch = 16, cex = 1.5)
    }
    ## Standard deviation vs. correlation plot:
    corr <- sapply(out$corr, function(x) x[1,2])
    plot(x = corr, y = log(out$sd[,1]), main = "SD vs. Corr", xlab = expression(rho[12])
       , ylab = expression(log(sigma[1])) )
    if( is.matrix(true.vcv) ){
        true.corr <- cov2cor(true.vcv)
        points(x = true.corr[1,2], y = log(true.sd[1]), col = "red", pch = 16, cex = 1.5)
    }
    ## Frequency of standard deviation:
    hist(log(out$sd[,1]), freq=TRUE, xlab = expression( log(sigma[1]) ), main = "Standard deviation [1]" )
    if( is.matrix(true.vcv) ){
        abline(v=log(true.sd[1]), col = "red", lwd = 1.5)
    }

    hist(log(out$sd[,2]), freq=TRUE, xlab = expression( log(sigma[2]) ), main = "Standard deviation [2]" )
    if( is.matrix(true.vcv) ){
        abline(v=log(true.sd[2]), col = "red", lwd = 1.5)
    }
    ## Frequency of correlation:
    hist(corr, freq=TRUE, xlab = expression(rho[12]), main = "Correlation" )
    if( is.matrix(true.vcv) ){
        abline(v=true.corr[1,2], col = "red", lwd = 1.5)
    }
    ## The ellipse plot. This should be the same as the isodensity plot.
    vcv.ellipse(mat, mu, burn, thin, sample.line=sample.line, n.points)
    if( is.matrix(true.vcv) ){
        RR <- chol(true.vcv)
        angles <- seq(0, 2*pi, length.out=n.points)
        ell <- 1 * cbind(cos(angles), sin(angles)) %*% RR
        ellCtr <- sweep(ell, 2, true.mu, "+")
        points(ellCtr, col = "red", type = "l", lwd = 1.5)
    }
}
