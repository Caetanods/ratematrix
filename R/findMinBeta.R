## Function that returns the rate categories for a discrete-Gamma distribution.
## NOTE: This is the same as the function implemented in C++. However, this functions is here only to compute the bounds for the search of Beta.
discreteGamma <- function(shape, ncats){
    ## Get the quantiles.    
    quantiles <- qgamma((1:(ncats - 1))/ncats, shape = shape, rate = shape)
    quantiles <- c(0, quantiles, Inf)

    ## Function to compute the average of the density function.
    meanGamma <- function(x){
        ## Using 'sapply' to make sure it works for vectors.
        res <- sapply(x, function(y) y * dgamma(y, shape = shape, rate = shape) )
        return( res )
    }

    rates <- vector( mode = "numeric", length = ncats )
    for( i in 1:ncats ){
        ## This is equation (10) in Yang (1994).
        int1 <- integrate(f = meanGamma, lower = quantiles[i], upper = quantiles[i+1])
        int2 <- integrate(f = dgamma, lower = quantiles[i], upper = quantiles[i+1], shape = shape, rate = shape)
        rates[i] <- int1$value / int2$value
    }

    return( rates )
}

findMinBeta <- function(ncats, precision = 0.01){
    ## Finds the lower and upper bounds for the beta parameter for the Gamma distribution.
    ## Prevents the generation of the 0 valued rates.
    ## This function searches for the min value for the Beta parameter.
    ## The max value is set to 1000. Seems to be a good guess.

    init.lb <- 0.01 ## VERY low value!
    ## Search for the minimum value of Beta that has all categories larger than 0.
    min.rate <- .Machine$double.eps * 10 # One order of magnitude larger than the minimum double.
    new <- discreteGamma(shape=init.lb, ncats=ncats)
    while( any(new < min.rate) ){
        init.lb <- init.lb + precision
        new <- discreteGamma(shape=init.lb, ncats=ncats)
    }

    lb <- init.lb
    ub <- 1000 ## Hard coded, but could be an option.
    
    bound <- setNames(c(lb, ub), c("lower","upper"))
    return( bound )
}
