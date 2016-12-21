##' Calculates the log Hastings ratio.
##'
##' Calculates the log Hastings ratio by dividing p(curr|prop)/p(prop|curr) densities. Internal function of the MCMC.
##' @title log Hastings ratio
##' @param curr.vcv matrix. The current evolutionary covariance matrix.
##' @param prop.vcv matrix. The proposed evolutionary covariance matrix.
##' @param p numeric. Number of traits in the analysis.
##' @param v numeric. Variance parameter of the inverse wishart distrubution.
##' @return numeric. The log Hastings ratio for the inverse wishart based proposal distribution.
##' @author Daniel S. Caetano, \email{caetanods1@@gmail.com}
##' @noRd
hastingsDensity <- function(curr.vcv, prop.vcv, p, v){
    ## Calculates the log Hastings ratio by dividing p(curr|prop)/p(prop|curr) densities.
    center.curr <- (v-p-1) * curr.vcv
    center.prop <- (v-p-1) * prop.vcv
    hh <- logDensityIWish(curr.vcv, v, center.prop) - logDensityIWish(prop.vcv, v, center.curr)
    return(hh)
}
