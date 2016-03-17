##' A custom function to round data.
##'
##' Custom function to round data. Internal function to be used by 'vcv.ellipse'.
##' @title Round data.
##' @param a numeric. Quantity to be rounded.
##' @param digits numeric. Number of digits.
##' @param up numeric. Control parameter.
##' @return The rounded quantity.
roundoff <- function(a, digits, up){
    ## Internal function for 'vcv.ellipse'.
    ## Calculates limits for the plots.
    if (up == -1) {
        add <- 10^(-digits)
        b <- round(a, digits = digits)
        if (a >= b) {
            x <- b
        }
        if (a < b) {
            x <- b - add
        }
    }
    if (up == 1) {
        add <- 10^(-digits)
        b <- round(a, digits = digits)
        if (a <= b) {
            x <- b
        }
        if (a > b) {
            x <- b + add
        }
    }
return(x)
}
