##' Check if all elements are in the HPD (Highest Posterior Density) interval. Internal function to be used by 'make.grid.plot'.
##'
##' This funtion is not optimized to be called by the user.
##' @title Check if matrix is within HPD
##' @param mat matrix or dataframe.
##' @param qq list. Elements are vectors of length 2 with range of the hpd.
##' @param dd numeric. The dimension of the matrix.
##' @return boolean. Return TRUE when all elements of 'mat' is within the interval of 'qq'.
##' @author Daniel S. Caetano, \email{caetanods1@@gmail.com}
##' @seealso \code{\link{make.grid.plot}}
##' @keywords utilities
eff.chart <- function(mat, burn=0.5, thin=1, true.R = NULL, file = NULL){
    ## Make a plot of the effective dependence and effective variance.
    e.var <- effective.var(mat, burn, thin)
    e.dep <- effective.dep(mat, burn, thin)
    
    ## Test whether the true matrix was provided.
    if( is.matrix(true.R) ){
        
        k <- nrow(true.R)
        true.e.var <- det(true.R)^(1/k)
        true.mcor <- cov2cor(true.R)
        true.e.dep <- 1 - det(true.mcor)^(1/k)
        
        if( !is.null(file) ){
            jpeg(file, width = 480*2, height = 480)
            par(mfrow = c(1,2))
            hist.rate(e.var, main = "Effective variance")
            abline(v = true.e.var, col = "red", lwd = 1.5)
            hist.rate(e.dep, main = "Effective dependence")
            abline(v = true.e.dep, col = "red", lwd = 1.5)
            dev.off()
            cat("Chart saved to ", file, "\n", sep="")
        } else{
            par(mfrow = c(1,2))
            hist.rate(e.var, main = "Effective variance")
            abline(v = true.e.var, col = "red", lwd = 1.5)
            hist.rate(e.dep, main = "Effective dependence")
            abline(v = true.e.dep, col = "red", lwd = 1.5)
        }
        return()
    }
    
    if( !is.null(file) ){
        jpeg(file, width = 480*2, height = 480)
        par(mfrow = c(1,2))
        hist.rate(e.var, main = "Effective variance")
        hist.rate(e.dep, main = "Effective dependence")
        dev.off()
        cat("Chart saved to ", file, "\n", sep="")
    } else{
        par(mfrow = c(1,2))
        hist.rate(e.var, main = "Effective variance")
        hist.rate(e.dep, main = "Effective dependence")
    }
}
