##' Print method for the "ratematrix_multi_mcmc" class.
##'
##' Print information about object.
##' @title Print method for the "ratematrix_multi_mcmc" class.
##' @param x The object.
##' @export
print.ratematrix_multi_mcmc <- function(x){
    cat("MCMC chain with multiple regimes","\n")
    cat("Number of traits: ", x$k, "\n")
    cat("Number of species: ", length( x$phy$tip.label ), "\n")
    cat("Number of regimes: ", x$p, "\n")
    cat("Number of generations: ", x$gen, "\n")
    cat("Output files: ", paste(x$outname, ".", x$ID, ".*", sep=""), "\n")
    if( x$dir == "." ){
        cat("Files directory: Same as analysis directory ('.')", "\n")
    } else{
        cat("Files directory: ", x$dir, "\n")
    }
    cat("\n")
    cat("Use 'read.multi.R.iwish' to load the MCMC chain.", "\n")
    cat("Check 'names' of x for more details.", "\n")
}


##' Print method for the "ratematrix_multi_chain" class.
##'
##' Print information about object.
##' @title Print method for the "ratematrix_multi_chain" class.
##' @param x The object.
##' @export
print.ratematrix_multi_chain <- function(x){
    ## First make some calculations:
    k <- ncol( x$root )
    post.length <- nrow( x$root )
    p <- length( x$matrix )
    cat("Posterior distribution with multiple regimes","\n")
    cat("Number of traits: ", k, "\n")
    cat("Number of regimes: ", p, "\n")
    cat("Number of posterior samples: ", post.length, "\n")
    cat("\n")
    cat("Use 'make.grid.plot' to plot the distribution.", "\n")
    cat("Use 'checkConvergence' to verify convergence.", "\n")
    cat("Use 'mergePosterior' to merge two or more posterior chains.", "\n")
    cat("Check 'names' of x for more details.", "\n")
}
