##' Print method for the "ratematrix_multi_mcmc" class.
##'
##' Print information about object.
##' @title Print method for the "ratematrix_multi_mcmc" class.
##' @param x The object.
##' @param ... Additional arguments. Not used here.
##' @export
print.ratematrix_multi_mcmc <- function(x, ...){
    cat("\n")
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
    cat("Use 'readMCMC' to load the MCMC chain.", "\n")
    cat("Check 'names' for more details.", "\n")
}
