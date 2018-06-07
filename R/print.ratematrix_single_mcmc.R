##' Print method for the "ratematrix_single_mcmc" class.
##'
##' Print information about object.
##' @title Print method for the "ratematrix_single_mcmc" class.
##' @param x The object.
##' @param ... Additional arguments. Not used here.
##' @export
print.ratematrix_single_mcmc <- function(x, ...){
    cat("\n")
    cat("MCMC chain with a single regime","\n")
    cat("Number of traits: ", x$k, "\n")
    cat("Number of species: ", length( x$phy$tip.label ), "\n")
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
