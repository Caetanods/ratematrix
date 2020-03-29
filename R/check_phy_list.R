## Implementing a more robust test to check if we have a single phylogeny or multiple ones:
## This test can be used across the package and it will make things simpler and more robust.
## Previous tests were getting tricked by the lack of class information.
check_phy_list <- function(x){
    ## Return a boolean.
    if( inherits(x = x, what = "multiPhylo") ){
        return( TRUE ) ## Easy case.
    }
    if( inherits(x = x[[1]], what = c("simmap","phylo")) ){
        return( TRUE ) ## List of phylo or simmaps
    }
    if( "edge" %in% names(x[[1]]) ){
        return( TRUE ) ## Should be a list of phylos.
    }
    if( inherits(x = x, what = c("phylo", "simmap") ) ){
        return( FALSE ) ## The object is a tree.
    }
    if( "edge" %in% names(x) ){
        return( FALSE ) ## Should be a single phylo, regardless of the class.
    } else{
        ## Don't know what it is!
        ## Function evaluation should not get to this point.
        stop( "Format of x is unclear!" )
    }
}
