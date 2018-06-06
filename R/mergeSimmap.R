##' Function will merge stochastic mapped regimes together to form a new regime. This can be used to decrease the number of regimes in the phylogeny. Additionally, the function can drop all regimes and return a phylogeny of the class 'phylo'.
##'
##' The distribution of the regimes across the tree will not change. The function only modify the labels of the regimes such that two or more regimes become one (with a new label). \cr
##' \cr
##' Function takes the elements of the 'merge.regimes' vector and collapse all those regimes into a single one. The branch length associated with 'merge.regimes' are summed and assigned to the regime correspondent to the first element of the 'merge.regimes' vector. Then this new regime is renamed as 'new.regime'.\cr
##'\cr
##' If the original phylogeny has only two regimes or if 'drop.regimes' is set to TRUE, then the output will be of class 'phylo' with no regime information.
##' @title Merge two or more regimes of a 'simmap' tree
##' @param phy a phylogeny of the 'simmap' format.
##' @param merge.regimes a vector with the names of the regimes to be merged.
##' @param new.regime the name of the new regime.
##' @param drop.regimes whether to simply drop all information about the regimes and return a phylogeny of class 'phylo'.
##' @return A phylogeny of the format 'simmap' with merged regimes or a phylogeny of class 'phylo' with no regime information.
##' @export
##' @author Daniel S. Caetano and Luke J. Harmon
##' @examples
##' \donttest{
##' library( phytools ) ## Need phytools for this example.
##' data(centrarchidae)
##' plot( centrarchidae$phy.map )
##' class( centrarchidae$phy.map )
##' ## Now drop all regime information:
##' no.regime.phy <- mergeSimmap(centrarchidae$phy.map, drop.regimes=TRUE)
##' plot( no.regime.phy )
##' class( no.regime.phy )
##' ## Create a new regime with three states:
##' dt <- c(rep( c("water","earth"), each=10 ), rep("fire", times=7))
##' names(dt) <- no.regime.phy$tip.label
##' map.phy <- phytools::make.simmap(tree=no.regime.phy, x=dt)
##' plot( map.phy )
##' ## Merge two regimes:
##' merged.phy <- mergeSimmap(phy=map.phy, merge.regimes=c("water","earth"), new.regime="mud")
##' plot( merged.phy )
##' }
mergeSimmap <- function(phy, merge.regimes=NULL, new.regime=NULL, drop.regimes=FALSE){
    ## Function to merge two or more regimes into a single regime.
    ## Need to check if the regime will be unique at the end.

    if( ncol(phy$mapped.edge) == 2 || drop.regimes){
        if( ncol(phy$mapped.edge) == 2 && !drop.regimes) warning("Original phylogeny has only two regimes. Returning a 'phy' object with no regimes.")
        phy$node.label <- NULL
        phy$maps <- NULL
        phy$mapped.edge <- NULL
        phy$Q <- NULL
        phy$logL <- NULL
        class( phy ) <- "phylo"
        return( phy )
    }

    if( is.null(merge.regimes) || is.null(new.regime) ) stop("Parameters 'merge.regimes' and/or 'new.regime' have no values.\n")

    ## Define some functions to be used here.

    collapseMaps <- function(x){ ## Collapse the regimes that occurr in the same branch.
        shared <- sum( merge.regimes %in% names(x) )
        if( shared == length(merge.regimes) ){
            new.length <- sum( x[merge.regimes %in% names(x)] )
            x[ which(names(x) == merge.regimes[1]) ] <- new.length
            x <- x[ !names(x) %in% merge.regimes[-1] ]
            names(x)[ which(names(x) == merge.regimes[1]) ] <- new.regime
        }
        return(x)
    }
    
    renameMaps <- function(x){ ## Rename all the occurrences of the regime to the new name.
        for( i in 1:length(merge.regimes) ){
            to.name <- merge.regimes[i] %in% names(x)
            if( to.name ){
                names(x)[ which(names(x) == merge.regimes[i]) ] <- new.regime
            }
        }
        return(x)
    }

    ## First work with the 'maps' element.
    ## Collapse the length of the branches which merge.regimes co-occurr and rename.
    new.maps <- lapply( phy$maps, function(x) collapseMaps(x) )
    new.maps <- lapply( new.maps, function(x) renameMaps(x) )

    ## Now change the 'mapped.edge' element.
    new.mapped.edge <- phy$mapped.edge
    for(i in 1:nrow(new.mapped.edge) ){
        shared <- sum( new.mapped.edge[,merge.regimes] > 0 )
        if( shared == length( merge.regimes ) ){
            new.length <- sum( new.mapped.edge[,merge.regimes] )
            new.mapped.edge[,merge.regimes[1]] <- new.length
        }
    }   
    new.mapped.edge <- new.mapped.edge[ ,-which(colnames(new.mapped.edge) == merge.regimes[-1]) ]
    colnames(new.mapped.edge)[ which(colnames(new.mapped.edge) == merge.regimes[1]) ] <- new.regime
    phy$maps <- new.maps
    phy$mapped.edge <- new.mapped.edge

    return(phy)
}
