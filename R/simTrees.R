##' Sample trees that show at least one clade with a given proportion of the total diversity.
##'
##' Function uses rejection sampling to generate phylogenetic trees with at least one clade that show a fixed proportion of the diversity. This can be used to perform simulations.
##' @title Generate phylogenies with focus clades.
##' @param tips numeric. The total number of tips in the phylogeny.
##' @param clade numeric. The proportion of the total tips in the focus clade.
##' @param sims numeric. Number of phylogenies to be generated.
##' @return A list with two elements. 'phy' is a "multiPhylo" object; 'nodes' is a list of the nodes in the same order as 'phy' pointing to the focus clade. The focus clade has the number of tips set by the argument 'clade'.
##' @noRd
simTrees <- function(tips, clade, sims){
    ## Function for sampling a phy with the desired number of tips that
    ##      have a clade with the exact diversity as 'clade'.
    ## tips = total number of tips of the phylogeny.
    ## clade = proportion of tips in the focus clade. numeric, min > 0, max < 1 .
    ## sims = number of simulations.
    node.list <- list()
    phy.list <- list()
    
    for(i in 1:sims){
        repeat{
            phy <- sim.bd.taxa(n = tips, numbsim = 1, lambda = 1, mu = 0)[[1]]
            nodes <- (tips+1):(tips+phy$Nnode)
            cc <- sapply( nodes, function(x) length( tips(phy, x) ) )
            nn <- nodes[ which( round(tips*clade) == cc ) ]
            if(!length(nn) == 0){
                node.list[[i]] <- nn
                phy.list[[i]] <- phy
                break
            }
        }
    }
    class(phy.list) <- "multiPhylo"
    return( list(phy = phy.list, nodes = node.list) )
}
