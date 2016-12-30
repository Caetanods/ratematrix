sampleTreeClade <- function(tips, clade, sims){
    ## Function for sampling a phy with the desired number of tips that
    ##      have a clade with the exact diversity as 'clade'.
    ## tips = total number of tips of the phylogeny.
    ## clade = proportion of tips in the focus clade. numeric, min > 0, max < 1 .
    ## sims = number of simulations.
    node.list <- list()
    phy.list <- list()
    
    for(i in 1:sims){
        repeat{
            phy <- geiger::sim.bdtree(stop="taxa", n=tips)
            nodes <- (tips+1):(tips+phy$Nnode)
            cc <- sapply( nodes, function(x) length( geiger::tips(phy, x) ) )
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
