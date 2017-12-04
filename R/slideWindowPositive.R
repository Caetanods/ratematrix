slideWindowPositive <- function(x, w){
    ## Sliding window proposal for positive real trait.
    ## This is using the reflection at the boundary.
    ## In this case is as easy as getting the absolute value of the quantity.
    ## x = the current value.
    ## w = the width parameter of the proposal.
    y <- abs( stats::runif(1, min = x - (w/2), max = x + (w/2) ) )
    return(y)
}
