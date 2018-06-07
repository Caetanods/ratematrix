logLikNode <- function(ss, sigma_len , sigma_len_inv, k){
    ## Internal function to calculate the log likelihood at each node.
    ## ss = contrast. a numeric vector. float
    ## sigma_len = Sigma * branch length. numeric matrix
    ## sigma_len_inv = solve( Sigma * branch length ). numeric matrix
    ## k = numeric (int). number of traits.
    ll <- -0.5 * ( k*log(2*pi) + log(det(sigma_len)) + (ss %*% sigma_len_inv %*% ss) )
    return(ll)
}
