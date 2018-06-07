##' Generates a random draw from a inverse-Wishart distribution.
##'
##' @title Sample from inverse-Wishart.
##' @param v numeric. Degrees of freedom.
##' @param S matrix. The scale matrix.
##' @return matrix. Returns a sample from a inverse-Wishart distribution.
##' @noRd
riwish <- function(v, S){
    S <- solve(S)
    if (!is.matrix(S)) S <- matrix(S)
    if (v < nrow(S)) {
        stop(message = "v is less than the dimension of S in rwish().\n")
    }
    p <- nrow(S)
    CC <- chol(S)
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(stats::rchisq(p, v:(v - p + 1)))
    if (p > 1) {
        pseq <- 1:(p - 1)
        Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- stats::rnorm(p *(p - 1)/2)
    }
    out <- crossprod(Z %*% CC)
    return(solve(out))
}
