##' Check if all elements are in the HPD (Highest Posterior Density) interval. Internal function to be used by 'make.grid.plot'.
##'
##' This funtion is not optimized to be called by the user.
##' @title Check if matrix is within HPD
##' @param mat matrix or dataframe.
##' @param qq list. Elements are vectors of length 2 with range of the hpd.
##' @param dd numeric. The dimension of the matrix.
##' @return boolean. Return TRUE when all elements of 'mat' is within the interval of 'qq'.
##' @author Daniel S. Caetano, \email{caetanods1@@gmail.com}
##' @keywords utilities
##' @noRd
checkHPD <- function(mat, qq, dd){
    ## Function to get a matrix and check if all elements are in the region.
    ## This should be an internal function to be used with make.grid.plot .
    count <- 1
    inside <- vector()
    for(i in 1:dd[1]){
        for(j in i:dd[1]){
            inside[count] <- mat[i,j] >= qq[[count]][1] & mat[i,j] <= qq[[count]][2]
            count <- count+1                    
        }
    }
    return( (sum(inside)+1) == count )
}
