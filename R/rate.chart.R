##' Make histograms of the rates for the evolutionary rate matrix.
##'
##' Internal function to be used with 'make.grid.plot'.
##' @title Histograms of the rates.
##' @param mat List of covariance matrices.
##' @param burn proportion of burnin
##' @param thin The thinning.
##' @param true.R A point estimate to be included in the plots. This will always be a vertical line in 'red'.
##' @param file When a file string is provided the function creates a jpeg figure file istead of creating the plot in the standard device.
##' @return Plot in the standard device or a 'jpeg' file.
rate.chart <- function(mat, burn=0.5, thin=1, true.R = NULL, file = NULL){
    ## Make histograms of the rates for the evolutionary rate matrix.
    ## Number of plots equal to the number of traits.
    ## Make plots for all correlations of a list of variance-covariance matrices.
    ## Number of plots equal to the unordered combination of traits in the matrix.
    ## Disposition of the plots is calculated to prefer spread in the columns while minimizing
    ##     number of rows.
    ## mat = list of covariance matrices.
    ## mu = matrix of mu values to use as centroid for each ellipse.
    ## burn = proportion of burnin.
    ## thin = thinning, same of 'by' in function 'seq'.
    ## true.R = The true evolutionary matrix used to generate the data.
    ## file = When a file string is provided the function creates a jpeg figure file istead of
    ##        creating the plot in the standard device.
    
    ll <- length(mat)
    if(burn == 0){
        out <- round( ll / thin )
        ss <- round(seq(from = 1, to = ll, length.out = out))
    } else{
        bb <- round( ll * burn )
        out <- round( ( ll - bb ) / thin )
        ss <- round(seq(from = bb, to = ll, length.out = out))
    }
    mat.list <- mat[ss]

    rate <- lapply(mat.list, diag)
    rate <- do.call(rbind, rate)
    chart <- ncol(rate)
    cat("Number of plots = ", chart, "\n", sep = "")

    ## Calculates size of the plot chart.
    ## This will optimize for less rows in function of columns.
    col <- row <- ceiling( sqrt( chart ) )
    while(col * row > chart){
        row <- row-1
        if(col * row < chart){
                row <- row+1
                break()
            }
    }

    ## Test whether the true matrix was provided:
    if( is.matrix(true.R) ){
        true.rate <- diag(true.R)
        ## If file has filename then save plot in jpeg format.
        if( !is.null(file) ){
            jpeg(file, width=480*col , height=480*row, pointsize = 6*col)
            par(mfrow = c(row,col))
            for(i in 1:chart){
                hist.rate(rate[,i], main = paste("Trait", i, sep=" "))
                abline(v = true.rate[i], col = "red", lwd = 1.5)
            }
            dev.off()
            cat("Rates plot chart saved to ", file, "\n", sep="")
        } else{
            ## Otherwise, print to default device.
            par(mfrow = c(row,col))
            for(i in 1:chart){
                hist.rate(rate[,i], main = paste("Trait", i, sep=" "))
                abline(v = true.rate[i], col = "red", lwd = 1.5)
            }
        }
        return()
    }

    ## If file has filename then save plot in jpeg format.
    if( !is.null(file) ){
        jpeg(file, width=480*col , height=480*row, pointsize = 6*col)
        par(mfrow = c(row,col))
        for(i in 1:chart){
            hist.rate(rate[,i], main = paste("Trait", i, sep=" "))
        }
        dev.off()
        cat("Rates plot chart saved to ", file, "\n", sep="")
    } else{
        ## Otherwise, print to default device.
        par(mfrow = c(row,col))
        for(i in 1:chart){
            hist.rate(rate[,i], main = paste("Trait", i, sep=" "))
        }
    }
}
