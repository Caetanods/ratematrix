## R package 'ratematrix'

R package for the estimation of the evolutionary rate matrix (R) using a Bayesian approach.

See biological background in "Revell, L. J., and L. J. Harmon. 2008. Testing quantitative genetic hypotheses about the evolutionary rate matrix for continuous characters. Evolutionary Ecology Research 10:311."

Please note that the package is under development. Contact the author (caetanods1[at]gmail.com) if you are interested in using this package. However, if you want to give it a try you can use the following lines of code to install the package via github:
```
install.packages("devtools")
library(devtools)
install_github("Caetanods/ratematrix")
library("ratematrix")
```

Short tutorial to show the usage of the package:
```
## This is a short example on how to estimate the evolutionary rate matrix (R) from continuous data and a
##    phylogenetic tree using MCMC.

## You can install the package from github. For this you need the 'devtools' package.
library(devtools)
install_github("Caetanods/ratematrix")
library(ratematrix)

## Load the dependency packages. Please install the packages if necessary.
library(geiger)
library(phytools)
library(corpcor)
library(phylolm)
library(readr)

## In this example we are going to run one analysis with a single rate matrix and another with two rate
##    matrices fitted to the tree.
## First just a single matrix fitted to the tree. This analysis will be much faster than estimating multiple
##    rate matrices in the same tree.

## The R matrices for the simulation. Here using 3 traits.
base.m <- diag(1, 3)
base.m[upper.tri(base.m)] <- 0.5
base.m[lower.tri(base.m)] <- 0.5

## The phylogeny with 200 tips:
phy <- compute.brlen( rtree(n=200) )

## Generate data.
dt.one <- sim.corrs.new( tree=phy, vcv=base.m, anc=c(5,5,5) )

## Set the prior and the starting point for the analysis.
## The Barnard prior uses a separation strategy which implements a way to set a uniform prior to the rates
##     of evolution of each trait (the diagonals) and a normal distribution centered on 0 for the
##     rates of correlated evolution (the off-diagonals).
emp.mu <- colMeans(dt.one)
prior.one <- make.prior.barnard(mu=emp.mu, sd=3, min=0, max=100)

## This will set a wide normal prior for the root value centered in the empirical mean of the trait data.
## The prior for the rate matrix is marginally uninformative.
## Plot the prior distribution for the root and the prior distribution for the rate matrix:
par(mfrow = c(2,2))
hist( rnorm(n=1000, mean=emp.mu[1], sd=3) )
hist( rnorm(n=1000, mean=emp.mu[2], sd=3) )
hist( rnorm(n=1000, mean=emp.mu[3], sd=3) )
dev.off()
R.sample <- lapply(1:1000, function(x) sample.vcv.barnard(k=3, min=0, max=100) )
make.grid.plot(mat1=R.sample, colDiag1 = "gray")

## For the start value of the MCMC we are going to sample from the prior distribution:
root <- sapply(1:3, function(x) rnorm(n=1, mean=emp.mu[x], sd=3) )
R <- sample.vcv.barnard(k=p, min=0, max=100)
start.one <- list(root = root, R = R )

## Start the MCMC analysis:
## Note that the analysis is only running 1000 generations. This is not enough to reach convergence.
## Simulations show that around 1000000 generations is enough to approximate convergence.
## This function will write the mcmc generation to the directory. So it might not work if you do not
##     have writing permissions to the current directory.
mcmc.out <- single.R.iwish.mcmc(X=dt.one, phy=phy, start=start.one, prior=prior.one,
                                gen=10000, v=1000, w=0.5, chunk=500, outname="single.mat.example")

## This function will use the 'mcmc.out' object as a 'handle' and will read the files with the mcmc output.
mcmc.chain <- read.single.R.iwish(out=mcmc.out, burn=0.5, thin=10)

## Now we can plot the posterior distribution for the ratematrix:
make.grid.plot(mat1=mcmc.chain[[2]], colDiag1="blue")

## Make a traceplot of the log( likelihood ) for the model:
## Note that this is the posterior distribution, after burning.
plot(1:length(mcmc.chain[[1]]), mcmc.chain[[1]], type="l", xlab = "generations", ylab = "log.lik")
```

![alt text](example.png)

Figure shows the posterior distribution of rates in the diagonal and posterior distribution of covariated evolution in the upper right. Lower left plot show an ellipse that describes the covariation structure. The vertical red lines and the red ellipse lines are the true values used to simulated data.
