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

To have an idea of what this package is made for:
```
## Load all required packages:
library(ratematrix)
library(ape)
library(geiger)
library(MASS)
library(corpcor)
library(mvMORPH)
library(readr)

## Simulate data and phylogeny:
phy <- compute.brlen( rtree(100) )
rate <- matrix( c(0.5,0.3,0.3,1), nrow = 2 )
data <- sim.char(phy = phy, par = rate, nsim = 1, model = "BM")[,,1]

## Set priors:
## (I will use arbitrary prior parameters just as an example.)
prior <- make.prior.barnard(mu=2, sd=3, min=0, max=100)

## Sample starting values:
start <- list(root = mvrnorm(n=1, mu=c(2,2), Sigma=diag(x=3,nrow=ncol(rate) ) )
            , R = sample.vcv.barnard(k=ncol(rate), min=0, max=100) )

## Run the MCMC (This may take a while):
## Reduce the number fo generations ('gen') for it to finish faster.
rate.chain <- single.R.iwish.mcmc(X=data, phy=phy, start=start, prior=prior, gen=5000, v=100, w=0.25
                                , chunk=500, outname = "example")

## Print results:
rate.mcmc <- ratematrix:::read.single.R.iwish(rate.chain, burn = 0.8, thin = 10)
names(rate.mcmc)
make.grid.plot(mat1=rate.mcmc$matrix, colDiag1="blue")
```

![alt text](example.png)

Figure shows the posterior distribution of rates in the diagonal and posterior distribution of covariated evolution in the upper right. Lower left plot show an ellipse that describes the covariation structure. The vertical red lines and the red ellipse lines are the true values used to simulated data.
