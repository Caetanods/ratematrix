## Example script for 'ratematrix'.
## Package version: 1.2.1
## Date: 10/11/2019

## Get the most up to date version of the package:
library( devtools )
devtools::install_github("Caetanods/ratematrix")

## Load packages.
library( ratematrix )
library( phytools )
library( geiger )

## Simulate a phylogenetic tree.
phy <- pbtree(n=100)
## Plot the tree:
plot( phy ); axisPhylo()

## Create three regimes and fit the regimes to the phylogeny using stochastic mapping.
reg_Q <- matrix(data = 0.1, nrow = 3, ncol = 3)
diag( reg_Q ) <- -1 * (0.1 + 0.1)
reg_trait <- sim.char(phy = phy, par = reg_Q, model = "discrete")[,,1]
reg_trait <- sapply(reg_trait, function(x) switch(x, "freshwater", "marine", "terrestrial"))
phy.map <- make.simmap(tree=phy, x=reg_trait, model="SYM")
plot( phy.map )

## Simulate two traits.
## The evolutionary correlation between the traits is larger on the marine species.
sea.rate <- rbind( c(1,0.5), c(0.5,1) )
land.rate <- rbind( c(1,0.1), c(0.1,1) )
fresh.rate <- rbind( c(1,0.1), c(0.1,1) )
rates <- list(sea.rate, land.rate, fresh.rate)
names(rates) <- c("marine","terrestrial", "freshwater")
data <- sim.corrs(tree=phy.map, vcv=rates, anc=c(2,5) )

## The data format is a matrix with a number of columns equal to the number of traits.
head( data )
## We can name our traits:
colnames( data ) <- c("trait_1", "trait_2")
head( data )

## ##################################################
## Simple analysis with 'ratematrix'
## ##################################################

## Set up and run a MCMC chain.
## First check an estimate of the time the MCMC chain will take to run:
## NOTE: This time estimate can vary in different computers and with the number of applications open at the same time. It is not very accurate at the moment.
estimateTimeMCMC(data=data, phy=phy.map, gen=100000)

## Note that this step will write files to the current working directory.
## You can change the path to another directory using the 'dir' option.
## If 'dir' is NULL, files will be written to the current directory.
handle <- ratematrixMCMC(data=data, phy=phy.map, gen=50000, dir="my_analysis")

## NOTE: The console will show "Starting MCMC ...", meaning that the MCMC is running.
## You can check the files with the posterior distribution to find if the analyses is close to an end.

## We can control the sampling frequency (thinning) and the burn-in proportion of the MCMC:
handle <- ratematrixMCMC(data=data, phy=phy.map, burn = 0.2, thin = 100, gen=100000, dir="my_analysis_thin")

## Now we can load the MCMC samples and make some plots:
## NOTE: We set the burn and thin before, so here we do not need to.
posterior <- readMCMC(handle)
logAnalyzer(handle)
plotRatematrix(posterior, colors=c("blue","brown", "orange"), , alphaDiag=0.7, alphaOff=0.7)
plotRootValue(posterior)

## A quick check of how good the posterior is compared to the true values:
plotRatematrix(posterior, colors=c("blue","brown","orange"), alphaDiag=0.7, alphaOff=0.7, point.matrix=rates, point.color=c("green","green"), point.wd=3)

## ##################################################
## Check convergence of the chains
## ##################################################

## Best way to check for convergence is to run multiple chains:
handle1 <- ratematrixMCMC(data=data, phy=phy.map, burn = 0.2, thin = 100, gen=100000, dir="my_analysis_conv")
handle2 <- ratematrixMCMC(data=data, phy=phy.map, burn = 0.2, thin = 100, gen=100000, dir="my_analysis_conv")
mcmc1 <- readMCMC( handle1 )
mcmc2 <- readMCMC( handle2 )

## Now we can use the Gelman and Rubin reduction factor.
checkConvergence( mcmc1, mcmc2 )
## Here we want to get as close as possible to values of 1 to all parameters.
## We also want large Effective Sample Sizes.

## ##################################################
## Merging the posterior distributions.
## ##################################################

## Because we ran two independent MCMCs with the same data we can merge the results in order to have a longer MCMC chain (i.e., more posterior samples to work with).
full.mcmc <- mergePosterior(mcmc1, mcmc2)
plotRatematrix( full.mcmc )
plotRootValue( full.mcmc )

## ##################################################
## Checking for difference in evolutionary integration between regimes.
## ##################################################

## We can use the posterior distribution to test for differences in the evolutionary correlation among the traits.

## First we can plot the correlation of the traits.
## NOTE: The 'plotRatematrix' function will show the covariances, not the correlations.

corr.samples <- extractCorrelation( full.mcmc )
names( corr.samples )
## This object is a list of matrices with the pairwise correlation values.

## Plot the posterior distribution of correlations for the marine species.
par(mfrow = c(1,3))
hist( x = corr.samples$marine, xlim = c(0,1), main = "Marine")
hist( x = corr.samples$terrestrial, xlim = c(0,1), main = "Terrestrial")
hist( x = corr.samples$freshwater, xlim = c(0,1), main = "Freshwater")

## We can compute a test statistics, derived from the posterior distribution to check for the difference between the evolutionary correlation and rates of evolution among the regimes.
## See more information here: Caetano, D. S., and L. J. Harmon. 2019. Estimating Correlated Rates of Trait Evolution with Uncertainty. Syst Biol 68:412–429.

## Test for changes in the evolutionary integration among traits:
testRatematrix( chain = full.mcmc, par = "correlation" )
testRatematrix( chain = full.mcmc, par = "correlation", median.test = TRUE)
testRatematrix( chain = full.mcmc, par = "correlation", plot = TRUE)

## Test for changes in the rates of evolution for the traits:
testRatematrix( chain = full.mcmc, par = "rates" )
testRatematrix( chain = full.mcmc, par = "rates", median.test = TRUE)
testRatematrix( chain = full.mcmc, par = "rates", plot = TRUE)

## Test for changes on all parameters of the evolutionary rate matrix:
testRatematrix( chain = full.mcmc, par = "all" )
testRatematrix( chain = full.mcmc, par = "all", median.test = TRUE)
testRatematrix( chain = full.mcmc, par = "all", plot = TRUE)

## NOTE: On the tests above, the 'median.test' option returns an average of the differences among all elements. For example, when testing for the correlation, then we get the median value of the overlap among all correlations.

## ##################################################
## Joint analyses of the predictor regimes and the evolutionary rate matrices.
## ##################################################

## The previous analyses were using a fixed predictor regime created using stochastic mapping. However, this approach is not taking into account the variation of the predictor regime.
## We can use another function to perform a MCMC to estimate both the evolutionary rate matrices and the predictor regimes.

## Now we will run the same flow with a joint estimate of the mvBM model and the discrete predictor.

joint1 <- ratematrixJointMCMC(data_BM = data, data_Mk = reg_trait, phy = phy, dir = "test_joint", gen = 10000, burn = 0.1, thin = 10, prop = c(1,1,1,1,1))
joint2 <- ratematrixJointMCMC(data_BM = data, data_Mk = reg_trait, phy = phy, dir = "test_joint", gen = 10000, burn = 0.1, thin = 10, prop = c(1,1,1,1,1))

## Read the chains
mcmc1 <- readMCMC( joint1 )
mcmc2 <- readMCMC( joint2 )

## Check the convergence
## IMPORTANT: Note that these chains DID NOT reach convergence. More generations are needed!
checkConvergence( mcmc1, mcmc2 )
## Merge the results
full.mcmc <- mergePosterior(mcmc1, mcmc2)

## Make some plots
plotRatematrix( full.mcmc )
plotRootValue( full.mcmc )

## Test for differences between the regimes (pairwise).
testRatematrix( chain = full.mcmc, par = "all" )
testRatematrix( chain = full.mcmc, par = "all", median.test = TRUE)
testRatematrix( chain = full.mcmc, par = "all", plot = TRUE)
