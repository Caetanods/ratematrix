# R package 'ratematrix'

*Daniel S. Caetano and Luke J. Harmon*

## Description

R package for the study of patterns of evolutionary correlation among two or more traits using phylogenetic trees. 'ratematrix' offers a suite of tools to estimate the evolutionary rate matrix (R) incorporating uncertainty in the form of a posterior distribution using Markov-chain Monte Carlo.

For more information on the kind of models implemented here and their performance with empirical and simulated data, please refer to our article published on Systematic Biology: "Caetano, D. S., and L. J. Harmon. 2019. Estimating Correlated Rates of Trait Evolution with Uncertainty. Syst Biol 68:412–429."

An overview of the package can be found at "Caetano, D. S., and L. J. Harmon. 2017. ratematrix: an R package for studying evolutionary integration among several traits on phylogenetic trees. Methods in Ecolology and Evolution http://dx.doi.org/10.1111%2F2041-210X.12826"

## Tutorials

Please check the package vignettes for tutorials on how to use the package. Vignettes can be found on the RCran page for the package "https://cran.r-project.org/web/packages/ratematrix/index.html". The pdf files can also be accessed on the directory `inst/doc/` of this repository.

## Examples from the literature

Hermansen, J. S., J. Starrfelt, K. L. Voje, and N. C. Stenseth. 2018. Macroevolutionary consequences of sexual conflict. Biology Letters 14:20180186. (http://rsbl.royalsocietypublishing.org/content/14/6/20180186)

Eliason, C. M., and J. A. Clarke. 2018. Metabolic physiology explains macroevolutionary trends in the melanic colour system across amniotes. Proccedings of the Royal Society B 285(1893) (https://doi.org/10.1098/rspb.2018.2014)

Slater, G. J., and A. R. Friscia. 2019. Hierarchy in adaptive radiation: A case study using the Carnivora (Mammalia). Evolution 73(3):524-539. (https://doi.org/10.1111/evo.13689)

Ledbetter, N. M., and R. M. Bonett. 2019. Terrestriality constrains salamander limb diversification: Implications for the evolution of pentadactyly. Journal of Evolutionary Biology 32(7):642-652. (https://doi.org/10.1111/jeb.13444)

Burress, E. D., C. M. Martinez, and P. C. Wainwright. 2020. Decoupled jaws promote trophic diversity in cichlid fishes. Evolution 74(5):950-961. (https://doi.org/10.1111/evo.13971)

Barden, P., V. Perrichot, and B. Wang. 2020. Specialized Predation Drives Aberrant Morphological Integration and Diversity in the Earliest Ants. Current Biology 30(19):3818:3824. (https://doi.org/10.1016/j.cub.2020.06.106)

Miles, M. C., E. R. Schuppe, and M. J. Fuxjager. 2020. Selection for Rhythm as a Trigger for Recursive Evolution in the Elaborate Display System of Woodpeckers. The American Naturalist 195(5) (https://doi.org/10.1086/707748)

## Install RCran (released) version

Package is now available on RCran! To install you can just type: `install.packages("ratematrix")`

Page for the package on RCran: https://cran.r-project.org/package=ratematrix

## Install development version

**For Linux and Mac:**

```{r,R.options=list(max.print=20)}
install.packages("devtools")
devtools::install_github("Caetanods/ratematrix")
```

NOTE: In case you are working with an older version of R and have problems to install the package using 'devtools', try to set this option as a workaround:
```{r,R.options=list(max.print=20)}
options(download.file.method = "wget")
```

## News and updates

**Feb-2021 (v 1.2.3): SUBMITTED TO RCran** New functions, improvements on help pages, updates to examples.

**Feb-2021 (v 1.2.3):** Adds function 'getStartPointFromPosterior' to extract the last MCMC sample from a previous run to use as the starting point of another MCMC. Adds function 'readMCMC_Mk' to read the posterior distribution of transition matrices estimated using the 'ratematrixJointMCMC' function.

**May-2019 (v 1.2.1): SUBMITTED PATCH TO RCran** Corrects issue with function 'ratematrixMCMC'. Simple (but fatal) mistake when indexing one matrix during the MCMC was breaking the chain.

**May-2019 (v 1.2): SUBMITTED TO RCran** Implements burnin and thinning during MCMC sampler. Implements function to compute the correlation among traits from the posterior samples. New version on RCran.

**Mar-2019 (v 1.2):** Updated citations for functions. Corrected BUG when merging a list of single regime MCMC chains. Corrected a potential conflict when using default priors under a 'mclapply' call to run multiple MCMC chains.

**Oct-2018 (v 1.1):** Major updates on the package. Adds new functions to perform the joint estimate of the evolutionary rate matrix and the stochastic mapping regimes together. Implements stochastic mapping algorithm on C++ (another major improvement on speed!).

**Jul-2018 (v 1.0):** Fix minor problem when using the "uniform_scaled" prior (default). The function was returning an error when the number of traits was different from the number of regimes.

**Jun-2018 (v 1.0):** MAJOR UPDATE! This is the list of changes: a) Implemented C++ code for the MCMC (Major speed improvement!), b) Now package works on Linux, Mac and Window systems!, c) Changed default prior (see help page for 'ratematrixMCMC', d) added new functions to facilitate computing the Effective Sample Size, extracting the evolutionary correlations among other things, e) updated the usability of many functions. (please check the help pages and examples for the updates), f) package now follow RCran policy. See change to 'dir' argument on function 'ratematrixMCMC'.

**May-2018 (v 0.3):** Fix minor bugs. Fix function to estimate the MCMC with no regime (single ratematrix fitted to the tree). Fix the call to the likelihood function outside the MCMC (interface issue only!). Corrects warnings from RCran checks. Almost ready for RCran.

**Feb-2018 (v 0.27.1):** Makes the x and y axes for the ellipse plots of 'plotRatematrix' isometric so it is easier to compare different traits. Imports package 'ellipse' to draw the ellipse lines (improves the code for the function).

**Aug-2017 (v 0.27):** Adds option to control the width of the proposal for sd and mean for each trait independently (see 'help(ratematrixMCMC)'). Drops option to force the MCMC to use a single rate. Fix bug that constrained argument 'w_sd' == 'w_mu'. Fix calculation of acceptance ratios performed by 'logAnalyzer'.

**Jun-2017 (v 0.25):** Add package 'vignette' with tutorial for setting a custom starting point for the MCMC.

**Jun-2017 (v 0.26):** Fix problem with 'plotRatematrix'. The ellipse lines had the x and y axes inverted. The problem is now fixed.

**Jun-2017 (v 0.25):** Add package 'vignette' with tutorial for setting a custom starting point for the MCMC.

Please contact the author (caetanods1[at]gmail.com) if you have any question, problem or suggestion for the package.
