# R package 'ratematrix'

*Daniel S. Caetano and Luke J. Harmon*

## THIS IS THE DEVELOPMENT VERSION!
This package version includes all the new things. Some are stable and some are not. Please send me a message if you need to use some of the development features right now. They will be soon in the stable version of the package.

'ratematrixMCMC' function in this version is implemented on C++ and has a MUCH BETTER performance.

To install you need to:

```{r,R.options=list(max.print=20)}
install.packages("devtools")
library(devtools)
install_github("Caetanods/ratematrix", ref="dev")
```

## News and updates

**May-2018 (v 0.3):** Fix minor bugs. Fix function to estimate the MCMC with no regime (single ratematrix fitted to the tree). Fix the call to the likelihood function outside the MCMC (interface issue only!). Corrects warnings from RCran checks. Almost ready for RCran.

**Jun-2017 (v 0.25):** Add package 'vignette' with tutorial for setting a custom starting point for the MCMC.

**Jun-2017 (v 0.26):** Fix problem with 'plotRatematrix'. The ellipse lines had the x and y axes inverted. The problem is now fixed.

**Aug-2017 (v 0.27):** Adds option to control the width of the proposal for sd and mean for each trait independently (see 'help(ratematrixMCMC)'). Drops option to force the MCMC to use a single rate. Fix bug that constrained argument 'w_sd' == 'w_mu'. Fix calculation of acceptance ratios performed by 'logAnalyzer'.

## Description

R package for the study of patterns of evolutionary correlation among two or more traits using phylogenetic trees. 'ratematrix' offers a suite of tools to estimate the evolutionary rate matrix (R) incorporating uncertainty in the form of a posterior distribution using Markov-chain Monte Carlo. The package allows for quick set-up and run of MCMC chain while also providing tools for users to customize their own MCMC chains.

For more information on the kind of models implemented here and their performance with empirical and simulated data, please refer to our article available at (http://biorxiv.org/content/early/2017/01/25/102939)

If you use the package please cite "Caetano, D. S., and L. J. Harmon. 2017. ratematrix: an R package for studying evolutionary integration among several traits on phylogenetic trees. Methods in Ecolology and Evolution http://dx.doi.org/10.1111%2F2041-210X.12826 "

Please contact the author (caetanods1[at]gmail.com) if you have any question, problem or suggestion for the package. The authors also watch and update the github 'issues' tab.

Installation depends on 'devtools'.
```{r,R.options=list(max.print=20)}
install.packages("devtools")
library(devtools)
install_github("Caetanods/ratematrix", build_vignettes = TRUE)
```

In case you are working with an older version of R and have problems to install the package using 'devtools', try to set this option as a workaround:
```{r,R.options=list(max.print=20)}
options(download.file.method = "wget")
```

The package offers some tutorials in form of 'vignettes'. To access use:
```{r,R.options=list(max.print=20)}
browseVignettes("ratematrix")
```

[Please check the wiki page for documentation and tutorials](https://github.com/Caetanods/ratematrix/wiki/Home).
