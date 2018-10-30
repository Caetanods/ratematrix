# R package 'ratematrix'

*Daniel S. Caetano and Luke J. Harmon*

## Description

R package for the study of patterns of evolutionary correlation among two or more traits using phylogenetic trees. 'ratematrix' offers a suite of tools to estimate the evolutionary rate matrix (R) incorporating uncertainty in the form of a posterior distribution using Markov-chain Monte Carlo.

For more information on the kind of models implemented here and their performance with empirical and simulated data, please refer to our article available at (http://biorxiv.org/content/early/2017/01/25/102939)

If you use the package please cite "Caetano, D. S., and L. J. Harmon. 2017. ratematrix: an R package for studying evolutionary integration among several traits on phylogenetic trees. Methods in Ecolology and Evolution http://dx.doi.org/10.1111%2F2041-210X.12826 "

## Examples in the literature

Hermansen, J. S., J. Starrfelt, K. L. Voje, and N. C. Stenseth. 2018. Macroevolutionary consequences of sexual conflict. Biology Letters 14:20180186. (http://rsbl.royalsocietypublishing.org/content/14/6/20180186)

Slater, G. J., and A. R. Friscia. 2018. Hierarchy, Morphology, and Adaptive Radiation: a Test of Osborn’s Law in the Carnivora. bioRxiv 285700. (https://www.biorxiv.org/content/early/2018/03/20/285700)

## Install RCran (released) version

Package is now available on RCran! To install you can just type: `install.packages("ratematrix")`

Page for the package on RCran: https://cran.r-project.org/web/packages/ratematrix/index.html

## Install development version

**For Linux and Mac:**

```{r,R.options=list(max.print=20)}
install.packages("devtools")
library(devtools)
install_github("Caetanods/ratematrix")
```

NOTE: In case you are working with an older version of R and have problems to install the package using 'devtools', try to set this option as a workaround:
```{r,R.options=list(max.print=20)}
options(download.file.method = "wget")
```

**For Windows:**

In this case you need to use the binary file. You can download it from this Dropbox link: https://www.dropbox.com/s/uvga0x04fakmiji/ratematrix_1.0.zip?dl=0 [please send me an email message if the file is not available]

Then install using the following lines. Of course, need to fix the path to the file to match the directory on your computer.

```{r,R.options=list(max.print=20)}
## First install the dependencies:
install.packages(pkgs=c("ape", "geiger", "coda", "corpcor", "MASS", "phylolm", "readr"
                        , "mvMORPH", "Rcpp", "ellipse"))
## Now install the 'ratematrix' package from binary:
install.packages("C:/Users/minion/Desktop/ratematrix_1.0.zip", repos = NULL, type = "win.binary")
## Try to load the package.
library( ratematrix )
## If you get an error message. Check if additional packages need to be installed.
```

## News and updates

**Oct-2018 (v 1.1.0):** Major updates on the package. Adds new functions to perform the joint estimate of the evolutionary rate matrix and the stochastic mapping regimes together. Implements stochastic mapping algorithm on C++ (another major improvement on speed!).

**Jul-2018 (v 1.0):** Fix minor problem when using the "uniform_scaled" prior (default). The function was returning an error when the number of traits was different from the number of regimes.

**Jun-2018 (v 1.0):** MAJOR UPDATE! This is the list of changes: a) Implemented C++ code for the MCMC (Major speed improvement!), b) Now package works on Linux, Mac and Window systems!, c) Changed default prior (see help page for 'ratematrixMCMC', d) added new functions to facilitate computing the Effective Sample Size, extracting the evolutionary correlations among other things, e) updated the usability of many functions. (please check the help pages and examples for the updates), f) package now follow RCran policy. See change to 'dir' argument on function 'ratematrixMCMC'.

**May-2018 (v 0.3):** Fix minor bugs. Fix function to estimate the MCMC with no regime (single ratematrix fitted to the tree). Fix the call to the likelihood function outside the MCMC (interface issue only!). Corrects warnings from RCran checks. Almost ready for RCran.

**Feb-2018 (v 0.27.1):** Makes the x and y axes for the ellipse plots of 'plotRatematrix' isometric so it is easier to compare different traits. Imports package 'ellipse' to draw the ellipse lines (improves the code for the function).

**Aug-2017 (v 0.27):** Adds option to control the width of the proposal for sd and mean for each trait independently (see 'help(ratematrixMCMC)'). Drops option to force the MCMC to use a single rate. Fix bug that constrained argument 'w_sd' == 'w_mu'. Fix calculation of acceptance ratios performed by 'logAnalyzer'.

**Jun-2017 (v 0.25):** Add package 'vignette' with tutorial for setting a custom starting point for the MCMC.

**Jun-2017 (v 0.26):** Fix problem with 'plotRatematrix'. The ellipse lines had the x and y axes inverted. The problem is now fixed.

**Jun-2017 (v 0.25):** Add package 'vignette' with tutorial for setting a custom starting point for the MCMC.

Please contact the author (caetanods1[at]gmail.com) if you have any question, problem or suggestion for the package.
