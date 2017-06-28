# R package 'ratematrix'

*Daniel S. Caetano and Luke J. Harmon*

## News and updates

**Jun-2017 (v 0.25):** Add package 'vignette' with tutorial for setting a custom starting point for the MCMC.

**Jun-2017 (v 0.26):** Fix problem with 'plotRatematrix'. The ellipse lines had the x and y axes inverted. The problem is now fixed.

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

[Please check the wiki page for documentation and tutorials](https://github.com/Caetanods/ratematrix/wiki/Home).
