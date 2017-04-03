## R package 'ratematrix'

*Daniel S. Caetano and Luke J. Harmon*

R package for the study of patterns of evolutionary correlation among two or more traits using phylogenetic trees. 'ratematrix' offers a suite of tools to estimate the evolutionary rate matrix (R) incorporating uncertainty in the form of a posterior distribution using Markov-chain Monte Carlo. The package allows for quick set-up and run of MCMC chain while also providing tools for users to customize their own MCMC chains.

For more information on the kind of models implemented here and their performance with empirical and simulated data, please refer to our article available at (http://biorxiv.org/content/early/2017/01/25/102939)

This package is under development. Although most of the functionality is already working without any problems, please contact the author (caetanods1[at]gmail.com) if you are interested in using this package for any publication. 

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
