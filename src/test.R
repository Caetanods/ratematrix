library( Rcpp )

sourceCpp("logLikPrunningMCMC.cpp")

logLikPrunningMCMC2(k=2, nodes=c(1,2,3) )
