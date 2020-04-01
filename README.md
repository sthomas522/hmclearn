# hmclearn

We developed the R package **hmclearn** to provide users with a framework to learn the intricacies of the HMC algorithm with hands-on experience by tuning and fitting their own models, with a focus on statistical modeling in particular.  The core functions in this package include the Hamiltonian Monte Carlo (HMC) algorithm itself, including functions for the leapfrog, as well as the Metropolis-Hastings (MH) algorithm. 

While the core functions are included for both **hmc** and **mh** algorithms, users must provide their own functions for the log posterior and, for HMC, the gradient of the log posterior.  Default values are provided for the tuning parameters.  However, users will likely need to adjust the parameters for their particular applications.  

These functions were developed as part of my work learning advanced MCMC methods.  I hope this work helps those who are interested in learning these powerful methods for statistical analysis.  

# Installation

The most recent **hmclearn** release can be installed from CRAN via

```r
install.packages("hmclearn")
```

