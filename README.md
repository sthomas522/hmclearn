# hmclearn

We developed the R package **hmclearn** to provide users with a framework to learn the intricacies of the HMC algorithm with hands-on experience by tuning and fitting their own models, with a focus on statistical modeling in particular.  The core functions in this package include the Hamiltonian Monte Carlo (HMC) algorithm itself, including functions for the leapfrog, as well as the Metropolis-Hastings (MH) algorithm. 

While the core functions are included for both **hmc** and **mh** algorithms, users must provide their own functions for the log posterior and, for HMC, the gradient of the log posterior.  Default values are provided for the tuning parameters.  However, users will likely need to adjust the parameters for their particular applications.  

# Installation

The most recent **hmclearn** package can be installed from github via

```r
devtools::install_github("sthomas522/hmclearn")
```

