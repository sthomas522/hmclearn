## Test environments
* local OS X install, R 3.6.3
* ubuntu 14.04 (on travis-ci), R 3.6.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Resubmission

Revised per feedback from CRAN volunteer, version 0.0.2

* Included references for theoretical background in DESCRIPTION.

* Alan Agresti is included as a theoretical reference, but not the cph.

* Added small examples to all exported functions.  Removed dontrun sections.

* Added 5 vignettes with additional examples. Tests are included for all of these vignettes to ensure accuracy. 

* Added tests using the testthat package for the examples, vignettes, paper (in process), and utility functions.

* Added \value to all .Rd files.

* Changed some of the parameters to match a paper in development. These changes are primarily in sample_posterior_gradients.R.

* First package submission to CRAN was on 4/1/2020, commit cc48b5283a

## Second resubmission

Revised from version 0.0.2 to version 0.0.3

* Converted all vignettes to use pre-compiled results.

* Fixed a few terminology errors in the documentation.

* Corrected some typos and special characters for the references in the vignettes. 

* Reduced the runtime of most of the tests by reducing the number of simulations.

* Second package submission to CRAN was on 6/2/2020, commmit 067724a7ef
