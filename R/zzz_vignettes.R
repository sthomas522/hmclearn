
# rebuild vignettes
if (FALSE) {
  # 1. run these commands to update vignettes with pre-compiled results
  # 2. move the image files from the figure to the vignettes directory
  # 3. delete the figure directory
  # 4. in each of the .Rmd files, copy/paste delete the figures/ folder reference

  knitr::knit("vignettes/linear_mixed_effects_hmclearn.Rmd.orig",
              output = "vignettes/linear_mixed_effects_hmclearn.Rmd")
  knitr::knit("vignettes/logistic_mixed_effects_hmclearn.Rmd.orig",
              output = "vignettes/logistic_mixed_effects_hmclearn.Rmd")
  knitr::knit("vignettes/linear_regression_hmclearn.Rmd.orig",
              output = "vignettes/linear_regression_hmclearn.Rmd")
  knitr::knit("vignettes/logistic_regression_hmclearn.Rmd.orig",
              output = "vignettes/logistic_regression_hmclearn.Rmd")
  knitr::knit("vignettes/poisson_regression_hmclearn.Rmd.orig",
              output = "vignettes/poisson_regression_hmclearn.Rmd")

}

