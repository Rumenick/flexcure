# `flexcure`
Flexible Parametric Survival Models with Cure Fraction

The `flexcure` package provides `R` implementations for analyze survival data, as from flexible parametric survival models with cure fraction. In this package, were considered the parametric accelerated failure time models: classics, the extended generalized gamma, the generalized F and Marshall Olkin extended Weibull, for to model time of susceptible elements. For the cured fraction, were considered the usual models of standard mixture and of promotion time. There are also tools that permit the exploitation of the fit and functions to generate a dataset with cure fraction, using one model specific.

## Install the development version from GitHub

To install the GitHub version you need to have the package `devtools` installed.

```
# install.packages("devtools") # run this to install the devtools package
devtools::install_github('rumenick/flexcure')
```