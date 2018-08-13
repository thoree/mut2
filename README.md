<!-- README.md is generated from README.Rmd. Please edit that file -->
mut2
====

The aim om `mut2` is to provide functions for mutation models and exact likelihoods for pairs of individuals. Current functions include

-   `lik2X()` : Pairwise likelihood for the X-chromosome
-   `lik2()` : Pairwise likelihood for the X-chromosome
-   `ibd1.parental` : Estimates for a pair of non-inbred individuals the probabilities of paternal origin when IBD is 1
-   `makeReversible` : An irreducible and aperiodic mutation matrix is tranformed to satisfy detailed balanced (DB) using the Metropolis-Hastings algorithm.

Installation
------------

To get the lastest version, install from GitHub as follows:

``` r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install mut2 from GitHub
devtools::install_github("thoree/mut2")
```
