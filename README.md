<!-- README.md is generated from README.Rmd. Please edit that file -->

# mut2

The aim om `mut2` is to provide functions that balance mutation
matrices.

## Installation

To get the last version, install from GitHub as follows:

``` r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install mut2 from GitHub
devtools::install_github("thoree/mut2")
```

## Calculations in mutation paper

We first load some libraries:

``` r
library(mut2, quietly = T)
library(forrel, quietly = T)
```

We state that, for the expected heterozygosity, “typical values for
multi-allelic forensic markers are in the range (0.60, 0.95)”. This is
based on the updated version of the database NorwegianFrequencies in

``` r
db = getFreqDatabase(KLINK::halfsib[[1]])
H = lapply(db, function(x) 1 - sum(x^2))
range(H)
#> [1] 0.6180022 0.9490541
```

### Table 1,2,3

``` r
library(mut2)
library(forrel)
```

``` r
db = getFreqDatabase(KLINK::halfsib[[1]])
tab1 = tabfRatio(db = db,  rate = 0.001, mutmodel = "equal", relabel = F, 
                 stationary = T, nr = 1)
tab2 = tabfRatio(db = db, rate = 0.001, mutmodel = "onestep", relabel = T, 
                 stationary = F, nr = 2)
tab3 = tabfRatio(db = db, rate = 0.001, rate2 = 0.00001, range = 0.1, 
                mutmodel = "stepwise", relabel = F, stationary = T, nr = 3)
```
