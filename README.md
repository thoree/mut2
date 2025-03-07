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

## Table 1,2,3

```r
library(mut2)
db = getFreqDatabase(KLINK::halfsib[[1]])
tab1 = tabfRatio(db = db,  rate = 0.001, mutmodel = "equal", relabel = F, 
                 stationary = T, nr = 1)
tab2 = tabfRatio(db = db, rate = 0.001, mutmodel = "onestep", relabel = T, 
                 stationary = F, nr = 2)
tab3 = tabfRatio(db = db, rate = 0.001, rate2 = 0.00001, range = 0.1, 
                mutmodel = "stepwise", relabel = F, stationary = T, nr = 3)
`

