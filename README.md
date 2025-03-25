<!-- README.md is generated from README.Rmd. Please edit that file -->

# R library mut2

The aim of `mut2` is to provide functions that balance mutation matrices
and evaluate and exemplify mutation models

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
library(pedmut)
#> 
#> Attaching package: 'pedmut'
#> The following object is masked from 'package:mut2':
#> 
#>     stepwiseReversible
```

We state that, for the expected heterozygosity, “typical values for
multi-allelic forensic markers are in the range (0.60, 0.95)”. This is
based on the updated version of the database NorwegianFrequencies:

``` r
db = getFreqDatabase(KLINK::halfsib[[1]])
H = lapply(db, function(x) 1 - sum(x^2))
range(H)
#> [1] 0.6180022 0.9490541
```

### Example 3.1 SNP markers

The mutation matrix is seen not to be stationary from

``` r
p = c("1" = 0.2, "2" = 0.8)
M = mutationMatrix("equal", alleles = names(p), rate =  0.003, afreq = p)
M
#>       1     2
#> 1 0.997 0.003
#> 2 0.003 0.997
#> 
#> Model: equal 
#> Rate: 0.003 
#> Frequencies: 0.2, 0.8 
#> 
#> Stationary: No 
#> Reversible: No 
#> Lumpable: Always
```

The \`PR’ transformations preserves the expected mytation rate as
claimed

``` r
MPR = findReversible(M, method = 'PR')
MPR
#>          1        2
#> 1 0.992500 0.007500
#> 2 0.001875 0.998125
#> 
#> Model: custom 
#> Rate: 0.003 
#> Frequencies: 0.2, 0.8 
#> 
#> Stationary: Yes 
#> Reversible: Yes 
#> Lumpable: Always
```

The *f* value 0.0075/0.003 = 2.5 can alternatively be found from

``` r
fratio(M, MPR) 
#> [1] 2.5
```

The mutation rates for the two other transformations are smaller

``` r
MMH = findReversible(M, method = 'MH')
MBA = findReversible(M, method = 'BA')
c("rateMH2" = attributes(MMH)$rate, "rateBA" = attributes(MBA)$rate)
#> rateMH2  rateBA 
#> 0.00120 0.00096
```

The last two transformations can be adjusted to have the original
mutation rate (output omitted)

``` r
MMHstar = adjustReversible(M, MMH)
MBAstar = adjustReversible(M, MBA)
```

For a SNP marker, the PR-reversible model coincides with the staionary
model and therefore we can also find the transformation from

``` r
pedmut:::stabilize(M)
#>          1        2
#> 1 0.992500 0.007500
#> 2 0.001875 0.998125
#> 
#> Model: custom 
#> Rate: 0.003 
#> Frequencies: 0.2, 0.8 
#> 
#> Stationary: Yes 
#> Reversible: Yes 
#> Lumpable: Always
```

which is an implementation of the PM transformation in Familias (output
omitted)

``` r
MPM = Familias:::FamiliasLocus(frequencies = p,
                            allelenames = names(p),
                            MutationModel = "Equal",
                            MutationRate = 0.003,
                            Stabilization = "PM")
```

### Table 1,2,3

``` r
tab1 = tabfRatio(db = db,  rate = 0.001, mutmodel = "equal", relabel = F, 
                 stationary = T, nr = 1)
tab2 = tabfRatio(db = db, rate = 0.001, mutmodel = "onestep", relabel = T, 
                 stationary = F, nr = 2)
tab3 = tabfRatio(db = db, rate = 0.001, rate2 = 0.00001, range = 0.1, 
                mutmodel = "stepwise", relabel = F, stationary = T, nr = 3)
```
