<!-- README.md is generated from README.Rmd. Please edit that file -->

# R library mut2

The aim of `mut2` is to provide functions that make mutation models
reversible and evaluate and exemplify the resulting mutation models

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
library(pedmut, quietly = T)
library(norSTR, quietly = T)
```

### Should we recommend PR?

A simulation experiment was performed to check if the recommendation of
PR transformation could be verified. We used the function
`mut::distLRR()`. Allele frequencies and mutation parameters were
simulated as explained in the documentation with specific values shown
below.

``` r
nA = 5 #No of alleles. Allele freqs simulated uniformly
verb = F
log = T #log2 transformation of LRR
nsim = 1000
rate = c(0.0001, 0.01) #Value simulated uniformly on this interval, also below
rate2 = 10^{c(-7,-5)}
range = c(0.05,0.5)
mutMod = "equal"
seed = 1729
BA = distLRR(nAls = nA, verbose = verb, method = "BA", mutMod = mutMod, 
             log2 = log, nsim = nsim, rate = rate, seed = seed)
MH = distLRR(nAls = nA, verbose = verb, method = "MH", mutMod = mutMod, 
             log2 = log, nsim = nsim, rate = rate, seed = seed)
PR = distLRR(nAls = nA, verbose = verb, method = "PR", mutMod = mutMod,
             log2 = log, nsim = nsim, rate = rate, seed = seed)
tab1 = rbind(BA,MH, PR)
mutMod = "stepwise"
BA = distLRR(nAls = nA, verbose = verb, method = "BA", mutMod = mutMod, 
             log2 = log, nsim = nsim, rate = rate, rate2 = rate2, 
             range = range, seed = seed)
MH = distLRR(nAls = nA, verbose = verb, method = "MH", mutMod = mutMod, 
             log2 = log, nsim = nsim, rate = rate, rate2 = rate2, 
             range = range,  seed = seed) 
PR = distLRR(nAls = nA, verbose = verb, method = "PR", mutMod = mutMod,
             log2 = log, nsim = nsim, rate = rate, rate2 = rate2,
             range = range, seed = seed)
tab2 =  rbind(BA,MH, PR)
tab = rbind(tab1, tab2)
```

| mutMod   | transform |        mu |        sd |        min |     pMin |       max |     pMax |
|:---------|:----------|----------:|----------:|-----------:|---------:|----------:|---------:|
| equal    | BA        | 0.0000420 | 0.0109724 | -0.3472051 | 1.52e-05 | 0.3154393 | 2.32e-05 |
| equal    | MH        | 0.0000487 | 0.0118155 | -0.3432213 | 1.53e-05 | 0.3259195 | 2.32e-05 |
| equal    | PR        | 0.0000303 | 0.0093190 | -0.3632738 | 1.30e-05 | 0.2930584 | 2.05e-05 |
| stepwise | BA        | 0.0000996 | 0.0170487 | -0.7842868 | 1.35e-05 | 0.7072154 | 1.13e-05 |
| stepwise | MH        | 0.0001064 | 0.0177835 | -0.7603802 | 1.35e-05 | 0.7845368 | 1.13e-05 |
| stepwise | PR        | 0.0000460 | 0.0113343 | -0.9413081 | 2.00e-06 | 0.5677689 | 5.70e-06 |

Summary statistics for log2(LRR)

THe PR transformations is superior in most respects: The expected
log2(LRR), mu, is closest to zero, the standard deviation SD and the max
is smallest. However, the minimum values is also smallest and this is
the only point disfavoring PR. However, overall, the exoeriment
indicates that PR remains the transformation of choice.

### Expected heterozygosity

We state that, for the expected heterozygosity, “typical values for
multi-allelic forensic markers are in the range (0.60, 0.95)”. This is
based on the updated version of the database NorwegianFrequencies:

``` r
db = norSTR::norwayDB
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
#> Model: Equal 
#> Rate: 0.003 
#> Frequencies: 0.2, 0.8 
#> 
#> Bounded: Yes 
#> Stationary: No 
#> Reversible: No 
#> Lumpable: Always 
#> Overall rate: 0.003
```

The \`PR’ transformations preserves the expected mutation rate as
claimed

``` r
MPR = makeReversible(M, method = 'PR')
MPR
#>          1        2
#> 1 0.992500 0.007500
#> 2 0.001875 0.998125
#> 
#> Model: Custom 
#> Rate: 0.003 
#> Frequencies: 0.2, 0.8 
#> 
#> Bounded: Yes 
#> Stationary: Yes 
#> Reversible: Yes 
#> Lumpable: Always 
#> Overall rate: 0.003
```

The *f* value 0.0075/0.003 = 2.5 can alternatively be found from

``` r
fratio(M, MPR) 
#> [1] 2.5
```

The mutation rates for the two other transformations are smaller

``` r
MMH = makeReversible(M, method = 'MH')
MBA = makeReversible(M, method = 'BA')
c("rateMH2" = attributes(MMH)$rate, "rateBA" = attributes(MBA)$rate)
#> rateMH2  rateBA 
#>   0.003   0.003
```

The last two transformations can be adjusted to have the original
mutation rate (output omitted)

``` r
MMHstar = adjustRate(MMH, newrate = 0.003)
MBAstar = adjustRate(MBA, newrate = 0.003)
```

For a SNP marker, the PR-reversible model coincides with the stationary
model and therefore we can also find the transformation from

``` r
pedmut:::makeStationary(M)
#>          1        2
#> 1 0.992500 0.007500
#> 2 0.001875 0.998125
#> 
#> Model: Custom 
#> Rate: 0.003 
#> Frequencies: 0.2, 0.8 
#> 
#> Bounded: Yes 
#> Stationary: Yes 
#> Reversible: Yes 
#> Lumpable: Always 
#> Overall rate: 0.003
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

## OBSOLETE: Section 3.2 Comparing …transformations

### Table 1,2,3 (3,4,5 in current version of ms)

We next consider the tables label{tab:comparison1} -
label{tab:comparison3} , referred to in Section \`Comparing mutation…’.

``` r
Table1 = tabfRatio(db = db,  rate = 0.001, mutmodel = "equal", relabel = F, 
                 stationary = T, nr = 1)
Table2 = tabfRatio(db = db, rate = 0.001, mutmodel = "onestep", relabel = T, 
                 stationary = F, nr = 2)
Table3 = tabfRatio(db = db, rate = 0.001, rate2 = 0.00001, range = 0.1, 
                mutmodel = "stepwise", relabel = F, stationary = T, nr = 3)
```

#### Table1

Regarding “the transformation BA has smallest *f* for all markers”

``` r
table(apply(Table1[,1:3], 1, function(x) (1:3)[which.min(x)]))
```

Regarding minimum value on diagonal

``` r
apply(Table1[,4:6], 2, function(x) min(x))
```

Regarding “All transformed matrices are bounded”, this is generally true
for PR, but not necessarily for the others as adjustment is performed.
Below we extract the columns indicating boundedness and show that all
are indeed bounded in this case:

``` r
foo = tabfRatio(db = db,  rate = 0.001, mutmodel = "equal", relabel = F, 
                 stationary = T, nr = 0)
foo = foo[, c("bMH", "bBA")]
table(foo)
```

#### Table2

Regarding f-values

``` r
table(apply(Table2[,1:3], 1, function(x) (1:3)[which.min(x)]))
```

Regarding minimum diagonal values

``` r
apply(Table2[,4:6], 2, function(x) min(x))
```

Regarding boundedness

``` r
foo =  tabfRatio(db = db, rate = 0.001, mutmodel = "onestep", relabel = T, 
                 stationary = F, nr = 0)
foo =  foo[foo$int == "Y",]
foo = foo[, c("bMH", "bBA")]
table(foo)
foo[foo[,1] == "N",]
```

#### Table3

``` r
int = Table3[Table3$int == "Y",] #integer
```

f-values and minima on diagonal:

``` r
table(apply(int[,1:4], 1, function(x) (1:4)[which.min(x)]))
```

``` r
apply(int[,5:8], 2, function(x) min(x))
```

Regarding boundedness, here’s a list with at least one of three
unbounded:

``` r
foo = tabfRatio(db = db, rate = 0.001, rate2 = 0.00001, range = 0.1, 
                mutmodel = "stepwise", relabel = F, stationary = T, nr = 0)
foo =  foo[foo$int == "Y",]
foo = foo[, c("bMH", "bBA", "bDA")]
foo[foo[,1] == "N" | foo[,2] == "N" | foo[,3] == "N" ,]
```

## Section 3.3 Comparing …simulation

Code for Figure:

``` r
M = lapply(db, function(x){
  mat = mutationMatrix("equal",
                       alleles = names(x),
                       rate =  0.01,
                       rate2 = NULL,
                       range = NULL,
                       afreq = x)
  mat
})

R = lapply(M, function(x){
  mat = findReversible(x, method  = "MH")
  mat = adjustReversible(x, mat)
  mat
})


#
db = getFreqDatabase(KLINK::halfsib[[1]])
# Increase nsim
TableLRR = tabLRR(db = db, nsim = 10, seed = NULL, maks = FALSE)
saveRDS(TableLRR, "TableLRR.rds")
navn = c("MH|M", "MH|R", "BA|M", "BA|R", "PR|M", "PR|R","MH|M", "MH|R", "BA|M", "BA|R", "PR|M", "PR|R")

q99 = apply(log2(TableLRR), 2, function(x) quantile(x, probs = 0.99))
q01 = apply(log2(TableLRR) , 2, function(x) quantile(x, probs = 0.01))

pdf("C:\\Users\\theg\\OneDrive - Norwegian University of Life Sciences\\manuscripts\\mutation\\v6\\figures\\TableLRR.pdf", width = 6, height = 6)
boxplot(log2(TableLRR), horizontal = F, names = navn, ylab = "log2(D)")
abline(v = 6.5)
text(1:12, q99, "-", col = "red", cex = 2)
text(1:12, q01, "-", col = "red", cex = 2)
dev.off()
```
