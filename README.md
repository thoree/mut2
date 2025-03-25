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

The \`PR’ transformations preserves the expected mutation rate as
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

For a SNP marker, the PR-reversible model coincides with the stationary
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

## Section 3.2 Comparing …transformations

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
#> 
#>  2 
#> 50
```

Regarding minimum value on diagonal

``` r
apply(Table1[,4:6], 2, function(x) min(x))
#>       mMH       mBA       mPR 
#> 0.9949567 0.9946379 0.9370625
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
#>    bBA
#> bMH  Y
#>   Y 50
```

#### Table2

Regarding f-values

``` r
table(apply(Table2[,1:3], 1, function(x) (1:3)[which.min(x)]))
#> 
#>  1  2 
#>  1 28
```

Regarding minimum diagonal values

``` r
apply(Table2[,4:6], 2, function(x) min(x))
#>       mMH       mBA       mPR 
#> 0.9960640 0.9948979 0.8611900
```

Regarding boundedness

``` r
foo =  tabfRatio(db = db, rate = 0.001, mutmodel = "onestep", relabel = T, 
                 stationary = F, nr = 0)
foo =  foo[foo$int == "Y",]
foo = foo[, c("bMH", "bBA")]
table(foo)
#>    bBA
#> bMH  N  Y
#>   N  4  0
#>   Y  0 25
foo[foo[,1] == "N",]
#>         bMH bBA
#> TPOX      N   N
#> D11S554   N   N
#> APOAI1    N   N
#> D17S906   N   N
```

#### Table3

``` r
int = Table3[Table3$int == "Y",] #integer
```

f-values and minima on diagonal:

``` r
table(apply(int[,1:4], 1, function(x) (1:4)[which.min(x)]))
#> 
#>  2  4 
#> 15 14
```

``` r
apply(int[,5:8], 2, function(x) min(x))
#>       mMH       mBA       mPR       mDA 
#> 0.9960221 0.9948494 0.8724023 0.8738733
```

Regarding boundedness, here’s a list with at least one of three
unbounded:

``` r
foo = tabfRatio(db = db, rate = 0.001, rate2 = 0.00001, range = 0.1, 
                mutmodel = "stepwise", relabel = F, stationary = T, nr = 0)
foo =  foo[foo$int == "Y",]
foo = foo[, c("bMH", "bBA", "bDA")]
foo[foo[,1] == "N" | foo[,2] == "N" | foo[,3] == "N" ,]
#>               bMH bBA bDA
#> TPOX            N   N   N
#> D2S1338         Y   Y   N
#> D3S1358         Y   Y   N
#> D5S2800         Y   Y   N
#> D5S818          Y   Y   N
#> D7S3048         Y   Y   N
#> D9S1122         Y   Y   N
#> D10S1248        Y   Y   N
#> D11S554         N   N   N
#> APOAI1          N   N   N
#> vWA             Y   Y   N
#> D14S1434        Y   Y   N
#> D16S539         Y   Y   N
#> D17S906         N   N   N
#> D17S1301        Y   Y   N
#> D18S1364        Y   Y   N
#> D22GATA198B05   Y   Y   N
#> D22S1045        Y   Y   N
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
