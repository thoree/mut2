\name{duo1}
\alias{duo1}

\docType{data}
\title{
Mutation example on Familias format
}
\description{
The data is suited for likelihood calculations in Familias.
There is genotype data for two persons, A and B, for 22 STR-s,
for two different mutation models, 'Equal' and 'Proportional'.
There are two pedigrees.

}
\usage{data("duo1")}
\format{   The format is a list of 5:
 \itemize{
   \item datamatrix dataframe with  two rows, for persons A and B, 
   and 44 columns, for the 22 markers.
   \item loci.equal A list of length 22 on Familias format for the 22 loci with 'Equal' mutation model, mutation rate 0.001.
   \item loci.proportional A list of length 22 on Familias format for the 22 loci with 'Proportional' mutation model, mutation rate 0.001.
   \item markers A character vector with the name of the markers
   \item pedigrees A list of length 2 with the pedigreeson Familias 
 }
}
\details{
The data can easiliy be transformed suitable for the R libraries
\code{pedtools, pedprobr, forrel}
}

\examples{
\dontrun{
data(duo1)
pedigrees = duo1$pedigrees
loci.proportional = duo1$loci.proportional 
datamatrix = duo1$datamatrix

# likelihood conditional on A father of B, first Familias, next pedrobr
library(Familias) 
AB = FamiliasPosterior(pedigrees, loci.proportional, datamatrix, ref = 2)
lik1 = AB$likelihoodsPerSystem[,1]

library(pedprobr)
library(forrel)
# Convert data with proportional model to ped format
x = forrel::Familias2ped(pedigrees, datamatrix, loci.proportional)[[1]]
nMarkers = 22
lik2 = rep(NA, nMarkers)
for(i in 1:nMarkers)
  lik2[i] = likelihood(x, i)
max(abs(lik2 - lik1)) #OK
}

}

