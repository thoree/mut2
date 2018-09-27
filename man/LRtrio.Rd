\name{LRtrio}
\alias{LRtrio}
\title{
Likelihood for trio, duo and unrelated with mutations
}
\description{
Exact formulae are implemented. result can optionally be checked in the pedtools suite.
}
\usage{
LRtrio(gMan, gMother, gChild, afreq, M, check = TRUE)
}

\arguments{
  \item{gMan}{
Character vector of length 2. Genotypes of 'Man'.
}
  \item{gMother}{
Character vector of length 2. Genotypes of 'Mother'.
}
  \item{gChild}{
Character vector of length 2. Genotypes of 'Child'.
}
  \item{afreq}{
Double vector. If named, names should be indentical to dimnames of M below.
}
  \item{M}{
A mutation matrix. If dimnames are not given, names "1","2", ... are added.
M could als0 be produced by 'mutationMatrix' or 'mutationModel' in 'pedtools
}
  \item{check}{
If TRUE, answers are checked using the likelihood function of 'pedprobr'
}
}
\details{
If check=TRUE, 'pedprobr' and 'pedmut' are required.
The likelihood ratios can be obtained easily from the likelihoods.
}
\references{
To appear
}
\value{
The likelihoods for H1: the three people are unrelated,
H2: Man is unrelated of the Mother-Child offspring pair and "H3: 'Man' and 'Mother' are the parents of
Child. 
}

\author{
Thore Egeland <thore.egeland@nmbu.no>
}

\examples{
als = 1:2
p = c(0.1, 0.9)
mm = c(0.985, 0.015, 1/600, 1-1/600)
M = matrix(ncol = 2, nrow = 2, mm, byrow = TRUE)
gMan = c("2","2")
gMother = c("1","1")
gChild = c("1","1")
LRtrio(gMan, gMother, gChild,  afreq=p, M, check = FALSE)
library(pedmut)
library(pedprobr)
M = mutationMatrix("custom", matrix = M, alleles = als)
LRtrio(gMan, gMother, gChild,  afreq=p, M, check = TRUE)
}

