\name{ibd1.parental}
\alias{ibd1.parental}
\alias{uniformMap}
\title{
Estimates for a pair of non-inbred individuals the probabilities of paternal origin when IBD is 1
}
\description{
Assume  IBD is 1 for a pair of non-inbred individuals.
The function estimates the probabilities of the four
possible combinations of parental origin:
(i) paternal-paternal, (ii) paternal-maternal, (iii) maternal-paternal
and (iv) maternal-maternal.
}
\usage{
ibd1.parental(x, id.pair, Nsim, cM = 10000, verbose = F, ...)
uniformMap(Mb = NULL, cM = NULL, M = NULL, cm.per.mb = 1, 
	        chromosome = 1) 
}
\arguments{
  \item{x}{
A pedigree in the form of a \code{\link[paramlink]{linkdat}} object.
}
  \item{id.pair}{
Integer vector of length 2 giving the pair of individuals.
}
  \item{Nsim}{
Integer. Number of simulations.
}
  \item{cM}{
NULL, or a non-negative number: The distance in centiMorgan between the two loci.
}
  \item{verbose}{
logical
}
  \item{\dots}{
Further arguments to be passed on to \code{\link[IBDsim]{IBDsim}}.
}
  \item{Mb}{
Double
}
  \item{M}{
Double
}
  \item{cm.per.mb}{
Double
}
  \item{chromosome}{
integer
}
}
\details{
This is a wrapper for \code{\link[ibdsim2]{ibdsim2}}.
}
\value{
\item{alpha.sample }{Each column gives the four alpha values.}
\item{alpha.hat }{Estimates of alpha.}
}

\author{
Magnus Dehli Vigeland and Thore Egeland
}

\seealso{
\code{\link[IBDsim]{IBDsim}}.
}
\examples{
library(ibdsim2)
x = pedtools::doubleFirstCousins()
ids = c(9,10)
ibd1.parental(x, ids, 10)$alpha.hat
x = pedtools::swapSex(x,3)
ibd1.parental(x, ids, 10)$alpha.hat
}
