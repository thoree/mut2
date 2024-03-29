% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjustReversible.R
\name{adjustReversible}
\alias{adjustReversible}
\title{Adjusts reversible mutation matrix}
\usage{
adjustReversible(
  mutmat,
  balancedMutmat,
  method = "MH",
  afreq = NULL,
  check = TRUE
)
}
\arguments{
\item{mutmat}{Original, non balanced,  mutation matrix.}

\item{balancedMutmat}{Balanced,  mutation matrix.}

\item{method}{Character. MH, PM or BA conversions.}

\item{afreq}{A vector with allele frequencies. 
of the same length as the size of \code{mutmat}.}

\item{check}{Logical. Checks if \code{balancedMutmat} is reversible.}
}
\value{
Adjusted mutation matrix.
}
\description{
The Metropolis - Hastings conversions may give a mutation matrix \code{balancedMutmat}
with too small expected mutation rate (gamma). The \code{balancedMutmat} matrix is adjusted
to have expected mutation rate a equal to that of the original mutation matrix \code{mutmat}.
}
\details{
If \code{balancedMutmat == NULL}, \code{mutmat} is first balanced.
The adjusted balanced matrix is

\code{alpha * balancedMutmat + (1-alpha) * I}

where

\code{alpha} is the ratio of the (expected mutation) rates of the original matrix,
\code{mutmat} to the balanced version \code{balancedMutmat} and \code{I} is the identity matrix.
}
\examples{
library(pedmut)
afreq = c(0.1, 0.3, 0.4, 0.2)
names(afreq) = 1:4
mutmat = mutationMatrix("onestep", rate = 0.02, alleles = 1:4)
adj = adjustReversible(mutmat, balancedMutmat = NULL, 
                       method = "BA", afreq = afreq,  check = TRUE)
attr(mutmat, "rate") - attr(adj, "rate")
}
\seealso{
[makeReversible()]
}
\author{
Thore Egeland.
}
