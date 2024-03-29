% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeWellBehaved.R
\name{makeWellBehaved}
\alias{makeWellBehaved}
\title{Make mutation model well behaved}
\usage{
makeWellBehaved(Q, afreq)
}
\arguments{
\item{Q}{a mutation model}

\item{afreq}{A numeric vector of allele frequencies.}
}
\value{
A reversible stepwise mutation model.
}
\description{
Replace q_ij by min(q_ij,p_j)
}
\examples{
\dontrun{
# Not well behaved:
p = c(0.4, 0.6)
Q = stepwiseReversible(alleles = 1:2, afreq = p, rate = 0.6, range = 0.1)
makeWellBehaved(Q, afreq = p)
}
}
\author{
Thore Egeland.
}
