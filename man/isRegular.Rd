% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isRegular.R
\name{isRegular}
\alias{isRegular}
\title{Checks if mutation model is well behaved}
\usage{
isRegular(M, afreq = NULL)
}
\arguments{
\item{M}{a mutation model}

\item{afreq}{A numeric vector of allele frequencies.}
}
\value{
A logical
}
\description{
Checks if m_ij <= p_j
}
\examples{
\dontrun{
# Not regular:
p = c(0.4, 0.6)
R = stepwiseReversible(alleles = 1:2, afreq = p, rate = 0.6, range = 0.1)
isRegular(R, afreq = p)
}
}
\author{
Thore Egeland.
}
