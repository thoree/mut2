\name{makeReversible}
\alias{makeReversible}
\alias{ER}

\title{
Makes Mutation Matrix Balanced
}
\description{
An irreducible and aperiodic mutation matrix
is tranformed to satisfy detailed balanced (DB) using the Metropolis-Hastings algorithm.
DB is important since swapping of genotypes can be performed
without changing LR if and only DB holds. The function ER calculates the expected mutation rate.
}
\usage{
makeReversible(Q, pi, alpha = NULL)
ER(Q, pi)
}

\arguments{
  \item{Q}{
 irreducible and aperiodic mutation matrix
}
  \item{pi}{
 allele frequency
}
  \item{alpha}{
Double, positive.
}
}
\details{
If alpha = NULL, alpha will be calculated to make alternative 3 matrix, see below, have same mutation rate as Q. Otherwise specified alpha is used.
}
\value{
\item{P1 }{Metropolis-Hastings balanced (MP)}
\item{P2 }{Averaged balanced}
\item{P3 }{alpha balanced}
\item{Minimum on diagonal }{Minimal diagonal values}
\item{Expected mutation rates }{Expected mutation rates for Q and P1,P2, and P3}
}

\references{
To appear
}
\author{
Thore Egeland <Thore.Egeland@nmbu.no>
}


\examples{
n = 4
pi = c(0.1, 0.3, 0.4, 0.2)
R = (n-1)*0.001
n = length(pi)
Q = matrix(ncol = 4, nrow = 4, 0.001)
diag(Q) = 0.997
makeReversible(Q, pi)
}

