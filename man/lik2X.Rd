\name{lik2X}
\alias{lik2X}
\title{
Likelihood of non-inbred male-female pair for X marker
}
\description{
Mutations are allowed.
}
\usage{
lik2X(g1, g2, n, p, M, kappa0)
}
\arguments{
  \item{g1}{
Integer. Allele of the male
}
  \item{g2}{
Integer vector of length 2. Alleles of the female
}
  \item{n}{
Integer. Number of transitions separating the pair
}
  \item{p}{
Double vector. Allele frequencies
}
  \item{M}{
Mutation matrix.
}
  \item{kappa0}{
Probability of no allele sharing
}
}
\details{
The formula is kappa0*l0 +  (1-kappa0)*l1 where l0 and l1 are the likelihoods conditioned
on sharing no or one allele IBD. The mutation matrix is assumed
to satisfy detailed balance.

}
\value{
The likelihood.
}

\author{
Thore Egeland <Thore.Egeland@nmbu.no>
}


\examples{
p = c(0.25,0.25,0.25, 0.25)
M = pedmut::mutationMatrix(alleles = 1:4, model = "proportional",
    afreq = p, rate = 0.05)
lik2X(1, 1:2, n = 1, p, M, 0) #Father-Daughter
}