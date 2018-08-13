\name{lik2}
\alias{lik2}
\title{
Pairwise likelihood with mutation
}
\description{
Detailed balance (DB) is assumed for the mutation model
and the likelihood is calculated for a pair of non-inbred individuals.
}
\usage{
lik2(g1, g2, n, p, M, kappa, alpha, theta, beta=0.5)
}
\arguments{
  \item{g1}{
Genotype, two integers  giving the alleles for individual 1.
}
  \item{g2}{
Genotype, two integers  giving the alleles for individual 2.
}
  \item{n}{
Integer vector of length 4 giving the distance between
paternal, paternal-maternal, maternal-paternal alleles and maternal.
}
  \item{p}{
Vector of real numbers. Allele frequency vector.
}
  \item{M}{
Matrix of real numbers. Mutation matrix.
}
  \item{kappa}{
Vector of real numbers describing relationship. 
IBD parameters for 0,1,2 IBD alleles.
}
  \item{alpha}{
Four probabilities, summing to 1, giving the probability,
in case IBD=1, that the alleles are paternal, paternal-maternal,
maternal-paternal, and maternal. 
}
  \item{theta}{
Real in [0,1]. Kinship coefficient. 
}
  \item{beta}{
Real in [0,1]. Probability of same parental origin when IBD=2. 
}
}
\value{
likelihood, real.
}
\details{
There are two non-inbred
individuals A and B, with genotypes a/b and c/d,
where the alleles may or may not differ.
We calculate
the likelihood assuming a relationship described by kappa.
}
\references{
Egeland, Pinto and Amorim, FSI:Genetics (2017), \code{http://dx.doi.org/10.1016/j.fsigen.2017.04.018}.
}
\author{
Thore Egeland <Thore.Egeland@nmbu.no>
}
\examples{
p <- c(0.5, 0.5)
M <- rbind(c(0.99, 0.01), c(0.01,0.99))
theta <- 0.0
n <- c(0, 1, 1, 0) 
kappa.num <- c(0, 1, 0)
kappa.den <- c(1, 0, 0)
alpha <- c(0, 0.5, 0.5, 0)
lik2(c(1,1), c(2,2), n, p, M, kappa.num, alpha, theta)
}
