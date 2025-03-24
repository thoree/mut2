#' Pairwise likelihood with mutation
#'
#' Reversibility is assumed for the mutation model and the likelihood is
#' calculated for a pair of non-inbred individuals.
#'
#' There are two non-inbred individuals A and B, with genotypes a/b and c/d,
#' where the alleles may or may not differ. We calculate the likelihood assuming
#' a relationship described by kappa.
#'
#' @param g1 Genotype, two integers giving the alleles for individual 1.
#' @param g2 Genotype, two integers giving the alleles for individual 2.
#' @param n Integer vector of length 4 giving the distance between paternal-paternal,
#'   paternal-maternal, maternal-paternal and maternal-maternal alleles.
#' @param p Vector of real numbers. Allele frequency vector.
#' @param M Matrix of real numbers. Mutation matrix.
#' @param kappa Vector of real numbers describing relationship.  IBD parameters
#'   for 0,1,2 IBD alleles.
#' @param alpha Four probabilities, summing to 1, giving the probability, in
#'   the case IBD = 1, that the alleles are paternal-paternal, paternal-maternal,
#'   maternal-paternal, and maternal-maternal.
#' @param theta Real in `[0,1]`. Kinship coefficient.
#' @return Likelihood.
#' @author Thore Egeland <Thore.Egeland@@nmbu.no>
#' @references Egeland, Pinto and Amorim, FSI:Genetics (2017),
#'   \doi{http://dx.doi.org/10.1016/j.fsigen.2017.04.018}.
#' @export
#' @importFrom expm %^%
#' @examples
#' library(pedtools)
#' library(pedmut)
#' # Example 1 Parent offspring relationship
#' p = c("1" = 0.2, "2" = 0.8)
#' M = mutationMatrix("proportional", afreq = p, rate = 0.01)
#' n = c(0, 1, 1, 0)
#' kappa =  c(0, 1, 0)
#' alpha =  c(0, 0.5, 0.5, 0)
#' l1 = lik2(c(1,1), c(2,2), n, p, M, kappa, alpha)
#' # Calculated using formula
#' gamma = mut2::expectedMutationRate(M, p)
#' K = gamma/(1- sum(p^2))
#' l1.formula = p[1]^2*K*p[2]^2
#' l1 - l1.formula
#'
#' # Example 2 Double first cousins
#' n = c(0, 4, 4, 0)
#' kappa = c(9,6,1)/16
#' alpha = c(0, 0.5,0.5, 0)
#' g1 = c(1,1); g2 = c(2,2)
#' lik1 = lik2(g1, g2, n, p, M, kappa, alpha)
#' lik0 = lik2(g1, g2, n = rep(0,4), p, M,
#'             kappa = c(1,0,0), alpha = rep(0,4))
#' LR = lik1/lik0
#' # Formula
#' LR.formula = kappa[1] +
#'            kappa[2]*(1-(1-K)^4)+
#'            kappa[3]*(1-(1-K)^4) * (1-(1-K)^4)
#' LR - LR.formula

lik2 <- function(g1, g2, n, p, M, kappa, alpha, theta = 0){

  l0 <- function(a, b, c, d, p, theta) {
    pa <- p[a]
    pb <- (b==a)*theta+(1-theta)*p[b]
    pc <- (((c==a)+(c==b))*theta+(1-theta)*p[c])/(1+theta)
    pd <- (((d==a)+(d==b)+(d==c))*theta+(1-theta)*p[d])/(1+2*theta)
    2^(-(a == b) - (c == d))*4*pa*pb*pc*pd
  }

  l1 <- function(a, b, c, d, g, p, M, theta) {
    Mg <- expm::'%^%'(M, g)
    pa <- p[a]
    pb <- (b==a)*theta+(1-theta)*p[b]
    pc <- (((c==a)+(c==b))*theta+(1-theta)*p[c])/(1+theta)
    pd <- (((d==a)+(d==b))*theta+(1-theta)*p[d])/(1+theta)
    2^(-(a == b) - (c == d))*pa*pb*(pd*(Mg[a, c]+Mg[b, c])+pc*(Mg[a, d] + Mg[b,d]))
  }
  l2 <- function(a, b, c, d, nP, nM, p, M, theta) {
    MnP <- expm::'%^%'(M, nP)
    MnD <- expm::'%^%'(M, nM)
    pa <- p[a]
    pb <- (b==a)*theta+(1-theta)*p[b]
    2^(-(a == b) - (c == d))*pa*pb*
      (MnP[a, c] * MnD[b, d] + MnP[b, c] * MnD[a, d] +
         MnP[a, d] * MnD[b, c] + MnP[b, d] * MnD[a, c])
  }
  a <- g1[1]
  b <- g1[2]
  c <- g2[1]
  d <- g2[2]
  beta <- alpha[1] + alpha[4]
  beta <- c(beta, 1-beta)
  lik <- kappa[1]*l0(a, b, c, d, p, theta)+
    kappa[2] * alpha[1] * l1(a, b, c, d, g = n[1], p, M, theta) +
    kappa[2] * alpha[2] * l1(a, b, c, d, g = n[2], p, M, theta) +
    kappa[2] * alpha[3] * l1(a, b, c, d, g = n[3], p, M, theta) +
    kappa[2] * alpha[4] * l1(a, b, c, d, g = n[4], p, M, theta) +
    kappa[3] * beta[1] * l2(a, b, c, d, n[1], n[4], p, M, theta)+
    kappa[3] * beta[2] * l2(a, b, c, d, n[2], n[3], p, M, theta)
  lik
}
