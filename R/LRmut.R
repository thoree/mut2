#' Pairwise LR with mutation
#'
#' Reversibility is assumed for the mutation model and the LR is
#' calculated for a pair of non-inbred individuals comparing a kappa to unrelated.
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
#'   case IBD=1, that the alleles are paternal-paternal, paternal-maternal,
#'   maternal-paternal, and maternal-maternal.
#' @param theta Real in `[0,1]`. Kinship coefficient.
#' @param K Real. Proportionality factor in proportional model
#' 
#' @return LR.
#' 
#' @author Thore Egeland <Thore.Egeland@nmbu.no>
#' @references Egeland, Pinto and Amorim, FSI:Genetics (2017),
#'   \doi{http://dx.doi.org/10.1016/j.fsigen.2017.04.018}.
#' @export
#' @importFrom expm %^%
#' @examples
#' 
#' # Parent offspring relationship LR
#' library(pedmut)
#' g1 = c(1,1)
#' g2 = c(2,2)
#' p = c("1" = 0.2, "2" = 0.8)
#' p = 1:10/sum(1:10)
#' names(p) = 1:10
#' M = mutationMatrix("proportional", afreq = p, rate = 0.003)
#' n = c(0, 1, 1, 0)
#' kappa = c(0, 1, 0)
#' alpha = c(0, 0.5, 0.5, 0)
#' gamma = mut2::expectedMutationRate(M, p)
#' K = gamma/(1- sum(p^2))
#' LR = LRmut(g1, g2, n, p, M, kappa = kappa, alpha, K = NULL)
#' LR - K # Difference implementation - exact
#' 
LRmut = function(g1, g2, n, p, M, kappa, alpha, theta = 0, K = 0){
  if(is.null(K))
    K = mut2::expectedMutationRate(M,p)/(1-sum(p^2))
  if(K > 1)
    stop("K > 1 not possible")
  likNumerator = lik2(g1 = g1, g2 = g2, n = n, p = p, M = M, 
                      kappa = kappa, alpha = alpha, theta = theta)
  a <- g1[1]
  b <- g1[2]
  c <- g2[1]
  d <- g2[2]
  likDenominator = mut2:::l0(a, b, c, d, p, theta)
  LR = likNumerator/likDenominator
  LR
}

