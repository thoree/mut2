#' Pairwise LR bounds for proportional mutation models
#'
#' Reversibility is assumed for the mutation model and the LR is
#' calculated for a pair of non-inbred individuals comparing a kappa to unrelated.
#' @param p Vector of real numbers. Allele frequency vector.
#' @param kappa Vector of real numbers describing relationship.  IBD parameters
#'   for 0,1,2 IBD alleles.
#' @param K Real. Proportionality factor in proportional model
#' @param share Logical. If TRUE, both are homozygous for the same allele, 
#' otherwise no akllele sharing.
#' 
#' @return LR and the lower bound,  assuming no allele sharing, and 
#' the upper bound assuming the individuals to be homozygous for the rarest allele.
#' 
#' @author Thore Egeland <Thore.Egeland@nmbu.no>
#' @references Egeland, Pinto and Amorim, FSI:Genetics (2017),
#'   \doi{http://dx.doi.org/10.1016/j.fsigen.2017.04.018}.
#' @export
#' @importFrom expm %^%
#' @examples
#' 
#' # Parent offspring relationship LR
#' p = c("1" = 0.2, "2" = 0.8)
#' M = mutationMatrix("proportional", afreq = p, rate = 0.003)
#' gamma = mut2::expectedMutationRate(M, p)
#' K = gamma/(1- sum(p^2))
#' kappa = c(0, 1, 0)
#' LRmutBounds(p, kappa, K, share = FALSE)
#' LRmutBounds(p, kappa, K, share = TRUE)
LRmutBounds = function(p, kappa, K = 0, share = TRUE){
  if(share){
    lowerBound = 1
    pa = min(p)
    term = (1-K*(1-pa))/pa
    upperBound = kappa[1] + kappa[2]*term + kappa[3]*term^2
    text = "Share both"
  }
  else{
    lowerBound = kappa[1] + kappa[2]*K + kappa[3]*K^2
    upperBound = 1
    text = "No sharing"
  }
      
  data.frame("Case" = text, "lower" = lowerBound, "upper" = upperBound)
}

