#' Reversible stepwise mutation model
#' 
#' A reversible stepwise mutation model is constructed.
#' 
#'  
#' @param alleles A character vector with allele labels.
#' @param afreq A numeric vector of allele frequencies.
#' @param rate  A numeric mutation rate.
#' @param range A positive number.
#' 
#' @return A reversible stepwise mutation model.
#' 
#' @details The approach of Dawid el al. (SJS, 2001) is implemented.
#' 
#' @author Thore Egeland.
#' 
#' @export
#' 
#' @examples
#' # Example 1. Small example
#' # Three alleles checked, equals stationary stepwise in R Familias
#' afreq = c(0.2, 0.3,  0.5)
#' alleles = 1:3
#' gamma = 0.05
#' delta = 0.1
#' R1 = stepwiseReversible(alleles = alleles, afreq = afreq, rate = gamma, range = delta)
#' library(Familias)
#' R2 = Familias:::FamiliasLocus(afreq, alleles, MutationModel = "Stepwise", 
#'              Stabilization = "PM", MutationRate = gamma, 
#'              MutationRange = delta)$maleMutationMatrix
#' max(abs(as.matrix(R1) - R2))
#' 
#' # Example 2 Larger example
#' n = 100
#' afreq = runif(n)
#' afreq =afreq/sum(afreq)
#' afreq = rep(1/n,n)
#' alleles = 1:n
#' gamma = min(afreq) 
#' gamma = 2*(n-1)*min(afreq)^2
#' delta = 0.1 
#' R1 = stepwiseReversible(alleles = alleles, afreq = afreq, rate = gamma, range = delta)
#' library(Familias)
#' R2 = Familias:::FamiliasLocus(afreq, alleles, MutationModel = "Stepwise", 
#'              Stabilization = "PM", MutationRate = gamma, 
#'              MutationRange = delta)$maleMutationMatrix
#' max(abs(as.matrix(R1) - R2))
#' \dontrun{
#' # Impossible parameter settings
#' p = c(0.01, 0.99)
#' stepwiseReversible(alleles = 1:2, afreq = p, rate = 0.05, range = 0.1)
#' # Not well behaved:
#' p = c(0.5, 0.5)
#' stepwiseReversible(alleles = 1:2, afreq = p, rate = 0.6, range = 0.1)
#' 
#' db = NorwegianFrequencies
#' m = db$D3S1358 # Only integers
#' alleles = as.integer(names(m))
#' afreq = as.double(m)
#' M = stepwiseReversible(alleles, afreq, rate = min(afreq), range = 0.1)
#' }

stepwiseReversible = function(alleles, afreq, rate, range){
  if(!is.integer(alleles))
    stop("alleles needs to be integers")
    # remaining checking will be taken care of by `mutationModel` below
    n = length(afreq)
    # Below the explicit formula for the transition probabilities is used with
    # rate = \gamma and range = \delta

    R = matrix(ncol = n, nrow = n, 0)

      for (i in 1:n){
        for(j in setdiff(1:n,i)){
            a = (1-range^n)/(1-range)
            R[i,j] = rate * (1 -range) * range^{abs(i-j)}/
                     (2*range*(n-a))*(1/afreq[i])
          }
        R[i,i] = 1 -sum(R[i,-i])   
      }
    
    dimnames(R) = list(alleles, alleles)
    R = mutationModel(matrix = R, model = "custom",  afreq = afreq)
    R
}
