#' Reversible stepwise mutation model
#'  
#' @param alleles A character vector with allele labels.
#' @param afreq A numeric vector of allele frequencies.
#' @param rate  A numeric mutation rate.
#' @param range A positive number.
#' 
#' @return A reversible stepwise mutation model with expected mutation rate equal input rate.
#' 
#' @details The approach of Dawid el al. (SJS, 2001) is implemented.
#' 
#' @author Thore Egeland.
#' 
#' @export
#' 
#' @examples
#' R = stepwiseReversible(alleles = 1:3, afreq = c(0.2, 0.3,  0.5), rate = 0.001, range = 0.1)
#' 
stepwiseReversible = function(alleles, afreq, rate, range){
  if(!is.integer(alleles))
    stop("alleles needs to be integers")
  # remaining checking will be taken care of by `mutationModel` below
  n = length(afreq)

  if(mut2::boundsGamma(alleles, afreq,  range)[[1]] < rate)
    return("Model not well defined")
  else{
  R = matrix(ncol = n, nrow = n, 0)

      for (i in 1:n){
        for(j in setdiff(1:n, i)){
            a = (1 - range^n)/(1 - range)
            R[i,j] = rate * (1 - range) * range^{abs(i-j)}/
                     (2*range*(n - a))*(1/afreq[i])
          }
        R[i,i] = 1 - sum(R[i,-i])   
      }
    
    dimnames(R) = list(alleles, alleles)
    R = mutationModel(matrix = R, model = "custom",  afreq = afreq, alleles = alleles)
    R
  }
}
