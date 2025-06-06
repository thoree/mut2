#' Upper limit for overall mutation rate
#'  
#' @param alleles A character vector with allele labels.
#' @param afreq A numeric vector of allele frequencies.
#' @param range A positive number.
#' 
#' @return Bounds on gamma
#' 
#' @author Thore Egeland.
#' 
#' @export
#' 
#' @examples
#' alleles = 1:3; afreq = c(0.2, 0.3,  0.5);  range = 0.1
#' bounds = boundsGamma(alleles, afreq ,range)
#'                         
boundsGamma = function(alleles, afreq,  range){
  n = length(afreq)
  
  R1 = matrix(ncol = n, nrow = n, 0)
  
  for (i in 1:n){
    for(j in setdiff(1:n, i)){
      a = (1 - range^n)/(1 - range)
      R1[i,j] = 1 * (1 - range) * range^{abs(i-j)}/
        (2*range*(n - a))*(1/afreq[i])
    }
  }
  linesums = apply(R1, 1, sum)
  boundDefined = 1/max(linesums)
  maks = apply(R1, 2, max)
  c(maxGammaDefined = boundDefined, 
    gammaRegular = min(afreq/maks))
}