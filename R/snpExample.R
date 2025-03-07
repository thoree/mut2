#' SNP example
#' 
#' The three reversible transformations and exact
#'
#' @param mutmat A mutation matrix.
#' @param afreq A vector with allele frequencies 
#' @param check Logical
#' @param adjust Logical
#' 
#' @return Balanced mutation matrix. The expected mutation rate
#' of the balanced matrix is returned as `rate`.
#' 
#' @author Thore Egeland
#' 
#' @export
#' 
#' @examples
#' library(pedmut)
#' n = 2
#' p = c("1" = 0.2, "2" = 0.8)#' 
#' m = rbind(c(0.997, 0.003), c(0.003, 0.997))
#' mutmat = mutationMatrix("custom", matrix  = m, afreq = p, alleles = 1:2)
#' res = snpExample(mutmat, afreq = p, adjust = FALSE)
#' 
snpExample = function(mutmat, afreq = NULL, check = TRUE, adjust = TRUE){ 
  if(check)
    validateMutationMatrix(mutmat)
  if (is.null(afreq))
    afreq = attr(mutmat, "afreq")
  if(is.null(afreq))
    stop("No allele frequencies provided.")
  n = length(afreq) 
  if( n!= 2)
    stop("Only for SNP-s")
  gammaM = expectedMutationRate(mutmat, afreq = afreq)
  reversed = list()
  reversed[[1]] = findReversible(mutmat, method = "MH", afreq = afreq)
  reversed[[2]] = findReversible(mutmat, method = "PR", afreq = afreq)
  reversed[[3]] = findReversible(mutmat, method = "BA", afreq = afreq)
  if(adjust){
    reversed[[1]]  = adjustReversible(mutmat, reversed[[1]], afreq = afreq)
    reversed[[3]]  = adjustReversible(mutmat, reversed[[3]], afreq = afreq)
  }
  P1 = matrix(c(1 - gammaM/(2*afreq[1]), gammaM/(2*afreq[1]),
              gammaM/(2*afreq[2]), c(1 - gammaM/(2*afreq[2]))),
              ncol = 2, nrow = 2, byrow = T)

  reversed[[4]] = pedmut:::newMutationMatrix(P1, afreq = afreq, model = "custom",
                             rate = expectedMutationRate(P1, afreq))
  names(reversed) = c("MH", "PR", "BA", "EQ")
  reversed
  }
