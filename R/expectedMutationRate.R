#' Expected mutation rate of a mutation matrix
#' 
#' @param mutmat A mutation matrix.
#' @param afreq A vector with allele frequencies 
#' of the same length as the size of mutmat.
#' @param check Logical.
#' 
#' @return Expected mutation rate.
#' 
#' @author Thore Egeland
#' 
#' @export
#' 
#' @examples
#' library(pedmut)
#' n = 4
#' p= 1:n/sum(1:n)
#' names(p) = 1:n
#' mutmat = mutationMatrix("onestep", rate = 0.02, afreq = p, alleles = 1:n)
#' expectedMutationRate(mutmat, afreq = p)
#' 

expectedMutationRate = function(mutmat, afreq = NULL, check = TRUE) {
  if(check)
    validateMutationMatrix(mutmat)
  
  if (is.null(afreq))
    afreq = attr(mutmat, "afreq")
  if(is.null(afreq))
    stop("Allele freuencies needed")
  sum(afreq * (1-diag(mutmat)))
}
