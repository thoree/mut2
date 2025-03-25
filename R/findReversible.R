#' Find reversible mutation matrix
#'
#' The Metropolis - Hastings algorithm is used to convert
#' a mutation matrix to a reversible one.
#'
#' @param mutmat A mutation matrix.
#' @param method Character. 'MH','PR' or 'BA'
#' @param afreq A vector with allele frequencies
#' of the same length as the size of mutmat.
#' @param check Logical.
#'
#' @return Reversible mutation matrix. The expected mutation rate
#' of the balanced matrix is returned as `rate`.
#'
#' @details Three different approaches are implemented.
#' The \code{method = "MH"}, based on the traditional proposal of MCMC,
#' gives a balanced matrix with off diagonal elements
#'
#' \code{q_{ij} min(1, p_j/p_i * q_{ji}/q_{ij})}
#'
#' where \code{q_{ij}} and \code{p_i} are the elements of the original mutation
#' matrix and the allele frequencies, respectively.
#'
#' The method \code{method = "BA"}, gives a balanced matrix with off diagonal elements
#'
#' \code{p_j q_{ij}q_{ji} / (p_i q_{ij} + (p_j q_{ji})}.
#'
#' The  method \code{method = "PR"} is given by
#'
#' \code{(p_i q_{ij} + p_j q_{ji} ) / (2p_i)}
#'
#' if \code{q_{ji} < p_i, i neq j} (and may otherwise fail to balance).
#' The expected mutation rate is preserved.
#'

#' @author Thore Egeland
#'
#' @export
#'
#' @examples
#' library(pedmut)
#' n = 4
#' p = 1:n/sum(1:n)
#' names(p) = 1:n
#' mutmat = mutationMatrix("onestep", rate = 0.02, afreq = p, alleles = 1:n)
#' R1 = findReversible(mutmat, method = "MH")
#' expectedMutationRate(R1)
#' isReversible(R1)
#' R2 = findReversible(mutmat, method = "PR")
#' isReversible(R2)
#'
#' R3 = findReversible(mutmat, method = "BA")
#' isReversible(R3)

findReversible = function(mutmat, method = "BA", afreq = NULL, check = TRUE){
  if(check)
    validateMutationMatrix(mutmat)
  if (is.null(afreq))
    afreq = attr(mutmat, "afreq")
  if(is.null(afreq))
    stop("No allele frequencies provided.")
  n = length(afreq)
  P1 =  matrix(ncol = n, nrow = n, 0)

  if(method == "MH"){
    # Metropolis - Hastings.
    for (i in 1:n)
     for (j in (1:n)[-i]){
       if (mutmat[i,j] > 0)
          P1[i,j]  = mutmat[i,j] * min(1,(afreq[j] * mutmat[j,i])/
                                         (afreq[i] * mutmat[i,j]))
     }
  }
  else if (method == "PR") {
    # PR Preserved expected Mutation rate method
    for (i in 1:n)
      for (j in (1:n)[-i])
        P1[i,j]  = (afreq[i] * mutmat[i,j] + afreq[j] * mutmat[j,i])/
                 (2 * afreq[i])
    }
  else if (method == "BA")
    # BArker transformation
    for (i in 1:n)
      for (j in (1:n)[-i]){
        if (mutmat[i,j] > 0 | mutmat[j,i] > 0)
          P1[i,j]  = afreq[j] * (mutmat[i,j] * mutmat[j,i])/
                                (afreq[i] * mutmat[i,j] +afreq[j] * mutmat[j,i])
        }
  else
    stop("Methods needs to be 'MH','PR' or 'BA'")

  # Find diagonal elements
  diag(P1) = 0
  s = apply(P1, 1, sum)
  diag(P1) = 1 - s

  dimnames(P1) = dimnames(mutmat)
  if(any(as.matrix(P1) < 0))
    P1 = NA
  else
    P1 =  pedmut:::newMutationMatrix(P1, afreq = afreq, model = "custom",
                             rate = expectedMutationRate(P1, afreq))
  P1
}
