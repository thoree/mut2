#' Find reversible mutation matrix
#' 
#' The Metropolis - Hastings algorithm is used to convert 
#' a mutation matrix to a reversible one, i.e., find the
#' balanced matrix matrix.
#'
#' @param mutmat A mutation matrix.
#' @param method Character.
#' @param afreq A vector with allele frequencies 
#' of the same length as the size of mutmat.
#' @param check Logical.
#' 
#' @return Balanced mutation matrix. The expected mutation rate
#' of the balanced matrix is returned as `rate`.
#' 
#' @details Three different approaches are implemented.
#' The default, \code{method = "MH"}, gives a balanced matrix with off diagonal elements
#' 
#' \code{q_{ij} min(1, p_j/p_i * q_{ji}/q_{ij})}
#' 
#' where \code{q_{ij}} and \code{p_i} are the elements of the original mutation
#' matrix and the allele frequencies, respectively. 
#' The average method \code{method = "AV"}, gives a balanced matrix with off diagonal elements
#' 
#' \code{(p_i q_{ij} + p_j q_{ji} ) / (2p_i)}
#' 
#' if \code{q_{ji} < p_i, i neq j} (and may otherwise fail to balance).
#' 
#' The method Barker \code{method = "BA"}, gives a balanced matrix with off diagonal elements
#' 
#' \code{ p_j q_{ji}/(p_i q_{ij} + p_j q_{ji})}
#' 
#' If \code{q_{ij} = 0} or \code{q_{ji} = 0}, the elements (i,j) and (j,i) of the transformed
#' matrix becomes 0.
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
#' findReversible(mutmat)
#' findReversible(mutmat, method = "AV")
#' findReversible(mutmat, method = "BA")
#' 
#' Q = matrix(ncol = 2, c(0.9,0.9, 0.1, 0.1))
#' p = c(0.01, 0.99)
#' mutmat = mutationMatrix("custom", matrix = Q, alleles = 1:2)
#' findReversible(mutmat, afreq = p)
#' # AVerage balancing not possible:
#' # findReversible(mutmat, method = "AV", afreq = p)


findReversible = function(mutmat, method = "MH", afreq = NULL, check = TRUE){ 
  if(check)
    validateMutationMatrix(mutmat)
  if (is.null(afreq))
    afreq = attr(mutmat, "afreq")

  n = length(afreq)   
  P1 =  matrix(ncol = n, nrow = n, 0)
  
  if(method == "MH"){
    # Metropolis - Hastings. 
    # Also diagonal elements are calculated for simplicity,
    # but they are not used.
    for (i in 1:n)
     for (j in 1:n){
       if (afreq[i] * mutmat[i,j] > 0)
          P1[i,j]  = mutmat[i,j] * min(1,(afreq[j] * mutmat[j,i])/
                                         (afreq[i] * mutmat[i,j]))
     }
  }
  else if (method == "AV") {
    #Try average balancing
    for (i in 1:n)
      for (j in 1:n)
        P1[i,j]  = (afreq[i] * mutmat[i,j] + afreq[j] * mutmat[j,i])/
                 (2 * afreq[i])
  }
  else if (method == "BA") {
    #Barker balancing: P1[i,j] = 0 if q[i,j] or q[j,i] = 0
    for (i in 1:n)
      for (j in 1:n)
        if(afreq[j]*mutmat[j,i]*mutmat[i,j] > 0)
        P1[i,j]  = afreq[j]*mutmat[j,i]*mutmat[i,j]/
                   (afreq[i] * mutmat[i,j] + afreq[j] * mutmat[j,i])
  }
  else
    stop("Methods needs to be 'MH', 'AV' or 'BA'")
  
  # Find diagonal elements
  diag(P1) = 0
  s = apply(P1, 1, sum)
  diag(P1) = 1 - s
  
  dimnames(P1) = dimnames(mutmat)
  valid = validateMutationMatrix(P1)
  pedmut:::newMutationMatrix(P1, afreq = afreq, model = "custom",
                             rate = expectedMutationRate(P1, afreq))
}
