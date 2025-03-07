#' Adjusts Metropolis - Hasting balanced mutation matrix
#' 
#' The Metropolis - Hastings conversion typically gives a mutation matrix with too small
#' small mutation rate. The previously balanced mutation matrix is adjusted
#' to have the same mutation rate as the original matrix. 
#'  
#' @param mutmat Original mutation matrix.
#' @param balancedMutmat Balanced, mutation matrix. If, NULL, balanced.
#' @param method Character.  'MH','PR' or 'BA'.
#' @param afreq A vector with allele frequencies. 
#' of the same length as the size of mutmat.
#' @param check Logical.
#' 
#' @return Adjusted mutation matrix
#' 
#' @details If \code{balancedMutmat == NULL}, \code{mutmat} is first balanced.
#' The adjusted balanced matrix is
#' 
#' \code{alpha * balancedMutmat + (1-alpha) * I}
#' 
#' where
#' 
#' \code{alpha} is the ratio of the (expected mutation) rates of the original matrix,
#' \code{mutmat} to the balanced version \code{balancedMutmat} and \code{I} is the identity matrix.
#'  
#' @seealso [findReveversible()].
#' 
#' @author Thore Egeland.
#' 
#' @export
#' 
#' @examples
#' library(pedmut)
#' afreq = c(0.1, 0.3, 0.4, 0.2)
#' names(afreq) = 1:4
#' mutmat = mutationMatrix("onestep", alleles = 1:4, rate = 0.02)
#' balancedMutmatMH = findReversible(mutmat, afreq = afreq, method = "MH")
#' adj = adjustReversible(mutmat, balancedMutmatMH, afreq = afreq,  check = TRUE)
#' attr(mutmat, "rate") - attr(adj, "rate")
#' 
#' # Does 'PR' give adjusted matrix directly?:
#' balancedMutmatPR = findReversible(mutmat, afreq = afreq, method = "PR")

adjustReversible = function(mutmat, balancedMutmat, method = "MH",
                            afreq = NULL,  check = TRUE){ 

  if (is.null(afreq))
    afreq = attr(mutmat, "afreq")

  if(is.null(balancedMutmat))
    balancedMutmat = findReversible(mutmat, method = method, 
                                     afreq = afreq, check = check)
    
  if(check){
    if(!isReversible(balancedMutmat, afreq))
      stop("Second argument needs to be a balanced mutation matrix")
  }
  
  RM = expectedMutationRate(mutmat, afreq)
  RP = expectedMutationRate(balancedMutmat, afreq)
  if(is.null(RM) | is.null(RP) | (RP == 0))
    alpha = 1
  else
    alpha = RM/RP
  newM = alpha * balancedMutmat + (1-alpha) * diag(rep(1,length(afreq)))
  pedmut:::newMutationMatrix(newM, afreq = afreq, model = "custom",
                             rate = expectedMutationRate(newM, afreq))
  }

