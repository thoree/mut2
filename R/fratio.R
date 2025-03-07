#' Calculates fratio
#'
#' Calculates fratio, distance between matrices
#'
#' @param A Mutation matrix
#' @param B Mutation matrix
#' @return Double
#' @author Thore Egeland <Thore.Egeland@@nmbu.no>
#' @export
#' @importFrom expm %^%
#' @examples
#' library(pedsuite)
#' library(pedmut)
#' n = 2
#' p = c(0.2, 0.8)
#' names(p) = paste(1:n)
#' mutmat = mutationMatrix("equal", rate = 0.02, alleles = 1:2, afreq = p)
#' res = snpExample(mutmat, afreq = p, adjust = FALSE)
#' fratio(mutmat, res[[1]])
#' 
fratio = function(A,B){
  if(any(is.na(A)))
    return(NA)
  if(any(is.na(B)))
    return(NA)
  if(!inherits(A, what = "mutationMatrix") | !inherits(A, what = "mutationMatrix"))
    return(NA)
  # Assumes that a_{ij} > 0 iff b_{ij}>0
  if(min(A) == 0 | min(B) == 0){
    A = as.vector(A)
    A = A[A > 0]
    B = as.vector(B)
    B = B[B > 0]
    if(length(A) != length(B))
      return(NA)
    }
  max(max(A/B),max(B/A))
}

