#' Ratio of LRs are compared for different transformations
#' 
#' The ratio Z = LR(M,p)/LR(R,p) is calculated one marker.
#'  
#' @param M List of mutations matrices, typically not reversed
#' @param R List of reversed mutations matrices
#' @param ped A \code{ped} object
#' @param ids A numeric with ID labels of one or more pedigree members.
#' @param ln character  Logical.
#' 
#' @return The smallest value (min), the expected value and standard deviation wrt the
#' numerator and the largest value.
#' 
#' @details This gies summary statistics, the complete distribution for one marker
#' is provided by mut2::exactLRR.
#' 
#' 
#' @seealso [exactLRR()].
#' 
#' @author Thore Egeland.
#' 
#' @export
#' 
#' @examples
#' # Example One marker. PO case.
#' library(pedmut)
#' library(pedprobr)
#' ped = nuclearPed(1)
#' p = c("1" = 0.2, "2" = 0.8)
#' M = mutationMatrix("equal", afreq = p, alleles = 1:length(p), 
#'                    rate = 0.01)
#' R = findReversible(M, method = "PR", afreq = p)
#' M = list(M)
#' names(M) = "L1"
#' R = list(R)
#' res1 = compareLRR(M, R)
#' 
compareLRR = function(M, R, ped = nuclearPed(1), ids = c(1,3), ln = FALSE){
  n = length(M)
  if (n == 1){
    res = exactLRR(M[[1]], R[[1]],  adjust = F, ln = ln)[c("min","muM", "sigma", "max")]
    res = matrix(res, ncol = 1)
    dimnames(res) = list(c("min","muM", "sigma", "max"), names(M))
  }
  else{
    res = matrix(nrow = 4, ncol = n)
    for(i in 1:n){
      foo = exactLRR(M[[i]], R[[i]],  adjust = F, ln = ln)[c("min","muM", "sigma", "max")]
      res[,i] = unlist(foo)
    }
    dimnames(res) = list(c("min","muM", "sigma", "max"), names(M))
  }
  res
}