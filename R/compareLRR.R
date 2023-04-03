#' Ratio of LRs are compared for different tranformations
#' 
#' The ratio Z = LR(M,p)/LR(R,p) is calculated. for a SNP marker
#'  
#' @param M List of mutations matrices, typically not reversed
#' @param R List of reversed mutations matrices
#' @param ped A \code{ped} object
#' @param ids A numeric with ID labels 
#' of one or more pedigree members.
#' @param method Character specifying reversing method.
#' @param character adjust Logical.
#' @param character ln Logical.
#' 
#' @return Statistcs evaluation transformation
#' 
#' 
#' @seealso [makeReversible()].
#' 
#' @author Thore Egeland.
#' 
#' @export
#' 
#' @examples
#' ped = nuclearPed(1)
#' p = c("1" = 0.2, "2" = 0.8)
#' M = mutationMatrix("equal", afreq = p, alleles = 1:length(p), rate = 0.01)
#' R = makeReversible(M, method = "PM", afreq = p)
#' M = list(M)
#' names(M) = "L1"
#' R = list(R)
#' res1 = compareLRR(M,R)
#' 
#' M =  R = list()
#' M[[1]] = mutationMatrix("equal", afreq = p, alleles = 1:length(p), rate = 0.01)
#' M[[2]] = mutationMatrix("equal", afreq = p, alleles = 1:length(p), rate = 0.05)
#' R[[1]] = makeReversible(M[[1]], method = "PM", afreq = p)
#' R[[2]] = makeReversible(M[[2]], method = "PM", afreq = p)
#' names(M) = names(R) = paste0("L", 1:2)
#' res2 = compareLRR(M, R)
#' res3 = compareLRR(M, R, ln = TRUE)
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
