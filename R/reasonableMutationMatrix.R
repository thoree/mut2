#' Checks if pj < mij
#' 
#' Entries where  pj < mij are unreasonable and are detected
#' @param mutmat A mutation matrix.
#' @param afreq A vector with allele frequencies 
#' of the same length as the size of mutmat.
#' @param verbose Logical.
#' 
#' @return Number of cases where pj < mij
#' 
#' @author Thore Egeland
#' 
#' @export
#' 
#' @examples
#' n = 4
#' p = 1:n/sum(1:n)
#' names(p) = 1:n
#' mutmat = mutationMatrix("equal", afreq = p, rate = 0.5)
#' reasonableMutationMatrix(mutmat, p, verbose = F)
#' 

reasonableMutationMatrix = function(mutmat, p, verbose = FALSE){

  if(is.null(p))
    p = attr(mutmat, "afreq")
  n = length(p)
  problem = matrix(ncol = n, nrow = n, FALSE)
  dimnames(problem) = list(names(p), names(p))
  for (i in 1:n){
    for(j in setdiff(1:n,i)){
      if (p[j] < mutmat[i,j])
          problem[i,j] = TRUE
    }
  }
  sump = sum(problem) 
  if(verbose){
    if(sump == 0)
      cat("No problems detected: mij <= pj")
    else{
      cat("TRUE indicates entries with pi < mij", "\n")
      print(problem)
      cat("\n")
    }
  }
  sump
}
