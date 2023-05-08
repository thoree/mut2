#' Checks if mutation model is well behaved
#' 
#' Checks if q_ij <= p_j
#' 
#'  
#' @param Q a mutation model
#' @param afreq A numeric vector of allele frequencies.
#' 
#' @return A logical
#' 
#' 
#' @author Thore Egeland.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' # Not well behaved:
#' p = c(0.4, 0.6)
#' Q = stepwiseReversible(alleles = 1:2, afreq = p, rate = 0.6, range = 0.1)
#' Q = Q$M
#' isWellBehaved(Q, afreq = p)
#' }

isWellBehaved = function(Q, afreq){
  n = length(afreq)
  Q = as.matrix(Q)
  lines = rep(NA, FALSE)
  for (i in 1:n)
    lines[i] = all(Q[i, -i] <= afreq[-i])
  all(lines)
}

