#' Make mutation model well behaved
#' 
#' Replace q_ij by min(q_ij,p_j)
#' 
#'  
#' @param Q a mutation model
#' @param afreq A numeric vector of allele frequencies.
#' 
#' @return A reversible stepwise mutation model.
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
#' makeWellBehaved(Q, afreq = p)
#' }

makeWellBehaved = function(Q, afreq){
n = length(afreq)
Q = as.matrix(Q)
for (i in 1:n)
  for(j in setdiff(1:n,i))
      Q[i,j] = min(Q[i,j], afreq[j])
for (i in 1:n)
  Q[i,i] = 1 - sum(Q[i,-i])

gamma = expectedMutationRate(Q, afreq)
Q = mutationModel(matrix = Q, model = "custom",  afreq = p)
list(Q = Q,  gamma = gamma)
}
