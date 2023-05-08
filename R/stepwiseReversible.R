#' Reversible stepwise mutation model
#' 
#' A reversible stepwise mutation model is constructed.
#' 
#'  
#' @param alleles A character vector with allele labels.
#' @param afreq A numeric vector of allele frequencies.
#' @param rate  A numeric mutation rate.
#' @param range A positive number.
#' 
#' @return A reversible stepwise mutation model.
#' 
#' @details The approach of Dawid el al. (SJS, 2001) is implemented.
#' 
#' @author Thore Egeland.
#' 
#' @export
#' 
#' @examples
#' p = c(0.2, 0.3,  0.5)
#' stepwiseReversible(alleles = 1:3, afreq = p, rate = 0.05, range = 0.1)
#' \dontrun{
#' # Impossible parameter settings
#' p = c(0.01, 0.99)
#' stepwiseReversible(alleles = 1:2, afreq = p, rate = 0.05, range = 0.1)
#' # Not well behaved:
#' p = c(0.5, 0.5)
#' stepwiseReversible(alleles = 1:2, afreq = p, rate = 0.6, range = 0.1)
#' 
#' db = NorwegianFrequencies
#' m = db$D3S1358 # Only integers
#' alleles = as.integer(names(m))
#' afreq = as.double(m)
#' M = stepwiseReversible(alleles, afreq, rate = min(afreq), range = 0.1)
#' }

stepwiseReversible = function(alleles, afreq, rate, range){
  if(!is.integer(alleles))
    stop("alleles needs to be integer")
    n = length(afreq)
    # Construct the matrix S used in the first part of the construction
    S = matrix(ncol = n, nrow =  n, 0)
    for (i in 1:n)
      for(j in setdiff(1:n,i))
        S[i,j] = range^{abs(i-j)}
    for (i in 1:n)
      S[i,i] = -sum(S[i,])
    diagVector = diag(S)
    lambda = rate/(-sum(diagVector))
    # Find largest possibe rate
    rateMax = min(-afreq/diagVector)*(-sum(diagVector))
    if(rateMax < rate)
      stop("Impossible parameter settings. Max rate = ", rateMax)
    
    Q = matrix(ncol = n, nrow = n, 0)
    for (i in 1:n)
      for(j in setdiff(1:n,i))
        Q[i,j] = lambda*S[i,j]/afreq[i]
    
    for (i in 1:n)
      Q[i,i] = 1 - sum(Q[i,-i])
    dimnames(Q) = list(alleles, alleles)
    gamma = expectedMutationRate(Q, afreq)
    Q = mutationModel(matrix = Q, model = "custom",  afreq = afreq)
    attr(Q, "rate") = gamma 
    list(M = Q, overallMutationRate = gamma, maxMutationRate = rateMax, )
}
