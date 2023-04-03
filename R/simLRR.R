#' Simulate LRRR
#' 
#' The ratio Z = LR(M,p)/LR(R,p) is simulated where 
#' M is a list of mutation matrices and R a reversed version
#'  
#' @param M List of mutation matrices, one for each marker.
#' @param R List of reversed mutation matrices, one for each marker.
#' @param markerNames Character vector, names of markers
#' @param method Character, reversing method.
#' @param ped1 A \code{ped} object.
#' @param ids A numeric with ID labels of one or more pedigree members.
#' @param nsim Integer.     
#' @param seed Integer.
#' 
#' @details
#' If \code{R == NULL}, mutation matrices are reversed.
#' 
#' @return LRR
#' 
#' @seealso [makeReversible()].
#' 
#' @author Thore Egeland.
#' 
#' @export
#' 
#' @examples
#' 
#' p = c("1" = 0.1, "2" = 0.9)
#' M = list(mutationMatrix("equal", rate = 0.001, afreq = p),
#'          mutationMatrix("equal", rate = 0.001, afreq = p))
#' mN = c("L1", "L2")
#' markerNames = mN; method = "PM"; ped1 = nuclearPed(1); ids = c(1,3); nsim = 2; seed = 17
#' res = simLRR(M, markerNames = mN, nsim = 100, seed = 177)
#' simLRR(M, markerNames = mN, nsim = 10, seed = 177)
#' res = simLRR(M[[1]], markerNames = mN[1], nsim = 10, seed = 177)
#' \dontrun{
#' M = lapply(NorwegianFrequencies, function(x){
#'   mat = mutationMatrix("stepwise",
#'   alleles = names(x),
#'   rate = 0.001, 
#'   rate2 = 1e-05,
#'   range = 0.1,
#'   afreq = x)
#'   mat
#' })
#' mN = names(NorwegianFrequencies)
#' simLRR(M, markerNames = mN, nsim = 10, seed = 177)
#' }


simLRR = function(M, R = NULL, markerNames, method = "PM", ped1 = nuclearPed(1), ids = c(1,3), 
                  nsim = 2, seed = NULL){
  set.seed(seed)
  # Find reversible matrices
  mN = markerNames
  n = length(mN)
  if(is.null(R)) {
    if(n == 1)
      R = makeReversible(M, method = method)
    else
      R = lapply(M, function(x) makeReversible(x, method = method))
  }

  
  #Add markers with mutation matrices M
  ped2 = ped1
  if (n == 1){
        ped1 = addMarker(ped1, mutmod = M, afreq = attr(M, "afreq"), name = mN)
        ped2 = addMarker(ped2, mutmod = R, afreq = attr(R, "afreq"), name = mN)
  } else{
    for (i in 1:n){
      ped1 = addMarker(ped1, mutmod = M[[i]], afreq = attr(M[[i]], "afreq"), name = mN[i])
      ped2 = addMarker(ped2, mutmod = R[[i]], afreq = attr(R[[i]], "afreq"), name = mN[i])
    }
  }
  
  # Simulations under numerator
  simM = profileSim(ped1, N = nsim, ids = ids)
  simR = lapply(simM, function(x) transferMarkers(x, ped2, erase = FALSE))
  likM = unlist(lapply(simM, function(x) prod(likelihood(x))))
  likR = unlist(lapply(simR, function(x) prod(likelihood(x))))
  LRR.M = likM/likR
  
  # Simulations under denominator
  simR = profileSim(ped2, N = nsim, ids = ids)
  simM = lapply(simR, function(x) transferMarkers(x, ped1, erase = FALSE))
  likM = unlist(lapply(simM, function(x) prod(likelihood(x))))
  likR = unlist(lapply(simR, function(x) prod(likelihood(x))))
  LRR.R = likM/likR

  data.frame(LRR.M = LRR.M, LRR.R = LRR.R)
}

