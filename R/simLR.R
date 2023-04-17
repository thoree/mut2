#' Simulate LRRR
#' 
#' The ratio Z = LR(M,p)/LR(R,p) is simulated where 
#' M is a list of mutation matrices and R a reversed version
#'  
#' @param M List of mutation matrices, one for each marker.
#' @param R List of reversed mutation matrices, one for each marker.
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
#' simLRR(M, markerNames = mN, nsim = 10, seed = 177)
#' }


simLR = function(M, R, ped1 = nuclearPed(1), ids = c(1,3), 
                  nsim = 2, seed = NULL){
  set.seed(seed)
  n = dim(M)[1]
  unr = list(singleton(ids[1]), singleton(ids[2]))
    for (i in 1:n){
      ped1 = addMarker(ped1, mutmod = M[[i]], afreq = attr(M[[i]], "afreq"), name = name(M[i]))
      ped2 = addMarker(ped2, mutmod = R[[i]], afreq = attr(R[[i]], "afreq"), name = name(R[i]))
    }
  # Simulations under numerator
  simM = profileSim(ped1, N = nsim, ids = ids)
  simR = lapply(simM, function(x) transferMarkers(x, ped2, erase = FALSE))
  likM = unlist(lapply(simM, function(x) prod(likelihood(x))))
  likR = unlist(lapply(simR, function(x) prod(likelihood(x))))
  LRR.M = likM/likR
  
  unr.M = lapply(simM, function(x) transferMarkers(x, unr))
  likunr.M = unlist(lapply(unr.M, function(x) prod(likelihood(x))))
  LR.M = likM/likunr.M
  
  # Simulations under denominator
  simR = profileSim(ped2, N = nsim, ids = ids)
  simM = lapply(simR, function(x) transferMarkers(x, ped1, erase = FALSE))
  likM = unlist(lapply(simM, function(x) prod(likelihood(x))))
  likR = unlist(lapply(simR, function(x) prod(likelihood(x))))
  LRR.R = likM/likR
  
  unr.R = lapply(simR, function(x) transferMarkers(x, unr))
  likunr.R = unlist(lapply(unr.R, function(x) prod(likelihood(x))))
  LR.R= likR/likunr.M
  data.frame(LRR.M = LRR.M, LRR.R = LRR.R, LR.M = LR.M, LR.R = LR.R)
}

