#' Expected ratio of LRs
#' 
#' The ratio Z = LR(Q,p)/LR(M,p) and its 
#' expected value is calculated.
#'  
#' @param Q Mutation matrix
#' @param M Mutation matrix with same dimension as Q
#' @param afreq A vector with allele frequencies, 
#' of the same length as the size of mutmat
#' @param ped A \code{ped} object
#' @param ids A numeric with ID labels 
#' of one or more pedigree members.
#' 
#' @return Expected ratio and ratios
#' 
#' @details If \code{M == NULL}, \code{M} is first made reversible with
#' \code{PM} option, preserving the expected mutation rate.
#' 
#' @seealso [findReversible()].
#' 
#' @author Thore Egeland.
#' 
#' @export
#' 
#' @examples
#' ids = c("id1", "id2")
#' ped = nuclearPed(father = ids[1], children = ids[2])
#' p = c("1" = 0.2, "2" = 0.8)
#' Q = mutationMatrix("equal", afreq = p, alleles = 1:length(p), rate = 0.01)
#' Q = mutationMatrix("onestep", afreq = p, alleles = 1:length(p), rate = 0.01)
#' M = findReversible(Q, method = "PM", afreq = p)
#' LRR(Q, M, p, ped, ids)


LRR = function(Q, M, afreq, ped, ids){
  m = marker(ped, afreq = afreq)
  x1 = setMarkers(ped, m)
  x1 = setMutationModel(x1, "custom", matrix = Q)
  likQ = oneMarkerDistribution(x1, partialmarker = 1, 
                              ids = ids, verbose = F)
  x2 = setMarkers(ped, m)
  if(is.null(M))
    M = findReversible(Q, method = "PM", afreq = afreq)
  x2 = setMutationModel(x2, "custom", matrix = M)
  likM = oneMarkerDistribution(x2, partialmarker = 1, 
                              ids = ids, verbose = F)
  Z = likQ/likM
  EZQ = sum(likQ*Z, na.rm = T)
  list(EZQ = EZQ, Z = Z)
}