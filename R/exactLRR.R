#' Ratio of LRs and likelihood exact for one marker
#' 
#' The ratio Z = LR(M,p)/LR(R,p) is calculated exactly.
#'  
#' @param M Mutation matrix
#' @param R Mutation matrix with same dimension as M or NULL
#' @param afreq A vector with allele frequencies, 
#' of the same length as the size of mutmat
#' @param ped A \code{ped} object
#' @param ids A numeric with ID labels 
#' of one or more pedigree members.
#' @param method Character specifying reversing method.
#' @param character adjust Logical.
#' @param character ln Logical.
#' 
#' @return Expected ratio and likelihoods
#' 
#' @details If \code{R == NULL}, \code{R} is first made reversible with
#' \code{PM} option, preserving the expected mutation rate.
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
#' exactLRR(M, R, ln = FALSE)
#' 
#' \dontrun{
#' p = NorwegianFrequencies[[1]]
#' M = mutationMatrix("onestep",
#'     alleles = names(p),
#'     rate = 0.001,
#'     afreq = p)
#' attr(M, "rate") = expectedMutationRate(M, p)
#' res = exactLRR(M, R = NULL, ln = T)
#' R = makeReversible(M, afreq = p, method = "MH")
#' R2 = adjustReversible(M, R, afreq = p)
#' res2 = exactLRR(M, R = R, p)
#' }

exactLRR = function(M, R = NULL, afreq = NULL, ped = nuclearPed(1), ids = c(1,3),
                    method = "PM", adjust = FALSE, ln = F){
  if(is.null(afreq))
    afreq = attr(M, "afreq")
  m = marker(ped, afreq = afreq)
  x1 = setMarkers(ped,m)
  x1 = setMutationModel(x1, "custom", matrix = M)
  likM = oneMarkerDistribution(x1, partialmarker = 1, 
                              ids = ids, verbose = F)

  if(is.null(R))
    R = makeReversible(M, method = method, afreq = afreq)
  if(adjust)
    R = adjustReversible(M, R)
  
  x2 = setMarkers(ped, m)
  x2 = setMutationModel(x2, "custom", matrix = R)
  likR = oneMarkerDistribution(x2, partialmarker = 1, 
                              ids = ids, verbose = F)
  Z = likM/likR
  if(ln)
    Z[!is.na(Z)] = log(Z[!is.na(Z)])
  EZM = sum(likM*Z, na.rm = T)
  EZR2 = sum(likR*Z^2, na.rm = T)
  list(Z = Z, likM = likM, likR = likR, 
       muM = EZM, sigma = sd(Z, na.rm = T), 
       min = min(Z, na.rm = T), 
       max = max(Z, na.rm = T), muR2 = EZR2)
}
