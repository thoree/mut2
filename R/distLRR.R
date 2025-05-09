#' Find summary statistics of LRR
#'
#' Finds distribution of LRR depending on allele frequencies.
#'
#'
#' @param adjust Logical. Transform to original mutation rate
#' @param method Transformation. See `pedmut::makeReversible()`
#' @param verbose Logical.
#' @param nsim integer.
#' @param seed Integer
#' @param log2 Logical. To get LR on log2
#' @param mutMod character.
#' @param rate Double. See `pedmut::muationModel()`.
#' @param rate2 Double. See `pedmut::muationModel()`
#' @param range Double. See `pedmut::muationModel()`
#' @param nAls  Integer. No of alleles.
#'
#'
#' @details The allele frequencies are simulated uniformly. The mutation
#' parameters are simulated uniformly between provided minimum and maximum.
#' The LRR distribution is then calculated exactlyusing
#' `oneMarkerDistribution()`. The simulations are conditional on PR not failing
#' (PR fails rarely).#'
#'
#' @return Summary of distribution of LRR, with more details if `verbose = T`
#'
#' @author Thore Egeland
#'
#' @export
#'
#' @examples
#' distLRR()
#'



distLRR = function(adjust = TRUE, method = "PR", verbose = FALSE, nsim = 2,
               seed = NULL, log2 = FALSE, mutMod = "equal", rate = c(0.005, 0.005),
               rate2 = c(0.000001, 0.000001), range = c(0.1, 0.1), nAls = 4){
  ids = c(1,3)
  ped1 = nuclearPed()
  transformed = TRUE # Changed to FALSE if PR transform fails
  set.seed(seed)

  rate = runif(nsim, min = rate[1], max = rate[2])
  rate2 = runif(nsim, min = rate2[1], max = rate2[2])
  range = runif(nsim, min = range[1], max = range[2])

  dd = nAls*(nAls + 1)/2 # No of genotypes
  # To contain joint genotype probabilities under M and R:
  dM = dR = matrix(ncol = dd, nrow = dd, 0)

  for (i in 1:nsim){
    p = runif(nAls)
    p = p/sum(p)
    names(p) = 1:nAls
    M = mutationMatrix(model = mutMod, rate = rate[i], rate2 = rate2[i],
                       range = range[i], afreq = p)

    if(method == "PR") {
      pm = M*p
      if (all(colSums(pm) <= p * (1 + 2 * diag(M))))
        R = makeReversible(M, method = "PR", adjust = F)
      else
        transformed = FALSE
    } else
      R = makeReversible(M, method = method, adjust = adjust)

    if(transformed){
      m = marker(ped1, afreq = p)
      ped1 = setMarkers(ped1, m)
      ped1 = setMutmod(ped1, model =  "custom", matrix = M)
      oneM = oneMarkerDistribution(ped1, ids, marker = 1, verbose = F)
      dM = dM + oneM
      ped1 = setMutmod(ped1, model =  "custom", matrix = R)
      oneR = oneMarkerDistribution(ped1, ids, marker = 1, verbose = F)
      dR = dR + oneR
    }
  }
  distM = dM/nsim
  distR = dR/nsim
  if(log2)
    LRR = log2(distM/distR)
  else
    LRR = distM/distR
  mu = sum(LRR*distM)
  mu2 = sum(LRR^2*distM)
  sd = sqrt(mu2-mu^2)
  index = which(LRR == min(LRR), arr.ind = TRUE)
  Min = unique(LRR[index])
  pMin = sum(distM[index])
  index = which(LRR == max(LRR), arr.ind = TRUE)
  Max = unique(LRR[index])
  pMax = sum(distM[index])
  sum2 = data.frame(mutMod = mutMod, transform = method, mu = mu, sd = sd, min = Min, pMin = pMin, max = Max, pMax = pMax)
  if(!verbose)
    sum2
  else
    list(summary = sum2,  LRR = LRR, distM = dM/nsim, distR = dR/nsim)
}


