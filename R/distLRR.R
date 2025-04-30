#' Find summary statistics of LRR
#'
#' Finds distribution of LRR depending on allele frequencies.
#'
#'
#' @param adjust Logical.
#' @param method Transformation.
#' @param verbose Logical.
#' @param nsim integer.
#' @param seed Integer
#' @param log2 Logical.
#' @param mutMod character.
#' @param Rate double.
#' @param Rate2 Double.
#' @param range2 Double.
#' @param nAls  Integer.
#'
#'
#' @details The allele frequencies are simulated. The LRR distribution is then
#' exactly calculated using oneMarkerDistribution. RMSD (Root Mean Square
#' Deviation) is the square root of
#' E(LRR - E(LRR))^2, very close to SD(LRR) in our applications.
#'
#' @return Summary of distribution of LRR
#'
#' @author Thore Egeland
#'
#' @export
#'
#' @examples
#' distLRR()
#'



distLRR = function(adjust = TRUE, method = "PR", verbose = FALSE, nsim = 2,
               seed = NULL, log2 = FALSE, mutMod = "equal", Rate = c(0.005, 0.005),
               Rate2 = c(0.000001, 0.000001), Range = c(0.1, 0.1), nAls = 4){
  ids = c(1,3)
  ped1 = nuclearPed()
  transformed = TRUE # Changed to FALSE if PR transform fails
  set.seed(seed)

  Rate = runif(nsim, min = Rate[1], max = Rate[2])
  Rate2 = runif(nsim, min = Rate2[1], max = Rate2[2])
  Range = runif(nsim, min = Range[1], max = Range[2])

  dd = nAls*(nAls + 1)/2 # No of genotypes
  # To contain joint genotype probabilities under M and R:
  dM = dR = matrix(ncol = dd, nrow = dd, 0)

  for (i in 1:nsim){
    p = runif(nAls)
    p = p/sum(p)
    names(p) = 1:nAls
    M = mutationMatrix(model = mutMod, rate = Rate[i], rate2 = Rate2[i],
                       range = Range[i], afreq = p)

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


