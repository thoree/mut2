#' Produce table comparing transformations based on simulations
#'
#' Produce table
#'
#' @param db allele frequencies
#' @param rate Double
#' @param rate2 Double
#' @param range Double
#' @param nsim  Integer
#' @param seed  Integer
#'
#' @return LRR table
#'
#' @details LRR table
#'
#' @author Thore Egeland
#'
#' @export
#' @examples
#' db = getFreqDatabase(KLINK::halfsib[[1]])
#' tabLRR(db = db)


tabLRR = function(db = NULL, rate = 0.001, rate2 = 0.00001,
                  range = 0.1, nsim = 2, seed = 1729, maks = FALSE){
  if(maks)
    ind = 3:4
  else
    ind = 1:2
  # Equal model

  M = lapply(db, function(x){
    mat = mutationMatrix("equal",
                         alleles = names(x),
                         rate = rate,
                         rate2 = NULL,
                         range = NULL,
                         afreq = x)
    mat
  })
  resMH = simLRR(M, markerNames = names(db), nsim =  nsim, seed = seed, method = "MH")[,ind]
  resBA = simLRR(M, markerNames = names(db), nsim =  nsim, seed = seed, method = "BA")[,ind]
  resPR = simLRR(M, markerNames = names(db), nsim =  nsim, seed = seed, method = "PR")[,ind]
  res1 = cbind(resMH, resBA, resPR)
  colnames(res1) = c("eq.MH.num", "eq.MH.den", "eq.BA.num", "eq.BA.den", "eq.PR.num", "eq.PR.den")
  # stepwise model
  M = lapply(db, function(x){
    mat = mutationMatrix("stepwise",
                         alleles = names(x),
                         rate = rate,
                         rate2 = rate2,
                         range = range,
                         afreq = x)
    mat
  })
  resMH = simLRR(M, markerNames = names(db), nsim =  nsim, seed = seed, method = "MH")[,ind]
  resBA = simLRR(M, markerNames = names(db), nsim =  nsim, seed = seed, method = "BA")[,ind]
  resPR = simLRR(M, markerNames = names(db), nsim =  nsim, seed = seed, method = "PR")[,ind]
  res2 = cbind(resMH, resBA, resPR)
  colnames(res2) = c("st.MH.num", "st.MH.den", "st.BA.num", "st.BA.den", "st.PR.num", "st.PR.den")
  cbind(res1, res2)
}

