#' Produce table comparing transformations based on exact calculations
#' 
#' Produce table
#'
#' @param method Character. 'MH','PR' or 'BA'
#' @param db Database
#' @param markers character vector
#' @param adjust logical 
#' @param rate Double
#' @param rate2 Double
#' @param range Double
#' 
#' @return Table (tab) and complete output (all)
#' 
#' @details None
#' 
#' @author Thore Egeland
#' 
#' @export
#' @examples
#' library(pedsuite)
#' res = tabExact(method = 'PR')$tab
#' 

tabExact = function(method = "PR", db = NorwegianFrequencies,
                      markers = c("TH01","D5S818","D13S317","D16S539","TPOX",
                                  "D10S1248","D22S1045","D3S1744","D2S1360",
                                  "D6S474","D4S2366","D5S2500","D10S2325"),
                    adjust = T, rate = 0.003, rate2 = 0.001, range = 0.5){
  n = length(markers)
  Min = Mu = Sd = Max = navn = fRatio = mindiag =  bM = bR = rep(NA, n)
  res = list()
  for ( i in 1:n){
    p = db[markers[i]]
    navn[i] = names(p)
    p = p[[1]]
    np = names(p)
    M = mutationMatrix("stepwise", alleles = np, rate = 0.001, rate2 = 0.0001, 
                       range = 0.1, afreq = p)
    R = findReversible(M,  method = method)
    if(!any(is.na(R))){
      if(adjust)
        R = adjustReversible(M, R)
      res[[i]] = exactLRR(M,R)
      Min[i] = res[[i]]$min
      Mu[i] = res[[i]]$muM
      Sd[i] = res[[i]]$sigma
      Max[i] = res[[i]]$max
      fRatio[i] = fratio(M,R)
      mindiag[i] = min(diag(R))
      bM[i] = isBounded(M)
      bR[i] = isBounded(R)
    }
  }
  tab = data.frame(min = Min, myM = Mu, sigmaM = Sd, max = Max, fRatio = fRatio, 
                   mindiag = mindiag, bM = bM, bR = bR)
  rownames(tab) = navn
  list(tab = tab, all = res)
}
