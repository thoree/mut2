#' Produce table comparing transformations
#'
#' Produce table
#'
#' @param db allele frequencies
#' @param rate Double
#' @param rate2 Double
#' @param range Double
#' @param mutmodel  Character
#' @param relabel Logical
#' @param stationary Logical
#' @param nr Integer 0 complete, 1,2,3 gives Table 1,2,3
#'
#' @return Table (tab) and complete output (all)
#'
#' @details If logical is TRUE, a stationary version is found
#'
#' @author Thore Egeland
#'
#' @export
#' @examples
#' library(pedtools)
#' library(pedmut)
#' library(forrel)
#' tabfRatio(db = NorwegianFrequencies[1:2], rate = 0.0001, mutmodel = 'equal')


tabfRatio = function(db = NULL, rate = 0.001, rate2 = 0.00001, range = 0.1,
                    mutmodel = "stepwise", relabel = F, stationary = T,
                    nr = 0){
  n = length(db)
  fRatioMH =  fRatioBA =  fRatioPR = fRatioDA = fRatioStationary = alsorev = rep(NA, n)
  mindiagMH =  mindiagBA =  mindiagPR = mindiagDA= rep(NA, n)
  bM = bMH =  bBA =  bPR = bDA = rep(NA, n)
  navn = rep(NA,n)
  integerMarker = rep(NA, n)
  for ( i in 1:n){
    p = db[i]
    navn[i] = names(p)
    p = p[[1]]
    np = names(p)
    p = as.double(p)
    dnp = as.double(np)
    integerMarker[i] = all(floor(dnp) == dnp)
    if(integerMarker[i]) #relabel to avoid gaps
      np = 1:length(np)
    if(relabel) #always relabel
      np = 1:length(np)
    if(mutmodel == 'stepwise')
      M = mutationMatrix(mutmodel, alleles = np, rate = rate, rate2 = rate2,
                       range = range, afreq = p)
    else
      M = mutationMatrix(mutmodel, alleles = np, rate = rate, afreq = p)

    MH = findReversible(M,  method = 'MH')
    MH = adjustReversible(M, MH)

    BA = findReversible(M,  method = 'BA')
    BA = adjustReversible(M, BA)

    PR = findReversible(M,  method = 'PR')
    if (integerMarker[i])
     DA = mut2::stepwiseReversible(alleles = as.integer(np), afreq = p, rate = rate, range = range)
    else
     DA = NA

    fRatioMH[i] = fratio(M, MH)
    fRatioBA[i] = fratio(M, BA)
    fRatioPR[i] = fratio(M, PR)
    fRatioDA[i] = fratio(M, DA)

    if(stationary){
     MSTAT = stabilize(M)
     fRatioStationary[i] = fratio(M, MSTAT)
     alsorev[i] = isReversible(MSTAT)
    }

    mindiagMH[i] = min(diag(MH))
    mindiagBA[i] = min(diag(BA))
    mindiagPR[i] = ifelse(class(PR)[1] == 'mutationMatrix', min(diag(PR)), NA)
    mindiagDA[i] = ifelse(class(DA)[1] == 'mutationMatrix', min(diag(DA)), NA)

    bM[i] = isBounded(M)
    bMH[i] = isBounded(MH)
    bBA[i] = isBounded(BA)
    bPR[i] = ifelse(class(PR)[1] == 'mutationMatrix', isBounded(PR), NA)
    bDA[i] = ifelse(class(DA)[1] == 'mutationMatrix', isBounded(DA), NA)
  }

  bM = ifelse(bM,"Y", "N")
  bMH = ifelse(bMH,"Y", "N")
  bBA = ifelse(bBA,"Y", "N")
  bPR = ifelse(bPR,"Y", "N")
  bDA = ifelse(bDA,"Y", "N")
  integerMarker = ifelse(integerMarker,"Y", "N")

  tab = data.frame(fMH =  fRatioMH, fBA = fRatioBA, fPR =  fRatioPR, fDA =  fRatioDA,
                   mMH = mindiagMH, mBA = mindiagBA, mPR = mindiagPR, mDA = mindiagDA,
                   bM = bM, bMH = bMH,  bBA = bBA,  bPR = bPR, bDA = bDA, int = integerMarker,
                   fRatioStationary = fRatioStationary, alsorev = alsorev)
  rownames(tab) = navn
  if(nr == 1)
    tab = tab[,c(1,2,3, 5:7)]
  else if (nr == 2){
    tab = tab[,c(1,2,3, 5:7, 10:11,14)]
    tab = tab[tab$int == "Y", ]
    tab = tab[,-9]
  }
  else if (nr == 3)
    tab = tab[,c(1:8, 10:14)]
  tab
}

checkBounded = function(dat = NULL){
  foo = dat[, c("bM", "bMH", "bBA", "bPR",  "bDA")]
  apply(apply(foo, 2, function(x) x=="Y"),2, sum)
}
