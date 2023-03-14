# added
LRtrio = function(gMan, gMother, gChild,  afreq, M, check = TRUE ){
  if(is.null(rownames(M)) & is.null(colnames(M)) & is.null(names(afreq))){
    navn = 1:length(afreq)
	dimnames(M) = list(navn, navn)
  }
  if(is.null(rownames(M)) | is.null(colnames(M)))
    stop("Row - and column names of M must be supplied")
  if(!(all(rownames(M) == colnames(M))))
    stop("Row - and column names of M must equal and equally sorted")
  als = rownames(M)
  afreq2 = as.double(afreq)
  names(afreq) = als
  a = gMother[1]
  b = gMother[2]
  cc = gChild[1] 
  d = gChild[2]
  e = gMan[1]
  f = gMan[2]
  mac = M[a, cc]
  mbc = M[b, cc]  
  med = M[e, d]
  mfd = M[f, d]
  mad = M[a, d]
  mbd = M[b, d]
  mec = M[e, cc]
  mfc = M[f, cc]
  pa = afreq[a]
  pb = afreq[b] 
  pc = afreq[cc]
  pd = afreq[d]
  pe = afreq[e]
  pf = afreq[f]
  pp = as.double(c(pa, pb, pc, pd, pe, pf))
  I1 = 1 + (a!=b)
  I2 = 1 + (cc!=d)
  I3 = 1 + (e!=f)
  II = c(I1,I2,I3)
  LH1 = prod(pp)*prod(II)
  LH2 = I3*pe*pf*I1*I2*(1/4)*pa*pb*((mac+mbc)*pd+(mad+mbd)*pc)
  term1 = pa*pb*pe*pf*prod(II)*(1/8)
  term2 =((mac+mbc)*(med+mfd)+(mad+mbd)*(mec+mfc))
  LH3 = term1*term2
  liks = as.double(c(LH1,LH2,LH3))
  if (!check){
    result = rbind(liks)
	rownames(result) = c("lik.formula")
  }
  else {
    M = mutationMatrix("custom", matrix = M, alleles = als)
    Man = singleton("MA")
    mMan = marker(Man, MA = gMan, alleles = als, afreq = afreq2, 
                  mutmod = M)
    Mother = singleton("MO")
    mMother = marker(Mother, MO = gMother, alleles = als, afreq = afreq2, 
                  mutmod = M)
    Child = singleton("CH")
    mChild = marker(Child, CH= gChild, alleles = als, afreq = afreq2, 
                     mutmod = M)
    LH1.suite = likelihood(list(Man, Mother, Child), list(mMan, mMother, mChild))
    
    H2 = nuclearPed(father = "NN", mother = "Mother", child = "Child")
    mH2 = marker(H2,  Mother = gMother, Child = gChild, 
                 alleles = als, afreq = afreq2, mutmod = M)
    
    LH2.suite = likelihood(list(H2, Man), list(mH2, mMan))

    H3 = nuclearPed(father ="Man", mother = "Mother", child = "Child")
    mH3 = marker(H3, Man = gMan, Mother = gMother, Child = gChild, 
                 alleles = als, afreq = afreq, mutmod = M)
    LH3.suite = likelihood(H3, mH3)
    liks.suite = as.double(c(LH1.suite, LH2.suite, LH3.suite))
    result = rbind(liks, liks.suite)
    rownames(result) = c("lik.formula", "lik.pedprobr")
  }
  colnames(result) = c("unrelated", "Mother-Child", "Tio")
  result
}
