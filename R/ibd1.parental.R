ibd1.parental = function(x, id.pair, Nsim, cM=10000, verbose=F, ...) {
    map = uniformMap(cM=cM)
    sim = ibdsim(x, sims=Nsim, map=map, model="haldane", verbose=verbose, ...)
    alpha.sample = vapply(sim, function(s) {
        a = alleleSummary(s, ids=id.pair, ibd.status=T)
        
        # Extract subset where IBD == 1
        ibd1 = a[a[, 'ibd'] == 1, , drop=FALSE]
        
        # If empty subset: return NA's. (These are later removed with na.rm=T.) 
        if(nrow(ibd1) == 0) return(c(NA,NA,NA,NA))
        
        # Columns 10:13 are ibd_pp,..., with exactly one 1 in each row (because ibd=1)
        # Hence the following are estimates of alpha_pp, ...
        len = ibd1[, 'length']
        colSums(len * ibd1[, 10:13, drop=FALSE])/sum(len)
    }, numeric(4))
    
    list(alpha.sample=alpha.sample, alpha.hat=rowMeans(alpha.sample, na.rm=T))   
}

uniformMap = function (Mb = NULL, cM = NULL, M = NULL, cm.per.mb = 1, chromosome = 1) 
{
  if (!is.null(cM) && !is.null(M)) 
    stop("Either 'cM' or 'M' must be NULL.")
  stopifnot(!is.null(cM) || !is.null(M) || !is.null(Mb))
  if (is.null(cM)) 
    if (!is.null(M)) 
      cM = M * 100
  else cM = cm.per.mb * Mb
  if (is.null(Mb)) 
    Mb = cM/cm.per.mb
  if (is.character(chromosome) && tolower(chromosome) == "x") 
    chromosome = 23
  map = switch(max(length(Mb), length(cM)), {
    m = cbind(Mb = c(0, Mb), cM = c(0, cM))
    list(male = m, female = m)
  }, {
    if (length(cM) == 1) cM = c(cM, cM)
    if (length(Mb) == 1) Mb = c(Mb, Mb)
    list(male = cbind(Mb = c(0, Mb[1]), cM = c(0, cM[1])), 
         female = cbind(Mb = c(0, Mb[2]), cM = c(0, cM[2])))
  })
  female_phys = as.numeric(map$female[2, 1])
  if (chromosome == 23) 
    map$male = NA
  else if (female_phys != map$male[2, 1]) 
    stop("Male and female chromosomes must have equal physical length.")
  structure(map, length_Mb = female_phys, chromosome = chromosome, 
            class = "chromosomeMap")
}

