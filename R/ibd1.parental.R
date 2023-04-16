#' Estimates for a pair of non-inbred individuals the probabilities of parental
#' origin when IBD is 1
#'
#' Assume IBD is 1 for a pair of non-inbred individuals. The function estimates
#' the probabilities of the four possible combinations of parental origin: (i)
#' paternal-paternal, (ii) paternal-maternal, (iii) maternal-paternal and (iv)
#' maternal-maternal.
#'
#'
#' @param x A `ped` object.
#' @param id.pair Integer vector of length 2 giving the pair of individuals.
#' @param Nsim Integer. The number of simulations.
#' @param cM Double. The length of the `uniformMap`.
#' @param verbose Logical.
#' @return The averaged values.
#' @author Magnus Dehli Vigeland and Thore Egeland.
#' @examples
#'
#' library(ibdsim2)
#' x = doubleFirstCousins()
#' ids = leaves(x)
#' ibd1.parental(x, ids, 10)

#' @importFrom ibdsim2 findPattern ibdsim uniformMap
#' @export
#' 
ibd1.parental = function(x, id.pair, Nsim, cM = 10000, verbose = FALSE) {
    map = ibdsim2::uniformMap(cM = cM)
    sim = ibdsim2::ibdsim(x, N = Nsim, map = map, model = "haldane", verbose = verbose)
    alpha.sample = vapply(sim, function(s) {
        a = ibdsim2::findPattern(s, list(carriers = id.pair))
        if(nrow(a) == 0) return(c(NA,NA,NA,NA))
        ibd1 = a[,6:9]

        if(!is.matrix(ibd1))
           ibd1 = matrix(as.integer(ibd1), nrow =1 , ncol = 4,byrow = T)

        # sum lengths patpat , patmat, matpat, matmat
        patpat = sum((ibd1[,1] == ibd1[,3]) & (ibd1[,2] != ibd1[,4]))
        patmat = sum((ibd1[,1] == ibd1[,4]) & (ibd1[,2] != ibd1[,3]))
        matpat = sum((ibd1[,2] == ibd1[,3]) & (ibd1[,1] != ibd1[,4]))       
        matmat = sum((ibd1[,2] == ibd1[,4]) & (ibd1[,1] != ibd1[,3]))
        len = patpat + patmat + matpat + matmat
        c(patpat, patmat, matpat, matmat)/len
    }, numeric(4))
    
    alpha.hat = rowMeans(alpha.sample, na.rm = T)
    names(alpha.hat) = c("patpat", "patmat", "matpat", "matmat")
    alpha.hat 
}
