lik2X <- function(g1, g2, n, p, M, kappa0){
  if(length(g1) != 1)
    stop("First argument must be an integer in 1,..., length(p)")
  if(length(g2) != 2)
    stop("Second argument must be two integers in 1,..., length(p)")
  if (!is.numeric(p) || any(p <= 0)) 
    stop("frequencies must be a vector of positive numbers.")
  if (round(sum(p), 6) != 1) 
    stop("The frequencies must sum to 1.")
  if(!is.matrix(M))
    stop("M must be a mutation matrix")
  if (any(as.vector(M) < 0)) 
    stop("The  mutation matrix cannot have negative entries in locus")
  if (any(round(apply(M, 1, sum), 6) != 1)) 
    stop("The rows in the  mutation matrix must sum to 1 in locus")
  if( kappa0 < 0 || kappa0 > 1)
    stop("kappa0 must be in[0,1]")

  a = g1
  b = g2[1]
  cc = g2[2]
  l0 = 2^(-(b == cc))*2*p[a]*p[b]*p[cc]
  Mg = M %^% n
  l1 = 2^(-(b == cc))*p[a]*(p[cc]*Mg[a, b] + p[b]*Mg[a, cc])
  lik = kappa0*l0 +  (1-kappa0)*l1
  lik
}
