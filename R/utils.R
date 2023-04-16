#Likelihood for the unrelated case for individuals
# with genotypes a/b and c/d, allele freq p and theta 
l0 <- function(a, b, c, d, p, theta) {
  pa <- p[a]
  pb <- (b==a)*theta+(1-theta)*p[b]
  pc <- (((c==a)+(c==b))*theta+(1-theta)*p[c])/(1+theta)
  pd <- (((d==a)+(d==b)+(d==c))*theta+(1-theta)*p[d])/(1+2*theta)
  2^(-(a == b) - (c == d))*4*pa*pb*pc*pd
}