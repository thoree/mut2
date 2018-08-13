lik2 <- function(g1, g2, n, p, M, kappa, alpha, theta, beta=0.5){
  l0 <- function(a, b, c, d, p, theta) {
    pa <- p[a]
    pb <- (b==a)*theta+(1-theta)*p[b]
    pc <- (((c==a)+(c==b))*theta+(1-theta)*p[c])/(1+theta)
    pd <- (((d==a)+(d==b)+(d==c))*theta+(1-theta)*p[d])/(1+2*theta)
    2^(-(a == b) - (c == d))*4*pa*pb*pc*pd
    
  }
  l1 <- function(a, b, c, d, g, p, M, theta) {
    Mg <- M %^% g
    pa <- p[a]
    pb <- (b==a)*theta+(1-theta)*p[b]
    pc <- (((c==a)+(c==b))*theta+(1-theta)*p[c])/(1+theta)
    pd <- (((d==a)+(d==b))*theta+(1-theta)*p[d])/(1+theta)
    2^(-(a == b) - (c == d))*pa*pb*(pd*(Mg[a, c]+Mg[b, c])+pc*(Mg[a, d] + Mg[b,d]))
  }
  l2 <- function(a, b, c, d, nP, nM, p, M, theta) {
    MnP <- M %^% nP
    MnD <- M %^% nM
    pa <- p[a]
    pb <- (b==a)*theta+(1-theta)*p[b]
    2^(-(a == b) - (c == d))*pa*pb*
      (MnP[a, c] * MnD[b, d] + MnP[b, c] * MnD[a, d] + 
         MnP[a, d] * MnD[b, c] + MnP[b, d] * MnD[a, c])
  }
    a <- g1[1]
    b <- g1[2]
    c <- g2[1]
    d <- g2[2]
    beta <- c(beta, 1-beta)
    lik <- kappa[1]* l0(a, b, c, d, p, theta)+ 
		kappa[2] * alpha[1] * l1(a, b, c, d, g = n[1], p, M, theta) +
		kappa[2] * alpha[2] * l1(a, b, c, d, g = n[2], p, M, theta) + 
		kappa[2] * alpha[3] * l1(a, b, c, d, g = n[3], p, M, theta) + 
		kappa[2] * alpha[4] * l1(a, b, c, d, g = n[4], p, M, theta) + 
		kappa[3] * beta[1] * l2(a, b, c, d, n[1], n[4], p, M, theta)+ 
		kappa[3] * beta[2] * l2(a, b, c, d, n[2], n[3], p, M, theta)
	lik
}

