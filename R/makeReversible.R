ER <- function(Q, pi) {
  # Calculates expected mutation rate of a mutation matrix X
  # Allele frequencies pi
  Dp = diag(pi)
  1-sum(diag(Dp%*%Q)) 
}

makeReversible <- function(Q, pi, alpha=NULL){ 
# Q proposed mutation matrix
# pi allele frequencies
# If alpha = NULL, alpha will be calculated
# to make alternative 3 matrix, see below, have same mutation rate as Q.
# Otherwise specified alpha is used.

# Three time reversible matrices, P1, P2, P3, are returned: 
# (1) Metropolis-Hastings (MP)
# (2) Averaged
# (3) alpha modification of 1, i.e., MP
# Minimimum diagonals and expected mutation rates are also returned
  

n <- length(pi)   
P1 <- P2 <- P3 <- matrix(ncol = n, nrow = n, 0)
#P1 
for (i in 1:n)
  for (j in 1:n){
    if (pi[i]*Q[i,j]>0)
        P1[i,j]  <- Q[i,j]*min(1,(pi[j]*Q[j,i])/(pi[i]*Q[i,j]))
  }
diag(P1) <- 0
s <- apply(P1,1,sum)
for (i in 1:n)
  P1[i,i] <- 1-s[i]
#P2
  for (i in 1:n)
    for (j in 1:n)
      P2[i,j]  <- (pi[i]*Q[i,j]+pi[j]*Q[j,i])/(2*pi[i])

diag(P2) <- 0
s <- apply(P2,1,sum)
for (i in 1:n)
  P2[i,i] <- 1-s[i]

#P3
R.Q = ER(Q, pi)
R.P1 = ER(P1, pi)
R.P2 = ER(P2, pi)
if(is.null(alpha))  alpha = R.Q/R.P1
P3 = alpha*P1+(1-alpha)*diag(rep(1,length(pi)))
R.P3 = ER(P3, pi)
list(P1 = P1, P2 = P2, P3 = P3,
     'Minimum on diagonal' = 
      c('Q' = min(diag(Q)), 'P1' = min(diag(P1)), 'P2' = min(diag(P2)), 'P3' = min(diag(P3)) ),
     'Expected mutation rates' = c('R.Q' = R.Q, 'R.P1' = R.P1, 'R.P2' = R.P2, 'R.P3' = R.P3))
}

