

testpar <- function (p1, p2, cov1, cov2, w1, w2)
{                            
  n1      <- as.numeric(sum(diag(w1))) 
  n2      <- as.numeric(sum(diag(w2)))
  pcov    <- (n1*cov1+n2*cov2)/(n1+n2)          #pooled var-cov beta
  t       <- t(p1-p2)%*%solve(pcov)%*%(p1-p2)
  t  
}