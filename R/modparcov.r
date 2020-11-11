
modparcov <- function (y, x, w)   #GWR/LWR
{ 
  xpw        <- t(x)%*%w                        #w matrice diagonale
  invxpwx    <- solve(xpw%*%x)         
  b          <- invxpwx%*%(xpw%*%y)             #kx1
  resid      <- y - x %*% b                     #nx1
  qresid     <- t(resid)%*%w%*%resid  #somma dei quadrati dei residui PESATI
  c          <- invxpwx%*%xpw
  s          <- x%*%c
  v1         <- sum(diag(s))
  v2         <- sum(diag(t(s)%*%s))
  varresid   <- qresid/(nrow(x)-2*v1+v2)        #sigma2
  ccp        <- c%*%t(c)
  res        <- list()
  res$par    <- b                               #kx1
  res$cov    <- ccp*as.numeric(varresid)        #kxk
  res
}