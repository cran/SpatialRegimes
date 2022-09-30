
sswf <- function(data, id, coly, colx, method=1,ind_col, lat, long, tau.ch) 
{
    dy = as.data.frame(data[,coly])
    colnames(dy) = c("dy")
    dx = data[,colx]
    dyx = as.data.frame(cbind(dy,dx))
    
    ### OLS
    
    if (method == 1)
       return(sum((lm(dy ~ ., data = dyx[id,])$residuals)^2))

    ### QUANTILE 

    if (method == 2)
      if (dim(dyx[id,])[1]<5) {
        return(sum(abs(lm(dy ~ ., data = dyx[id,])$residuals)))
      }
    else {
      rhoq <- function(u,tau=.5)u*(tau - (u < 0))
      resid = rq(dy ~ ., data = dyx[id,], tau=tau.ch)$residuals
      return(sum(rhoq(resid, tau.ch)))
    }  
    
    ### ELSE
    
    else 
      
        return(0)
}

