
sswf <- function(data, id, coly, colx, method=1,ind_col) 
{
    dy = as.data.frame(data[,coly])
    colnames(dy) = c("dy")
    dx = data[,colx]
    dyx = as.data.frame(cbind(dy,dx))
    
    ### OLS
    
    if (method == 1)
       return(sum(abs(lm(dy ~ ., data = dyx[id,])$residuals)))

  
    ### ELSE
    
    else 
      
        return(0)
}

