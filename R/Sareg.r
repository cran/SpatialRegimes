


Sareg <- function (data, coly, colx, cont, intemp=1, rho, niter, subit, ncl, bcont)
{
  
  x = as.matrix(data[,colx])
  y = as.matrix(data[,coly]) 
  
  codic<-1:ncl
  size<-length(y)
  cl<-sample(1:ncl,size,replace=T)
  res<-resregby(y,x,cl)
  nco<-ncont(cont,cl)
  v<-0
  nch<-1
  while (nch > 0 && (v <- v+1)<niter)
  {
    nch<-0
#    print(c(floor(v),res,nco,intemp))
    
    message(paste0("Loop #: ", floor(v), " - Res: ",res, " - Nco: ", nco, " - Intemp: ", intemp))
    
    
    
    for (k in 1:subit)
    {
      for (i in 1:size)
      {
        oldcl<-cl[i]
        cl[i]<-sample(codic[-cl[i]],1)
        newres<-resregby(y,x,cl)
        newnco<-ncont(cont,cl)
        genr<-log(runif(1))
        objfun<-(((res-newres)+bcont*(nco-newnco))/intemp)
        if (genr < objfun)
        {
          res<-newres
          nco<-newnco
          nch<-nch+1
        }
        else
        {
          cl[i]<-oldcl
        }
      }
    }
    intemp<-intemp*rho
  }
  ress <- NA
  ress$groups = cl
  attr(ress, "class") <- "Sareg"
  return(ress)
  
}

