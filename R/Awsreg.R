

Awsreg <- function (data, coly, colx, kernel="bisquare",kernel2="gaussian",coords,bw,tau,niter,conv,eta,numout,sout)  
{                     
  x = as.matrix(data[,colx])
  y = as.matrix(data[,coly])  
  
  codic  <- 1:nrow(x)
  clas   <- rep(0,nrow(x))
  woutl  <- rep(0,nrow(x))
  
  vdist  <- as.matrix(gw.dist(coords,focus=0,p=2,theta=0,longlat=F)) #dmat  
  knn    <- knn2nb(knearneigh(coords,k=bw))
  nbdist <- nbdists(knn,coords)
  
  if (kernel=="exponential") {
    w<-matrix(0,nrow(vdist),ncol(vdist))
    for (i in 1:nrow(vdist)){
      for (j in 1:ncol(vdist)){
        w[i,j] <- exp(-vdist[i,j]/max(nbdist[[i]]))           #esponenziale
      }}
  }
  
  if (kernel=="gaussian") {
    w<-matrix(0,nrow(vdist),ncol(vdist))
    for (i in 1:nrow(vdist)){
      for (j in 1:ncol(vdist)){
        w[i,j] <- exp(-0.5*(vdist[i,j]/max(nbdist[[i]]))^2)   #gaussiana
      }}
  }
  
  if (kernel=="bisquare") {
    w<-matrix(0,nrow(vdist),ncol(vdist))                     #bisquare
    for (i in 1:nrow(vdist)){
      for (j in 1:ncol(vdist)){
        if (vdist[i,j]<max(nbdist[[i]])) w[i,j] <- (1-(vdist[i,j]/max(nbdist[[i]]))^2)^2   
        else                             w[i,j] <- 0.00001
      }}
  }
  
  if (kernel=="tricube") {
    w<-matrix(0,nrow(vdist),ncol(vdist))                     #tricube
    for (i in 1:nrow(vdist)){
      for (j in 1:ncol(vdist)){
        if (vdist[i,j]<max(nbdist[[i]])) w[i,j] <- (1-(vdist[i,j]/max(nbdist[[i]]))^3)^3   
        else                             w[i,j] <- 0.00001
      }}
  }
  
  wnew   <- matrix(0,nrow(x),nrow(x))
  parm   <- matrix(0,nrow(x),ncol(x))   #matrice parametri nxk
  qresid <- rep(0,length=nrow(x))       #vett.somme quad. resid. 
  listcov <- list()                     #lista di cov matrices
  
  nout   <- 0
  j      <- 0
  differ <- 1
  while (differ > conv && (j <- j+1)<=niter)
  {
    listw<-list()
    for (i in 1:ncol(w))
    {
      listw[[i]]<-diag(w[,i])
    }
    for (i in 1:nrow(x))
    {
      wls         <- modparcov(y, x, listw[[i]]) 
      parm[i,]    <- wls$par                     #parm[i,] i-esima riga di beta 1xk 
      listcov[[i]]<- wls$cov
    }
    for (u in 1:(nrow(x)-1))                     #test
    {
      wnew[u,u] <- 1
      for (v in (u+1):nrow(x))
      {
        if (w[u,v] < conv)
        {
          wnew[u,v] <- 0
          wnew[v,u] <- wnew[u,v]
        }
        else
        {
          if (kernel2=="exponential"){
            wnew[u,v] <- exp(-testpar(parm[u,],parm[v,],listcov[[u]],listcov[[v]],listw[[u]],listw[[v]])*tau)        #esponenziale
          }
          if (kernel2=="gaussian"){
            wnew[u,v] <- exp(-0.5*(testpar(parm[u,],parm[v,],listcov[[u]],listcov[[v]],listw[[u]],listw[[v]])*tau)^2) #gaussiana  
          }
          wnew[v,u] <- wnew[u,v] #per simmetria
        }
      }
    }
    wnew[nrow(x),nrow(x)] <- 1
    
    if (kernel2=="exponential") {
      wnew2<-matrix(0,nrow(vdist),ncol(vdist))                 #esponenziale
      for (i in 1:nrow(vdist)){
        for (h in 1:ncol(vdist)){
          wnew2[i,h]<-wnew[i,h]*exp(-vdist[i,h]/(max(nbdist[[i]])*j))  
        }}
      wnew<-wnew2
    }
    
    if (kernel2=="gaussian") {
      wnew2<-matrix(0,nrow(vdist),ncol(vdist))                  #gaussiana
      for (i in 1:nrow(vdist)){
        for (h in 1:ncol(vdist)){
          wnew2[i,h]<-wnew[i,h]*exp(-0.5*(vdist[i,h]/(max(nbdist[[i]])*j))^2)  
        }}
      wnew<-wnew2
    }
    
    differ <- max(abs(w-wnew))
    w      <- eta*w+(1-eta)*wnew
    
    for (i in 1:nrow(x))            #controllo outlier
    {
      wcontr    <- w[,i]
      woutl[i]  <- length(wcontr[wcontr > sout])
    }
    selez  <- which(woutl[1:nrow(x)] < numout)
    if (length(selez) > 0)
    {
      nout   <- nout + length(selez)
      x      <- x[-selez,]
      y      <- y[-selez,]
      codic  <- codic[-selez]
      vdist  <- vdist[-selez,-selez]
      w      <- w[-selez,-selez]
      wnew   <- wnew[-selez,-selez]             
    }
    
    message(paste0("Loop #: ", floor(j), " - Differ: ",differ, " - nout:", nout))

  }
  ncl <- 0
  for (i in 1:(nrow(x)-1))
  {
    if (clas[codic[i]] < 1)
    {
      ncl     <- ncl + 1
      clas[codic[i]] <- ncl
      for (j in (i+1):nrow(x))
      {
        if (w[i,j] > 0.9) clas[codic[j]] <- clas[codic[i]]
      }
    }
  }
  if (clas[nrow(x)] < 1) clas[nrow(x)] <- ncl + 1
  clas[clas == 0] <- NA
  
  res <- NA
  res$groups <- clas
  attr(res, "class") <- "Awsreg"
  return(res)
  
}





