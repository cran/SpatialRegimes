
resregby <- function (y, x, cod) # solo lineare semplice e pericolosissimo
{
  codic<-as.numeric(levels(as.factor(cod)))
  res<-0
  for (i in codic)
  {
    b<-cov(y[cod==i],x[cod==i])/var(x[cod==i])
    a<-mean(y[cod==i])-b*mean(x[cod==i])
    res<-res+sum((y[cod==i]-a-b*x[cod==i])^2)
  }
  res
}
