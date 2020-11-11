


ncont <- function (cont, cod)
{
  codic<-as.numeric(levels(as.factor(cod)))
  nc<-0
  for (i in codic)
  {
    nc<-nc+sum(cont[cod==i,cod==i])
  }
  nc
}