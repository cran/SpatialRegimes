\name{Awsreg}
\alias{Awsreg}

\title{Spatial clusterwise regression by an iterated spatially weighted regression algorithm}

\description{This function implements a spatial clusterwise regression based on the procedure suggested by Andreano \emph{et al.} (2017) and Bille' \emph{et al.} (2017).}

\usage{Awsreg(data, coly,colx,kernel,kernel2,coords,bw,tau,niter,conv,eta,numout,sout)}

\arguments{
  \item{data}{A data.frame.}
  \item{coly}{The dependent variable in the c("y_ols") form.}
  \item{colx}{The covariates in the c("x1","x2") form.}
  \item{kernel}{Kernel function used to calculate distances between units (default is "bisquare", other values: "exponential", "gaussian","tricube").}
  \item{kernel2}{Kernel function used to calculate distances between units in the second step (default is "gaussian", other values: "exponential").}
  \item{coords}{The coordinates in terms of longitude and latitude.}
  \item{bw}{The bandwidth parameter of the initial weights.}
  \item{tau}{The confidence test parameter of the difference between regression parameters.}
  \item{niter}{The maximum number of iterations.}
  \item{conv}{The smallest accepted difference between the weights in two successive iterations.}
  \item{eta}{The parameter that regulates which is the weight of the weights of the previous iteration in the moving average that updates the new weights.}
  \item{numout}{The minimum number of areal units accepted for each cluster.}
  \item{sout}{Minimum value of weights such as to be considered equal to zero. Parameter used essentially to control clusters consisting of too few areal units.}
}

\value{
A object of Awsreg class with:
  \item{groups}{Estimated clusters.}
}

\details{Author really thanks Bille' A.G. for her contribution to revising the original code.}

\author{R. Benedetti}

\references{
Andreano, M.S., Benedetti, R., and Postiglione, P. (2017). "\emph{Spatial regimes in regional European growth: an iterated spatially weighted regression approach}", Quality & Quantity. 51, 6, 2665-2684.

Bille', A.G., Benedetti, R., and Postiglione, P. (2017). "\emph{A two-step approach to account for unobserved spatial heterogeneity}", Spatial Economic Analysis, 12, 4, 452-471. 
}



\examples{

data(SimData)
SimData = SimData[1:50,]
coords = cbind(SimData$long, SimData$lat)

#######################

dmat<-gw.dist(coords,focus=0,p=2,theta=0,longlat=FALSE)
bw<-bw.gwr(y_ols~A+L+K,
           data=SpatialPointsDataFrame(coords,SimData),
           approach="AIC",kernel="bisquare",
           adaptive=TRUE,p=2,theta=0,longlat=FALSE,dMat=dmat)

#######################

aws<-Awsreg(data=SimData,
           coly=c("y_ols"),
           colx=c("A","L","K"),
           kernel="bisquare",
           kernel2="gaussian",
           coords=coords,
           bw=bw,
           tau=0.001,
           niter=200,
           conv=0.001,
           eta=0.5,
           numout=15,
           sout=1e-05)

SimData$regimes = aws$groups
plot(lat~long,SimData,col=regimes,pch=16)
}
