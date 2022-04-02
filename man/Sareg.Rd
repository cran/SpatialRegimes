\name{Sareg}
\alias{Sareg}

\title{Spatial clusterwise regression by a constrained version of the Simulated Annealing}

\description{This function implements a spatial clusterwise regression based on the procedure suggested by Postiglione \emph{et al.} (2013).}

\usage{Sareg(data, coly,colx, cont, intemp, rho, niter, subit, ncl, bcont)}

\arguments{
  \item{data}{A data.frame}
  \item{coly}{The dependent variable in the c("y_ols") form.}
  \item{colx}{The covariates in the c("x1","x2") form.}
  \item{cont}{The contiguity matrix.}
  \item{intemp}{The initial temperature.}
  \item{rho}{The temperature decay rate parameter.}
  \item{niter}{The maximum number of iterations.}
  \item{subit}{The number of sub-iterations for each iteration.}
  \item{ncl}{The number of clusters.}
  \item{bcont}{A parameter that regulates the penalty of simulated annealing in non-contiguous configurations of the clusters.}
  }


\value{
A object of Sareg class with:
  \item{groups}{Estimated clusters.}
}

%\details{....}

\author{R. Benedetti}

\references{
Postiglione, P., Benedetti, R., and Andreano, M.S. (2013). "\emph{Using Constrained Optimization for the Identification of Convergence Clubs}", Computational Economics, 42, 151-174.
}

%\seealso{ \code{\link{ci_bod}}, \code{\link{ci_mpi}}}

\examples{

data(SimData)
SimData = SimData[1:50,]
coords = cbind(SimData$long, SimData$lat)

#######################

dmat   <-gw.dist(coords,focus=0,p=2,theta=0,longlat=FALSE)
W      <- matrix(0,nrow(dmat),ncol(dmat))
W[dmat < 0.2] <- 1
diag(W)       <- 0

#######################

sa <- Sareg(data=SimData,
               coly = c("y_ols"),
               colx = c("A"),
               W,
               intemp=0.5,
               rho=0.96,
               niter=30,
               subit=3,
               ncl=2,
               bcont=-4)

SimData$regimes = sa$groups
plot(lat~long,SimData,col=regimes,pch=16)
}
