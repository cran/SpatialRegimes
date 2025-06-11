\name{SkaterF}
\alias{SkaterF}

\title{Spatial constrained clusterwise regression by Spatial 'K'luster Analysis by Tree Edge Removal}

\description{This function implements a spatial constrained clusterwise regression based on the Skater procedure by Assuncao et al. (2002).}

\usage{SkaterF(edges,data,coly,colx,ncuts,crit,method=1,ind_col,lat,long,tau.ch)}

\arguments{
  \item{edges}{A matrix with 2 colums with each row is an edge.}
  \item{data}{A data.frame with the informations over nodes.}
  \item{coly}{The dependent variable in the c("y_ols") form.}
  \item{colx}{The covariates in the c("x1","x2") form.}
  \item{ncuts}{The number of cuts.}
  \item{crit}{A scalar or two dimensional vector with with criteria for groups. Examples: limits of group size or limits of population size. If scalar, is the minimum criteria for groups.}
  \item{method}{1 (default) for OLS, 2 for Quantile regression, 3 for logit}
  \item{ind_col}{Parameter still not used in this version.}
  \item{lat}{Parameter still not used in this version.}
  \item{long}{Parameter still not used in this version.}
  \item{tau.ch}{Chosen quantile (for method = 2).}
}


\value{
A object of skaterF class with:
  \item{groups}{A vector with length equal the number of nodes. Each position identifies the group of node.}
  \item{edges.groups}{A list of length equal the number of groups with each element is a set of edges}
  \item{not.prune}{A vector identifying the groups with are not candidates to partition.}
  \item{candidates}{A vector identifying the groups with are candidates to partition.}
  \item{ssto}{The total dissimilarity in each step of edge removal.}
}

\details{Author really thanks Renato M. Assuncao and Elias T. Krainski for their original code (skater, library spdep).}

\author{F. Vidoli}

\references{

For method = 1: F. Vidoli, G. Pignataro, and R. Benedetti. (2022) "Identification of spatial regimes of the production function of italian hospitals through spatially constrained cluster-wise regression. In: Socio-Economic Planning Sciences, page 101223, doi: https://doi.org/10.1016/j.seps.2022.101223

For method = 2: Vidoli, F., Sacchi A. & Sanchez Carrera E. (2025) "Spatial regimes in heterogeneous territories: The efficiency of local public spending" In: Economic modelling https://doi.org/10.1016/j.econmod.2025.107139

}

%\seealso{ \code{\link{ci_bod}}, \code{\link{ci_mpi}}}

\examples{

data(SimData)
coords = cbind(SimData$long, SimData$lat)

#######################

neighbours = tri2nb(coords, row.names = NULL)
bh.nb <-  neighbours
lcosts <- nbcosts(bh.nb, SimData)
nb <- nb2listw(bh.nb, lcosts, style="B")
mst.bh <- mstree(nb,5)
edges1 = mst.bh[,1:2]

#######################

ncuts1 = 4
crit1 = 10
coly1 = c("y_ols")
colx1 = c("A","L","K")

# OLS
sk = SkaterF(edges = edges1, 
              data= SimData, 
              coly = coly1,
              colx= colx1, 
              ncuts=ncuts1, 
              crit=crit1,
              method=1)
			  
SimData$regimes = sk$groups
# plot(lat~long,SimData,col=regimes,pch=16)   

## quantile 0.8
# sk2 = SkaterF(edges = edges1, 
#              data= SimData, 
#              coly = coly1,
#              colx= colx1, 
#              ncuts=ncuts1, 
#              crit=crit1,
#              method=2,tau.ch=0.8)
#			  
# SimData$regimes_q = sk2$groups
# plot(lat~long,SimData,col=regimes_q,pch=16)   


}
