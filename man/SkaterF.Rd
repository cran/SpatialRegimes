\name{SkaterF}
\alias{SkaterF}

\title{Spatial constrained clusterwise regression by Spatial 'K'luster Analysis by Tree Edge Removal}

\description{This function implements a spatial constrained clusterwise regression based on the Skater procedure by Assuncao et al. (2002).}

\usage{SkaterF(edges,data,coly,colx,ncuts,crit,method,ind_col)}

\arguments{
  \item{edges}{A matrix with 2 colums with each row is an edge.}
  \item{data}{A data.frame with the informations over nodes.}
  \item{coly}{The dependent variable in the c("y_ols") form.}
  \item{colx}{The covariates in the c("x1","x2") form.}
  \item{ncuts}{The number of cuts.}
  \item{crit}{A scalar or two dimensional vector with with criteria for groups. Examples: limits of group size or limits of population size. If scalar, is the minimum criteria for groups.}
  \item{method}{1 (default) for OLS.}
  \item{ind_col}{Parameter still not used in this version.}
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
Assuncao, R.M., Lage J.P., and Reis, E.A. (2002). Analise de conglomerados espaciais via arvore geradora minima. Revista Brasileira de Estatistica, 62, 1-23.

Assuncao, R. M, Neves, M. C., Camara, G. and Freitas, C. da C. (2006). Efficient regionalization techniques for socio-economic geographical units using minimum spanning trees. International Journal of Geographical Information Science Vol. 20, No. 7, August 2006, 797-811
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

sk = SkaterF(edges = edges1, 
              data= SimData, 
              coly = coly1,
              colx= colx1, 
              ncuts=ncuts1, 
              crit=crit1,
              method=1)
			  
SimData$regimes = sk$groups
plot(lat~long,SimData,col=regimes,pch=16)              
}
