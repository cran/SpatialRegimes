\name{SimData}
\alias{SimData}
\title{Simulated data for estimating spatial regimes.}
\description{Simulated production function like data for estimating spatial regimes; data has been generated for the paper "F. Vidoli, G. Pignataro, R. Benedetti, F. Pammolli, "\emph{Spatially constrained cluster-wise regression: optimal territorial areas in Italian health care}", forthcoming.}
\usage{data(SimData)}

\format{SimData is a simulated dataset with 500 observations and 7 variables.

\describe{
\item{long}{Longitude}	
\item{lat}{Latitude}	
\item{A}{Land input}	
\item{L}{Labour input}	
\item{K}{Capital input}	
\item{clu}{Real regime}	
\item{y_ols}{Production output}	
}

500 units (100 units for each of the 5 regimes) are generated and, for each unit, the longitude and latitude coordinates are randomly drawn by using two Uniform distributions from 0 to 50 and from -70 to 20, i.e. U(0,50) and U(-70,20), respectively.
Consequently, we set the matrix of covariates which include the constant, A, L and K variables by drawing from U(1.5,4).
For each regime, finally, a different (in the coefficients) spatial function is set assuming a linear functional form. 
More in particular, we set 5 different vectors of parameters (including the intercept): beta1 = (13,0.5,0.3,0.2), beta2 = (11,0.8,0.1,0.1), beta3 = (9,0.3,0.2,0.5), beta4 = (7,0.4,0.3,0.3) and beta5 = (5,0.2,0.6,0.2) and a normally distributed error term in N(0,1).
}

\references{F. Vidoli, G. Pignataro, R. Benedetti, F. Pammolli, "\emph{Spatially constrained cluster-wise regression: optimal territorial areas in Italian health care}", forthcoming
}

\author{Vidoli F.}

\examples{
data(SimData)
}