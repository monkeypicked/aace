\name{aace-package}
\alias{aace-package}
\alias{aace}
\docType{package}
\title{
Covariance Estimation
}
\description{
Covariance estimation using modified PCA, based on BurStFin

Recap of the factor model:  

R         n x 1       returns; 
a         n x 1       intercept;  
B         n x m       loading;  
G         n x m       'factor-mimicking portfolios';  
f         m x 1       factor score;  
S         n x n       covariance of R;  

R = a + B.f + e     (1)  

E[f.f'] = I(m)  orthonormal factors;  
E[e.e'] = psi   specific risk (diagonal)  

S = B.B' + psi      (2)  

f = R.G             (3)  

The function fms2 estimates B,G; returns a list which is informally of class ce.  

The remaining functions are accessors.  
}
\details{
\tabular{ll}{
Package: \tab aace\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-04-16\cr
License: \GPL (>=2)\cr
}

}
\author{
Giles Heywood

Maintainer: gilesheywood at gmail.com
}
\references{
}
\keyword{ package }
\seealso{
}
\examples{
}
