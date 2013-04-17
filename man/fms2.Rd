\name{fms2}
\alias{fms2}
\title{Covariance Estimation by modified PCA}
\usage{
  fms2(x, weight = seq(1, 3, length = nobs), center = TRUE,
    frac.var = 0.5, iter.max = 1, nfac.miss = 1,
    full.min = 20, reg.min = 40, sd.min = 20,
    quan.sd = 0.9, tol = 0.001, zero.load = FALSE,
    range.factors = c(20, 20), lambda = 0,
    minunique = 0.02, shrinkb = 0.3, shrinkv = shrinkb,
    shrinkr = 0.9)
}
\arguments{
  \item{x}{matrix or dataframe of timeseries returns}

  \item{weight}{weights in estimation}

  \item{center}{flag to center}

  \item{frac.var}{controls auto-selection of number of
  factord}

  \item{iter.max}{maximum number of iterations}

  \item{nfac.miss}{number of factors to estimate if data is
  missing}

  \item{full.min}{minimum acceptable number of NA-free
  columns}

  \item{reg.min}{minimum dates to do regression}

  \item{sd.min}{minimum dates to estimate vol}

  \item{quan.sd}{missing vol assigned this quantile}

  \item{tol}{estimation tolerance}

  \item{zero.load}{flag to use zero loadings for columns
  with missing}

  \item{range.factors}{range of factors to estimate, as a
  function of valid data length}

  \item{lambda}{exponent on eigenvalue for shrinkage}

  \item{minunique}{minimum uniqueness}

  \item{shrinkb}{shrinkage for factor 1}

  \item{shrinkv}{shrinkage for vol}

  \item{shrinkr}{shrinkage for regressed loadings}
}
\value{
  list(loadings fmp hpl method full uniqueness sdev qua
  weight call)
}
\description{
  Covariance Estimation by modified PCA
}
\section{Details}{
  more detail on the underlying algorithm may be found in
  documentation for BurStFin
}
\author{
  Giles Heywood from Pat Burns original
}

