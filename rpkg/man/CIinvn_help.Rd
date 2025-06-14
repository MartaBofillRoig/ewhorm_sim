% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIinvn_help.R
\name{CIinvn_help}
\alias{CIinvn_help}
\title{help function to calculate quantiles (e.g., 1-alpha lower confidence interval or median) of concordance with inverse normal method to combine concordance of stage 1 and stage 2.}
\usage{
CIinvn_help(theta, thetahat1, thetahat2, value, N1, N2)
}
\arguments{
\item{theta}{unknown confidence limit - to be found by uniroot function}

\item{thetahat1}{estimated empirical concordance of stage 1}

\item{thetahat2}{estimated empirical concordance of stage 2}

\item{value}{quantile - significance level or median}

\item{N1}{sample size of group 1}

\item{N2}{sample size of group 2}
}
\value{
variance of concordance
}
\description{
help function to calculate 1-alpha confidence interval of concordance with pooled variances (modified Hanley-McNeil approach) according to Newcombe, 2006.
}
\details{
eWHORM simulations
}
\author{
Sonja Zehetmayer
}
\keyword{internal}
