% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mu_sigma.R
\name{get_mu_sigma}
\alias{get_mu_sigma}
\title{Function to transform mean and standard deviation from "original" baseline values to parameters of the normal distribution of the log values.}
\usage{
get_mu_sigma(mu_raw_0, sd_raw_0, reductrate_6, reductrate_12, rho)
}
\arguments{
\item{mu_raw_0}{num mean value of original data at baseline}

\item{sd_raw_0}{num standard deviation of original data at baseline}

\item{reductrate_6}{vector of reduction rates after 6 months for each dose (between 0 and 1)}

\item{reductrate_12}{vector of reduction rates after 12 months for each dose (between 0 and 1)}

\item{rho}{correlation between baseline and 6 months follow-up observations and between 6 months and 12 months follow-up observations}
}
\value{
list of mean values and variance covariance matrix
}
\description{
Function to transform mean and standard deviation from "original" baseline values to parameters of the normal distribution of the log values.
}
\details{
eWHORM simulations
}
\author{
Marta Bofill Roig
}
\keyword{internal}
