% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_dataind.R
\name{sim_dataind}
\alias{sim_dataind}
\title{Simulate individual data from a multi-arm trial with shared control}
\usage{
sim_dataind(n_arms, N, mu_0m, mu_6m, mu_12m, sg, rr, bound)
}
\arguments{
\item{n_arms}{number of arms (including control)}

\item{N}{total sample size}

\item{mu_6m}{6-month mean response per arm (vector of length \code{n_arm})}

\item{mu_12m}{12-month mean response per arm (vector of length \code{n_arm})}

\item{sg}{covariance matrix between 6- and 12-month responses assumed equal across arms (matrix of dim 2x2)}

\item{rr}{responder rate for each dose, which gives the proportion of patients with value 0 at follow-up}

\item{bound}{lower bound to define total responder in simulation study}
}
\value{
simulated data consisting of the responses at 6 and 12 months, treatment arm, for each subject.
}
\description{
Function to simulate trial data (1-stage, multiple arms)
}
\details{
eWHORM simulations
}
\author{
Marta Bofill Roig, Sonja Zehetmayer
}
\keyword{internal}
