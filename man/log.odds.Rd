% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Auxiliary_Functions.R
\name{log.odds}
\alias{log.odds}
\title{Compute log-odds from probability}
\usage{
\method{log}{odds}(p)
}
\arguments{
\item{p}{}
}
\value{
log(p/(1-p))
}
\description{
Compute log-odds from probability
}
\examples{

library(mcleod)
log.odds(0.5)

}
