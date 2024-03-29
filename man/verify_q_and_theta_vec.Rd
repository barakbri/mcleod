% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CI_mixture_distribution.R
\name{verify_q_and_theta_vec}
\alias{verify_q_and_theta_vec}
\title{An internal function used to verify logic of q_vec and theta_vec}
\usage{
verify_q_and_theta_vec(
  q_vec,
  q_vec_for_computation,
  theta_vec,
  theta_vec_for_computation,
  sampling_distribution
)
}
\arguments{
\item{q_vec}{same parameter as in mcleod.CI.estimation.parameters}

\item{q_vec_for_computation}{same parameter as in mcleod.CI.estimation.parameters}

\item{theta_vec}{same parameter as in mcleod.CI.estimation.parameters}

\item{theta_vec_for_computation}{same parameter as in mcleod.CI.estimation.parameters}

\item{sampling_distribution}{same parameter as in mcleod.CI.estimation.parameters}
}
\value{
returns a list with q_vec,q_vec_for_computation,theta_vec,theta_vec_for_computation, n_q and n_theta, to be updated in CI_params object
}
\description{
An internal function used to verify that: 1) q_vec contains the same values as 1-q_vec (for binomial sampling distribution);
2) that all values in q_vec_for_computation and 1-q_vec_for_computation are found in q_vec;
3) q_vec values are between 0 and 1; 4) theta_vec contains all values found in theta_vec_for_computation; 5) for binomial sampling - theta_vec contains all values also found in -1*theta_vec.
}
\examples{
see internal code usages
}
\keyword{internal}
