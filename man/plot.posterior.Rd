% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcleod_main_function.R
\name{plot.posterior}
\alias{plot.posterior}
\title{Plot the posterior mixing distribution, estimated for the data from a mcleod model}
\usage{
\method{plot}{posterior}(mcleod.obj, plot_only_point_estimate = F)
}
\arguments{
\item{mcleod.obj}{Object received from \code{\link{mcleod}}.}

\item{plot_only_point_estimate}{Logical value indicating if only the posterior mean should be plotted (using a red line). If set to F (default value), the posterior density will also be plotted.}
}
\value{
NONE
}
\description{
The function receives an object returned from \code{\link{mcleod}} and plots the posterior estimate for the mixing distribution. The posterior mean is shown in red (means for each point in the support of the mixing distribution). A point cloud shows the posterior distribution at each point of the support.
}
\examples{
 # For full description of package model and workflow,
 # including this function, Type browseVignettes(package = 'mcleod') 
 # in the R console and check the package vignette

library(mcleod)
##################################
# Example 1: Binomial sampling distribution:
# The function plots the posterior mixing distribution
##################################
N = 30
K = 300
set.seed(1)
u = sample(c(0,1),size = K,replace = T)
x = rbinom(K,size = N,prob =inv.log.odds(rnorm(K,-1+3*u,sd = 0.3)))
n = rep(N,K)

#fit model
res = mcleod(x, n)
# Plot CDF with distribution of posterior samples
plot.posterior(res)

#############################
# Example 2: Binomial regression with a random normal intercept
# The function plots the posterior distribution of the random intercept
#############################
# Generate data:
N = 30 #Number of draws per binomial observations
K = 200 #Number of samples
set.seed(1)
covariates = matrix(rnorm(K*2,sd = 0.5),nrow = K) #Generate covariates
colnames(covariates) = c('covariate 1','covariate 2')
#define slopes:
real_beta_1 = -1
real_beta_2 = 1
#sample
x = rbinom(K,size = N,
prob = inv.log.odds(rcauchy(K,location = 0,scale = 0.5) +
real_beta_1*covariates[,1] + real_beta_2*covariates[,2]))
n = rep(N,K)

# Fit model:
res = mcleod(x, n, covariates = covariates)

#Posterior distribution of random intercept:
plot.posterior(res)
}
