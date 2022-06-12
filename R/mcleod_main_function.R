
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#load libraries required:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Constants and definitions:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CLASS.NAME.MCLEOD = 'mcleod.obj' # Class for the objectr returned from the main function. This is the result of deconvolution
CLASS.NAME.PRIOR.DEFINITION = 'mcleod.prior.def.obj' #object defining the prior: Hbeta (polya tree) or dirichlet tree, and hyperparameters
CLASS.NAME.COMPUTATIONAL.PARAMETERS.DEFINITION = 'mcleod.computational.parameters.obj' #Object holding definitions for the MCMC algorithm
CLASS.NAME.COVARIATES.ESTIMATION.PARAMETERS.DEFINITION = 'mcleod.covariates.estimation.parameters.obj' #object holding parameters to when covariates are also supplied to the deconvolition problem.

MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL = 1L
MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET = 0L

MCLEOD.BINOMIAL.ERRORS = 0L
MCLEOD.POISSON.ERRORS = 1L
MCLEOD.NORMAL.MEAN.IS.VAR.ERRORS = 2L


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constructors for auxiliary objects:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Generates an object of type 'mcleod.prior.def.obj' defining the prior for the random effect in the functions \code{mcleod}, \code{mcleod.estimate.CI}, \code{mcleod.estimate.CI.single.q} and \code{mcleod.estimate.CI.single.theta}
#'
#' Details the type of prior for mixing distribution. For a \code{\link{mcleod}} model with covariates, this function is used to define the distribution of the random intercept term. The prior for the distribution of random effects can be either Hierarchical Beta (a Polya Tree, i.e. a full binary tree with L levels, with each leaf corresponding to a segment on the real line, and each internal node associated with a beta random variable) or a two level dirichlet tree. See additional details in the package vignette.
#'
#' @details For a prior of type \code{MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL}, the Polya tree defines a piecewise constant prior with 2^(Beta.Heirarchical.Levels) segments on the real line. For priors defined using \code{MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET}, the number of segments on the real line is defined using \code{Two.Layer.Dirichlet.Intervals}, and the number of nodes in the upper later of the dirichlet tree is defined by \code{Two.Layer.Dirichlet.Nodes.in.First.Layer}. .The support for the prior is given by the parameter \code{a.limits} in the function \code{\link{mcleod}}, with "jumps" in the piecewise constant density function equally placed equally across the support defined by \code{a.limits}.
#' 
#' All intensity hyperparameters for either Beta variables (in the Polya tree), or gamma variables (in the 2-level Dirichlet tree) have a value of 1 by default. This can be changed using the parameters \code{Prior_Hyper_Parameters_BetaH_L}, \code{Prior_Hyper_Parameters_BetaH_U} and \code{Prior_Hyper_Parameters_2LDT}, see example in the package vignette.
#' 
#' @param prior.type The type of prior used. Values are either \code{MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL} or \code{MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET}
#' @param Beta.Heirarchical.Levels The number of levels used (for prior of type MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL). The first level (root of the tree) is 1. A full tree with L=4 levels has 2^(4-1) segments for the prior. 
#' @param Two.Layer.Dirichlet.Intervals When using \code{MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET} as the prior - the total number of leafs in the two-layer dirichet tree.
#' @param Two.Layer.Dirichlet.Nodes.in.First.Layer Number of nodes in the first level of the Dirichlet tree. \code{Two.Layer.Dirichlet.Intervals} must be an integer multiple of \code{Two.Layer.Dirichlet.Nodes.in.First.Layer}
#' @param Prior_Hyper_Parameters_BetaH_L A matrix of size (L, (2^(L-1))). Each row represents a level in the tree. Values in the lth row of the matrix, in entries 1 to 2^(l-1), correspond to the alpha_1 parameters of beta random variables, in the lth level of the polya tree.
#' @param Prior_Hyper_Parameters_BetaH_U matrix of size (L, (2^(L-1))). Each row represents a level in the tree. Values in the lth row of the matrix, in entries 1 to 2^(l-1), correspond to the alpha_2 parameters of beta random variables, in the lth level of the polya tree.
#' @param Prior_Hyper_Parameters_2LDT a matrix with dimensions (2,Two.Layer.Dirichlet.Intervals). The leftmost \code{Two.Layer.Dirichlet.Intervals} entries in the first row give the intensity parameters for the Dirichlet variable at the root of the tree. The values for the second row give the intensity parameters for the dirichlet random variables at the middle level: each \code{Two.Layer.Dirichlet.Intervals / Two.Layer.Dirichlet.Nodes.in.First.Layer} define the intensity parameters of a single Dirichlet random variable, from the leftmost to rightmost variables. See example on how to use in the package vignette.
#' @return An object of type \code{mcleod.prior.def.obj}
#' @export
#'
#' @examples
#' See the package vignette
mcleod.prior.parameters = function(prior.type = MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL,
                                   Beta.Heirarchical.Levels = 6,
                                   Two.Layer.Dirichlet.Intervals = 64,
                                   Two.Layer.Dirichlet.Nodes.in.First.Layer = 8,
                                   Prior_Hyper_Parameters_BetaH_L = NULL,
                                   Prior_Hyper_Parameters_BetaH_U = NULL,
                                   Prior_Hyper_Parameters_2LDT = NULL){
  
  #checks:
  if(!(prior.type%in% c(MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL,MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET))){
    stop('prior type must be MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL or MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET')
  }
  if(Beta.Heirarchical.Levels < 3){
    stop('Beta.Heirarchical.Levels must be at least 3') 
  }
  temp = Two.Layer.Dirichlet.Intervals/Two.Layer.Dirichlet.Nodes.in.First.Layer
  if(temp != as.integer(temp)){
    stop('Two.Layer.Dirichlet.Intervals must be a complete multiple of Two.Layer.Dirichlet.Nodes.in.First.Layer')
  }
  
  
  #if no priors are given, we generate default parameters:
  #for the BetaH case:
  if(prior.type == MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL & (is.null(Prior_Hyper_Parameters_BetaH_L) | is.null(Prior_Hyper_Parameters_BetaH_U)) ){
    Prior_Hyper_Parameters_BetaH_L = matrix(NA,nrow = Beta.Heirarchical.Levels,ncol = 2^(Beta.Heirarchical.Levels-1))
    Prior_Hyper_Parameters_BetaH_U = matrix(NA,nrow = Beta.Heirarchical.Levels,ncol = 2^(Beta.Heirarchical.Levels-1))
    for(l in 1:Beta.Heirarchical.Levels){
      for(u in 1:(2^(l-1))){
        Prior_Hyper_Parameters_BetaH_L[l,u] = 1
        Prior_Hyper_Parameters_BetaH_U[l,u] = 1
      }
    }
  }
  #for the 2 layer Dirichlet tree case:
  if(prior.type == MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET & is.null(Prior_Hyper_Parameters_2LDT)){
    Prior_Hyper_Parameters_2LDT = matrix(NA,nrow = 2,ncol = Two.Layer.Dirichlet.Intervals)
    for(u in 1:Two.Layer.Dirichlet.Nodes.in.First.Layer){
      Prior_Hyper_Parameters_2LDT[1,u] = 1
    }
    for(u in 1:Two.Layer.Dirichlet.Intervals){
      Prior_Hyper_Parameters_2LDT[2,u] = 1
    }
  }
  
  #check dimensions:
  if(prior.type == MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL){
    dimensions_should_be = c(Beta.Heirarchical.Levels,2^(Beta.Heirarchical.Levels-1))
    if(!all.equal(dim(Prior_Hyper_Parameters_BetaH_L),dimensions_should_be) | 
                  !all.equal(dim(Prior_Hyper_Parameters_BetaH_U),dimensions_should_be)){
      stop('Prior type selected to be Heirarchical Beta, but supplied values for Prior_Hyper_Parameters_BetaH_L or Prior_Hyper_Parameters_BetaH_U are not a (L,2^(L-1)) matrix')
    }
    Prior_Hyper_Parameters_2LDT = matrix(NA,nrow = 1,ncol = 1)
  }
  if(prior.type == MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET){
    dimensions_should_be = c(2,Two.Layer.Dirichlet.Intervals)
    if(!all.equal(dim(Prior_Hyper_Parameters_2LDT),dimensions_should_be)){
      stop('Prior type selected to be 2 Layer Dirichlet Tree, but supplied value for Prior_Hyper_Parameters_2LDT is not a (2,Two.Layer.Dirichlet.Intervals) matrix')
    }
    Prior_Hyper_Parameters_BetaH_L = matrix(NA,nrow = 1,ncol = 1)
    Prior_Hyper_Parameters_BetaH_U = matrix(NA,nrow = 1,ncol = 1)
  }
  
  
  ret = list()
  ret$prior.type = prior.type
  ret$Beta.Heirarchical.Levels = Beta.Heirarchical.Levels
  ret$Two.Layer.Dirichlet.Intervals = Two.Layer.Dirichlet.Intervals
  ret$Two.Layer.Dirichlet.Nodes.in.First.Layer = Two.Layer.Dirichlet.Nodes.in.First.Layer
  ret$Prior_Hyper_Parameters_BetaH_L = Prior_Hyper_Parameters_BetaH_L
  ret$Prior_Hyper_Parameters_BetaH_U = Prior_Hyper_Parameters_BetaH_U
  ret$Prior_Hyper_Parameters_2LDT = Prior_Hyper_Parameters_2LDT
  class(ret) = CLASS.NAME.PRIOR.DEFINITION
  return(ret)
}


#' Returns an object of type 'mcleod.computational.parameters.obj' used to set the computational parameters for the Gibbs sampler, used for estimating the random effect distribution.
#'
#' @param nr.gibbs Number of gibbs samples used for the MCMC chain
#' @param nr.gibbs.burnin Number of gibbs samples taken as burn-in (i.e., samples at the start of the chain that are disregarded when estimating model parameters)
#' @param integration_step_size Used for numerical integration when computing the probability for an observed value to originate from a given segment of the underlying density (specifically, from the segments of a.vec, the segments of the piecewise constant density function for the prior)
#' @param Fast.Gamma.Used Should a "fast-and-inaccurate" sampling procedure from a gamma density be used. Not recommended and disabled by default.
#' @param Fast.Gamma.Bank.Size Size of pregenerated bank of gamma variables: for each rate (0.1,1,10), this number of random deviates is sampled from each rate. Gamma variables are generated by sampling from the pregenerated banks, for any custom rate, by sampling deviates from each "rate-bank" and summing the results up to the required rate.
#'
#' @return An object of size 'mcleod.computational.parameters.obj'
#' @export
#'
#' @examples See the package vignette
mcleod.computational.parameters = function(nr.gibbs = 500,
                                                   nr.gibbs.burnin = 250,
                                                   integration_step_size = 0.01,
                                                   Fast.Gamma.Used = F,
                                                   Fast.Gamma.Bank.Size = 1000L){
  #checks:
  #nr.gibbs > nr.gibbs burnin
  if(nr.gibbs.burnin >= nr.gibbs){
    stop('nr.gibbs must be strictly larger than nr.gibbs.burnin')
  }
  #warnings nr.gibbs and burnins
  if(nr.gibbs.burnin < 200)
    warning('warning: nr.gibbs.burnin lower than 200 iterations')
  if(nr.gibbs < 500)
    warning('warning: nr.gibbs lower than 500 iterations')
  #warning on integration_step_size
  if(integration_step_size>0.01)
    warning('warning: integration_step_size in natural parameter scale larger than 0.01')
  
  if(Fast.Gamma.Bank.Size < 1000)
    warning('warning: Fast.Gamma.Bank.Size smaller than 1000')
  
  ret = list()
  ret$nr.gibbs = nr.gibbs
  ret$nr.gibbs.burnin = nr.gibbs.burnin
  ret$integration_step_size = integration_step_size
  ret$Fast.Gamma.Used = Fast.Gamma.Used
  ret$Fast.Gamma.Bank.Size = Fast.Gamma.Bank.Size
  class(ret) = CLASS.NAME.COMPUTATIONAL.PARAMETERS.DEFINITION
  return(ret)
}


#' Creates object defining computational and statistical parameters, when estimating a random effect model with covariates.
#' 
#' For binomial response, the model with covariates is \eqn{log(P_i/(1-P_i)) = \gamma_i + \vec{\beta}^T \vec{Z}_i}, where the density of \eqn{\gamma_i}'s is estimated from the data, \eqn{\vec{\beta}} is a vector of slopes, and \eqn{\vec{Z}_i} is a vector of covariates. For Poisson errors, the model equation is \eqn{log(\lambda_i) = \gamma_i + \vec{\beta}^T \vec{Z}_i}. The function is used to create an object the defines the prior distribution of \eqn{\vec{\beta}}, as well as the proposal distribution used for drawing new values of \eqn{\vec{\beta}} in the MCMC algorithm.
#' 
#' @details The default prior assumed for \eqn{\vec{\beta}} is a multivariate normal prior, with mean of 0, and a covariance matrix given by \code{diag(beta_prior_sd)}. Setting one of the standard deviations to a negative value means a non-informative prior is assumed for the slope coefficient.
#'
#' Alternatively, the user may define a non normal prior for each slope coefficient. This is done by supplying the the function with a series piece-wise constant density functions serving as priors for the different slope coefficients.
#' A piecewise constant prior density function with M segments has M+1 partitioning points. These are supplied using the argument \code{Manual_Prior_Values}. The density values for each segment are given by \code{Manual_Prior_Probs}, which has M rows (one for each segment), and number of columns equal to the number of covariates.
#' The parameter \code{beta_init} sets the initial values for \eqn{\vec{\beta}} in the MCMC sampler. The default value is \code{NULL}, meaning a normal random intercept regression will be used to guess an initial estimate. The type of regression used is either Binomial or Poisson, depending on the sampling distribution defined when calling \code{\link{mcleod}}.
#' 
#' @param proposal_sd vector of standard deviations for the proposal distribution, for the values of beta. Proposals are generated by adding a random value from N(0,diag(proposal_sd)) to the current MCMC value of beta.
#' @param beta_prior_sd vector of standard deviations for the prior. The prior for beta is assumed to be N(0,diag(beta_prior_sd)). Setting an entry to a negative value (e.g. -1) will use an uninformative prior (constant) for that component of beta.
#' @param beta_init Initial values for \eqn{\vec{\beta}} in the MCMC sampler. If set to \code{NULL}, will use a normal random intercept regression for initial estimate.
#' @param Manual_Prior_Values When using a predefined prior for beta, a vector of sorted values defining the paritioning points for a piece-wise constant prior for beta.
#' @param Manual_Prior_Probs A matrix with number of rows equal to = length(Manual_Prior_Values)-1, and number of columns equal to the number of covariates, giving the probablities for the different segments defined by  \code{Manual_Prior_Values}. Within each segment, it is assumed the prior is uniformly distributed.
#' @param do_P_k_i_hashing A computational workaround for faster computation of bayesian probabilites. See additional details in the appendix describing the computational method, in the dissertation of Brill (2022)
#' @param P_k_i_hashing_resolution_by_theta    Resolution parameter for the computational workaround activated by do_P_k_i_hashing
#'
#' @return an object of type 'mcleod.covariates.estimation.parameters.obj'
#' @export
#'
#' @examples See the package vignette
mcleod.covariates.estimation.parameters = function(proposal_sd = c(0.05),
                                           beta_prior_sd = c(1),
                                           beta_init = NULL,
                                           Manual_Prior_Values = NULL,
                                           Manual_Prior_Probs = NULL,
                                           do_P_k_i_hashing = F,
                                           P_k_i_hashing_resolution_by_theta = 0.0001){
  
  #checks:
  if(any(beta_prior_sd <=0))
     warning('not all entries of beta_prior_sd are strictly larger than zero, non positive entries will be given a non-informative prior in sampling')
  if(any(proposal_sd <=0))
    stop('all entries of proposal_sd must be strictly larger than zero')
  
  if(!is.null(Manual_Prior_Values)){
    if(length(Manual_Prior_Values) - 1 != dim(Manual_Prior_Probs)[1]){
      stop('Manual_Prior_Values must be longer by one entry than nr rows of Manual_Prior_Probs')      
    }
    if(!all.equal(Manual_Prior_Values,sort(Manual_Prior_Values))){
      stop('Manual_Prior_Values must be sorted')
    }
    if(any(Manual_Prior_Probs<0)){
      stop('All values of Manual_Prior_Probs must be positive')
    }
    for(j in 1:ncol(Manual_Prior_Probs)){
      if(abs(sum(Manual_Prior_Probs[,j]) - 1.0)> 0.00001){
        stop(paste0('probabilities in column ',j,' of Manual_Prior_Probs do not sum up to 1'))
      }  
    }
    
  }
  
  ret = list()
  ret$proposal_sd = proposal_sd
  ret$beta_prior_sd = beta_prior_sd
  ret$beta_init = beta_init
  ret$Manual_Prior_Values = Manual_Prior_Values
  ret$Manual_Prior_Probs = Manual_Prior_Probs
  ret$do_P_k_i_hashing = do_P_k_i_hashing
  ret$P_k_i_hashing_resolution_by_theta = P_k_i_hashing_resolution_by_theta
  class(ret) = CLASS.NAME.COVARIATES.ESTIMATION.PARAMETERS.DEFINITION
  return(ret)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# main function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Function estimates the mixing distribution, for data with Binomial/Poisson samples with varying success probabilities/rates.
#' 
#' Given a dataset with either binomial or Poisson samples, the function estimates the mixing distribution of \eqn{P} or \eqn{\lambda}{lambda}, respectively. The function assumes the following hierarchical Bayesian model for the mixing distribution: For binomial measurements, the model is \eqn{X_i\sim bin(N,P_i)}{Xi~bin(N,Pi)}, where  \eqn{log(P_i/(1-P_i))}{log(Pi/(1-Pi))} is random deviate from a piecewise-constant function, whose density values are modeled using either a Polya-Tree or Dirichlet Tree distribution. For Poisson errors, the model equation is \eqn{X_i\sim Pois(lambda_i)}{X_i~Pois(lambda_i)}, where \eqn{log(lambda_i)} is modeled using the heirarchical Bayesian approach (using a piecewise-constant density function).
#' If covariates are included in the model, the model is \eqn{log(P_i/(1-P_i)) = \gamma_i + \vec{\beta}^T \vec{Z}_i}, where the density of \eqn{\gamma_i}'s is estimated using the heirarchical Bayesian approach, \eqn{\vec{\beta}} is a vector of slopes, and \eqn{\vec{Z}_i} is a vector of sample covariates. For Poisson errors, the model equation is \eqn{log(\lambda_i) = \gamma_i + \vec{\beta}^T \vec{Z}_i}.
#' 
#' @details The mixing distribution ( or the density of \eqn{\gamma_i}'s for the setting with covariates) is estimated using an MCMC sampler conditioning on the observed data. The hyperparameters for prior for the mixing distribution are defined using the function \code{\link{mcleod.prior.parameters}}, and passed as an object to the argument \code{prior_parameters}, together with the support for the mixing distribution (random intercept term for the case with covariates), and the number of breaks in the piecewise constant density function (discontinuity points in the density function are equally placed across the support).
#' The returned object contains both the parameters and data given by the user, together with results for the MCMC sampler. The main result returned by the function is \code{pi_smp}: a series of density functions sampled from the posterior distribution of the mixing distribution (density of \eqn{\gamma_i}'s for the case with covariates). The function \code{\link{plot.posterior}} plots the samples from the posterior distribution of the random effect. See the package vignette on how to extract results from the returned object.
#' See the package vignette and cited papers for additional information on the model and statistical estimation approach.
#' 
#' If covariates are also supplied by the user, the algorithm performs two types of MCMC steps: one for the density of the random intercept term (density of  \eqn{\gamma_i's}) and one for the vector of slope coefficients (given by \eqn{\vec{\beta}}). MCMC steps for the slope coefficients are done using an accept-reject Metropolis-Hastings step. The function \code{\link{results.covariate.coefficients.posterior}} allows the mean posterior estiamtes of \eqn{\vec{\beta}} to be extracted from the output, together with plots and statistics for the posterior distribution of \eqn{\vec{\beta}}, and the acceptance ratio of the MH step.
#' The prior distribution for \eqn{\vec{\beta}}, together with parameters for the proposal distribution (suggesting the MCMC algorithm new values of \eqn{\vec{\beta}} based on the current iteration's values) can be set using the parameter \code{covariates_estimation_parameters}.
#'
#' @param x.smp a vector of N measurements. For the binomial sampling distribution, measurement must be smaller or equal to n.smp (entry-wise)
#' @param n.smp number of draws, per observation, for the binomial model. For poisson sampling errors, set this value to NULL.
#' @param a.limits a vector with two entries, setting the range for the support of log(p/(1-p)) or log(lambda) for the model with no covariates. When using a model with covariates, these two values set the range for the support of gamma.
#' @param Noise_Type Noise type: either \code{MCLEOD.BINOMIAL.ERRORS} or \code{MCLEOD.POISSON.ERRORS}
#' @param covariates (optional) a matrix of covariates. Must as have a number of rows equal to length(x.smp). Different colmn represent different covariates.
#' @param prior_parameters Result of \code{\link{mcleod.prior.parameters}}. Sets parameter for the mixing / random intercept distribution.
#' @param computational_parameters Result of \code{\link{mcleod.computational.parameters}} Sets computational parameters for the MCMC sampler used for estimating the mixing / random intercept distribution.
#' @param covariates_estimation_parameters Result of \code{\link{mcleod.covariates.estimation.parameters}}. Sets statistical and computational 
#' @param input_P_k_i A parallel approach for specifying the sampling distribution (instead of using binomial/poisson sampling). This parameter is used to define a matrix specifying each observations probability of being sampled, from each segment of the prior. The rows of the matrix need to correspond to observations, and the columns need to be associated with the segments of the piece-wise constant mixing distribution. The number of segments is defined by the prior object given in the argument \code{prior_parameters}. For example, for the Heirarchical Beta (Polya tree) prior there are 2^(L-1) segments placed evenly in the range defined by the argument \code{a.limits}, see additional details in the function \code{\link{mcleod.prior.parameters}}
#' @param exact.numeric.integration should exact numeric integration be used. See additional details in the doc page for \code{\link{mcleod.prior.parameters}}.
#' @param offset_vec A vector of covariates added with a slope of 1 to the model. This can be used for example to add terms of the form "+log(N_i)" to a poisson model, in order to incorporate knowledge about the extensive size of the sampled units.
#' @param nr_threads The function allows multiple x.smp to be a list of data samples, with each entry of the list being a vector of the same size as \code{n.smp}. When x.smp is set to be a list of vectors, this parameter sets the number of threads used for parallel computation of the deconvolution estimates. There results are returned as a list: \code{original_stat_res} (described below) become a list of the deconvolution estimates for the different data samples (with entries corresponding to the entries of the list supplied under \code{x.vec}). If \code{x.smp} is a list, then \code{input_P_k_i} can also be a list (but not vice versa)
#'
#' @return an object of type 'mcleod.obj' with the two following lists.
#' A list named \code{parameters_list} with the following entries:
#' \itemize{
#' \item{a.vec}{ - The mixing distribution for the binomial/poisson samples is estimated as a piece-wise constant function, over these breaking point.}
#' \item{nr.gibbs}{ - Number of gibbs samples}
#' \item{nr.gibbs.burnin}{ - Number of gibbs samples taken as burnin (excluded when computing estimators)}
#' \item{prior_parameters}{ - Input object by same name.}
#' \item{computational_parameters}{ - Input object by same name.}
#' \item{covariates_estimation_parameters}{ - Input object by same name.}
#' \item{covariates}{ - The matrix of inserted covariates}
#' \item{x.smp}{ - Input data by same name}
#' \item{n.smp}{ - Input data by same name}
#' \item{Noise_Type}{ - 0 for binomial, 1 for Poisson.}
#' \item{covariates_given}{ - Logical value indicating if covariates were given.}
#' }
#'
#' And a second list named \code{original_stat_res} with the following entries:
#' \itemize{
#' \item{p_k_i}{ - A matrix with dimensions (length(x.smp),length(a.vec)-1). Entries of the matrix are the probabilities of each observation (by row), to be sampled when the random effect is sampled uniformly from each of the segments of a.vec (segments corresponding to columns).}
#' \item{n_smp}{ - A matrix of size (length(a.vec)-1,nr.gibbs) giving for each MCMC iteration (by column), the number of observations generated from each segment of the support of the mixing distribution. }
#' \item{pi_smp}{ - A matrix of size (length(a.vec)-1,nr.gibbs) giving for each MCMC iteration (by column), a posterior sample of the mixing distribution. computing the mean of each row gives the posterior mean for each segment of the piecewise constant mixing distribution. }
#' \item{beta_smp}{ - When covariates are included, a matrix of size (ncol(covariates),nr.gibbs), giving the MCMC samples for the vector of slope coefficients.}
#' \item{beta_suggestion}{ - When covariates are included,a matrix of size (ncol(covariates),nr.gibbs), giving the MCMC proposals for the vector of slope coffeficients. At each iteration, the proposal is suggested based on a normally distributed step from the previous MCMC step, and chosen usin a Metropolis Hastings rule.}
#' \item{proposal_approved}{ - When covariates are included: for each step, was the proposal for the vector of slope coefficients approved. The parameters proposal_approved, ll_proposal and ll_current are used for interogating the acceptance rate for the sampler, when covariates are used.}
#' \item{elapsed_secs}{ - Running time in seconds.}
#' \item{ll_proposal}{ - When covariates are included: for each MCMC step, the loglikelihood value of proposal of the new parameter values. Values computed without the probability of the polya tree prior, see CPP code for details.}
#' \item{ll_current}{ - When covariates are included: for each MCMC step, the loglikelihood value of the current parameter values. Values computed without the probability of the polya tree prior, see CPP code for details.}
#' }
#' @export
#'
#' @examples
#' see package vignette
mcleod	<- function( x.smp,
                     n.smp,
                     a.limits = c(-4,4),
                     Noise_Type = MCLEOD.BINOMIAL.ERRORS,
                     covariates = NULL,
                     prior_parameters = NULL,
                     computational_parameters = NULL,
                     covariates_estimation_parameters = NULL,
                     input_P_k_i = NULL,
                     exact.numeric.integration = TRUE,
                     offset_vec = rep(0,length(x.smp)),
                     nr_threads = 1
                     )
{
  
#  exact.numeric.integration = TRUE # We force exact numeric integration over dbinom for computation of P_k_i. Normal approximation is not sufficiant.
  
  
  #%%% Retreive prior parameters
  
  
  if(is.null(prior_parameters)){
    prior_parameters = mcleod.prior.parameters()
  }
  #check object type
  if(class(prior_parameters) != CLASS.NAME.PRIOR.DEFINITION){
    stop('argument prior_parameters must be generated by function mcleod.prior.parameters(...)')
  }
  L = prior_parameters$Beta.Heirarchical.Levels
  I1 = prior_parameters$Two.Layer.Dirichlet.Nodes.in.First.Layer
  Prior_Type = prior_parameters$prior.type
  I_specificy_parameter = prior_parameters$Two.Layer.Dirichlet.Intervals
  
  #Hyper parameters for the prior:
  Prior_Hyper_Parameters_BetaH_L = prior_parameters$Prior_Hyper_Parameters_BetaH_L
  Prior_Hyper_Parameters_BetaH_U = prior_parameters$Prior_Hyper_Parameters_BetaH_U
  Prior_Hyper_Parameters_2LDT = prior_parameters$Prior_Hyper_Parameters_2LDT
  
  
  #%%% Retreive computational parameters
  if(is.null(computational_parameters)){
    computational_parameters = mcleod.computational.parameters()
  }
  #check object type
  if(class(computational_parameters) != CLASS.NAME.COMPUTATIONAL.PARAMETERS.DEFINITION){
    stop('argument computational_parameters must be generated by function mcleod.computational.parameters(...)')
  }
  nr.gibbs = computational_parameters$nr.gibbs
  nr.gibbs.burnin = computational_parameters$nr.gibbs.burnin
  integration_step_size = computational_parameters$integration_step_size
  Fast.Gamma.Used = computational_parameters$Fast.Gamma.Used
  Fast.Gamma.Bank.Size = computational_parameters$Fast.Gamma.Bank.Size
  
  #%%% Retreive covariate estimation parameters
  if(is.null(covariates_estimation_parameters)){
    covariates_estimation_parameters = mcleod.covariates.estimation.parameters()
  }
  if(class(covariates_estimation_parameters) != CLASS.NAME.COVARIATES.ESTIMATION.PARAMETERS.DEFINITION){
    stop('argument covariates_estimation_parameters must be generated by function mcleod.covariates.estimation.parameters(...)')
  }
  #check object type
  proposal_sd = covariates_estimation_parameters$proposal_sd
  beta_prior_sd = covariates_estimation_parameters$beta_prior_sd
  beta_init = covariates_estimation_parameters$beta_init
  
  #%%% Retreive covariates:
  if(is.null(covariates)){
    covariates = matrix(c(1),nrow = 1)
    covariates_given = 0
  }else{
    covariates_given = 1
  }
  
  # adjust length of covariate estimation parameters, if they are given as a single value. Else, throw an error
  if(length(proposal_sd) != ncol(covariates)){
    if(length(proposal_sd) == 1)
      proposal_sd = rep(proposal_sd,ncol(covariates))
    else
      stop('proposal_sd must be of length ncol(covariates)')
  }
  if(length(beta_prior_sd) != ncol(covariates)){
    if(length(beta_prior_sd) == 1)
      beta_prior_sd = rep(beta_prior_sd,ncol(covariates))
    else
      stop('beta_prior_sd must be of length ncol(covariates)')
  }
  
  if(!is.null(beta_init))
    if(length(beta_init) != ncol(covariates)){
      if(length(beta_init) == 1)
        beta_init = rep(beta_init,ncol(covariates))
      else
        stop('beta_init must be of length ncol(covariates)')
    }
  
  if(!all(is.numeric(offset_vec)) & length(offset_vec)!=length(n.smp)){
    stop('offset_vec must be numric and same length as n.smp')
  }
  
  #%%% Function settings
  
  #handle if input_P_k_i is given. Also if given, is it a matrix.
  if(is.null(input_P_k_i)){
    P_k_i_is_given = 0L
    temp_matrix = matrix(1,nrow = 1)
    if(!is.list(x.smp)){
      P_k_i_precomputed = temp_matrix
    }else{
      P_k_i_precomputed = list()
      for(i in 1:length(x.smp)){
        P_k_i_precomputed[[i]] = temp_matrix
      }
    }
    K			<- length(x.smp)
  }else{
    P_k_i_is_given = 1L
    P_k_i_precomputed = input_P_k_i
    if(is.list(P_k_i_precomputed)){
      K			<- nrow(P_k_i_precomputed[[1]])
    }else{
      K			<- nrow(P_k_i_precomputed)  
    }
    
  }
  
  #checks
  #Noise type legal
  if(!(Noise_Type %in% c(MCLEOD.BINOMIAL.ERRORS,MCLEOD.POISSON.ERRORS))){
    stop('Noise_Type must be either set to MCLEOD.BINOMIAL.ERRORS or MCLEOD.POISSON.ERRORS')
  }
  if(!P_k_i_is_given & Noise_Type == MCLEOD.BINOMIAL.ERRORS){
    if(is.list(x.smp)){
      length_x_smp = length(x.smp[[1]])
    }else{
      length_x_smp = length(x.smp)
    }
   if(length_x_smp != length(n.smp)){
     stop('For binomial errors, x.smp and n.smp must be of equal length')
   } 
  }
  if(!P_k_i_is_given & Noise_Type == MCLEOD.POISSON.ERRORS){
    if(!is.null(n.smp)){
      stop('For poisson errors, n.smp must be set to NULL')
    }
    if(is.list(x.smp)){
      n.smp = x.smp[[1]]  
    }else{
      n.smp = x.smp
    }
    
  }
  if(a.limits[2]<=a.limits[1]){
    stop('a.limits must be vector of size 2: c(lower,upper) for natural parameter')
  }
  
  
  
  #for poisson - check if a.limits are in appropriate
  
  
  a.max = a.limits[2]
  a.min = a.limits[1]
  
  
  Fast.Gamma.Bank = matrix(1,nrow = 1)
  Fast.Gamma.Used.p = 0
  if(Fast.Gamma.Used){
    Fast.Gamma.Used.p = 1
    Fast.Gamma.Bank = rcpp_Generate_Fast_Gamma_Bank(Fast.Gamma.Bank.Size)
  }
  
  #compute the number of segments for the prior, based on the prior type
  I = 2^(L)
  if(Prior_Type == MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET){
    if(!is.null(I_specificy_parameter))
        I = I_specificy_parameter
    if(I/I1 != as.integer(I/I1)){
     stop('I/I1 not an integer!') 
    }
  }
  
  #the partitioning point for the prior
  a.vec.used		<- seq(a.min,a.max,length = I+1)  
  
  if(P_k_i_is_given){
    if(!is.list(P_k_i_precomputed)){
      #check precomputed P_k_i is same size as a.vec - this can be done in Wrapper_rcpp_Gibbs
      if(ncol(P_k_i_precomputed) != I){
        stop(paste0(' P_k_i contains ',ncol(P_k_i_precomputed),' columns, but prior definitions has ',I,' intervals'))
      }
      if(nrow(P_k_i_precomputed) != K){
        stop(paste0(' P_k_i contains ',nrow(P_k_i_precomputed),' rows, but x.smp is of length ',K))
      }
      #check that if Precomputed P_k_i - there are no covariates given!
      if(P_k_i_is_given & covariates_given){
        stop('covariates can be given only for Binomial or Poisson errors. A custom P_k_i must be given without covariates')
      }
    }else{
      #need to implement check for list precomputed P_k_i
    }
    
  }
  
  #handle manual priors:
  Manual_Prior_Given = c(0L)
  Manual_Prior_Values = c(-4,4)
  Manual_Prior_Probs = matrix(c(1))
  if(!is.null(covariates_estimation_parameters$Manual_Prior_Values)){
    Manual_Prior_Given = c(1L)
    Manual_Prior_Values = covariates_estimation_parameters$Manual_Prior_Values
    Manual_Prior_Probs = covariates_estimation_parameters$Manual_Prior_Probs
    if(ncol(covariates) != dim(Manual_Prior_Probs)[2]){
      stop('ncols of  Manual_Prior_Probs must be the same as the number of covariates')
    }
    
  }
  
  
  if(covariates_given & is.null(beta_init)){
    if(P_k_i_is_given){
      beta_init = rep(0,ncol(covariates))
    }else{
      x_to_pass = x.smp
      if(is.list(x.smp)){
        x_to_pass = x.smp[[1]]
      }
      beta_init = init_mcleod_random_intercept_regression(x = x_to_pass,
                                                          n = n.smp,
                                                          covariates = covariates,
                                                          offset_p = data.frame('offset' = offset_vec),
                                                          family = ifelse(Noise_Type == MCLEOD.BINOMIAL.ERRORS,'binomial','poisson'))
    }
  }
  if(!covariates_given){
    beta_init = c(0)
  }  
  if(!is.list(x.smp)){
    #%%% call Rcpp wrapper, for the case x.smp IS NOT a list
    res = Wrapper_rcpp_Gibbs(x.smp,
                             n.smp,
                             a.vec.used,
                             nr.gibbs,
                             nr.gibbs.burnin,
                             as.integer(exact.numeric.integration),
                             as.integer(0), #verbose - turned off
                             L,
                             Prior_Hyper_Parameters_BetaH_L,
                             Prior_Hyper_Parameters_BetaH_U,
                             Prior_Hyper_Parameters_2LDT,
                             Fast.Gamma.Used.p,
                             Fast.Gamma.Bank,
                             PriorType = Prior_Type,
                             I1 = I1,
                             covariates_given = covariates_given,
                             covariates = covariates,
                             proposal_sd = proposal_sd,
                             beta_prior_sd = beta_prior_sd,
                             beta_init = beta_init,
                             integration_step_size = integration_step_size,
                             Noise_Type = Noise_Type,
                             P_k_i_is_given = P_k_i_is_given,
                             P_k_i_precomputed = P_k_i_precomputed,
                             Manual_Prior_Given = Manual_Prior_Given,
                             Manual_Prior_Values = Manual_Prior_Values,
                             Manual_Prior_Probs = Manual_Prior_Probs,
                             do_P_k_i_hashing = covariates_estimation_parameters$do_P_k_i_hashing,
                             P_k_i_hashing_resolution = covariates_estimation_parameters$P_k_i_hashing_resolution_by_theta,
                             offset_vec = offset_vec
    )
  }else{
    #%%% call Rcpp wrapper,  for the case x.smp IS not a list
    res = Wrapper_rcpp_Gibbs_list(x.smp,
                             n.smp,
                             a.vec.used,
                             nr.gibbs,
                             nr.gibbs.burnin,
                             as.integer(exact.numeric.integration),
                             as.integer(0), #verbose - turned off
                             L,
                             Prior_Hyper_Parameters_BetaH_L,
                             Prior_Hyper_Parameters_BetaH_U,
                             Prior_Hyper_Parameters_2LDT,
                             Fast.Gamma.Used.p,
                             Fast.Gamma.Bank,
                             PriorType = Prior_Type,
                             I1 = I1,
                             covariates_given = covariates_given,
                             covariates = covariates,
                             proposal_sd = proposal_sd,
                             beta_prior_sd = beta_prior_sd,
                             beta_init = beta_init,
                             integration_step_size = integration_step_size,
                             Noise_Type = Noise_Type,
                             P_k_i_is_given = P_k_i_is_given,
                             P_k_i_precomputed = P_k_i_precomputed,
                             Manual_Prior_Given = Manual_Prior_Given,
                             Manual_Prior_Values = Manual_Prior_Values,
                             Manual_Prior_Probs = Manual_Prior_Probs,
                             do_P_k_i_hashing = covariates_estimation_parameters$do_P_k_i_hashing,
                             P_k_i_hashing_resolution = covariates_estimation_parameters$P_k_i_hashing_resolution_by_theta,
                             offset_vec = offset_vec
    )
  }
  
  
  #%%% Wrap results
  
  ret = list()
  class(ret) = CLASS.NAME.MCLEOD
  if(is.list(x.smp))
    class(ret) = paste0(class(ret),'_MULTIPLE') # this is so print functions wont send an error
  
  ret$parameters_list = list(
    a.vec = a.vec.used,
    nr.gibbs = nr.gibbs,
    nr.gibbs.burnin = nr.gibbs.burnin,
    prior_parameters = prior_parameters,
    computational_parameters = computational_parameters,
    covariates_estimation_parameters = covariates_estimation_parameters,
    covariates = covariates,
    x.smp = x.smp,
    n.smp = n.smp,
    Noise_Type = Noise_Type,
    covariates_given = covariates_given
  )
  
  # note that original_stat_res is a list itself, if x.smp was a list
  ret$additional = list(original_stat_res = res)
  return(ret)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions for plotting results
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Plot the posterior mixing distribution, estimated for the data from a mcleod model
#' 
#' The function receives an object returned from \code{\link{mcleod}} and plots the posterior estimate for the mixing distribution. The posterior mean is shown in red (means for each point in the support of the mixing distribution). A point cloud shows the posterior distribution at each point of the support.
#' 
#' @param mcleod.obj Object received from \code{\link{mcleod}}.
#' @param plot_only_point_estimate Logical value indicating if only the posterior mean should be plotted (using a red line). If set to F (default value), the posterior density will also be plotted.
#'
#' @return NONE
#' @export
#'
#' @examples
#' See package vignette
plot.posterior	<- function(mcleod.obj, plot_only_point_estimate = F)
{
  library(ggplot2)
 if(class(mcleod.obj) != CLASS.NAME.MCLEOD){
   stop('input argument for plot.posterior must be an object returned from function mcleod')
 }

  burnin = mcleod.obj$parameters_list$nr.gibbs.burnin
  gibbs_pi_plot = cbind(0,t(apply(t(mcleod.obj$additional$original_stat_res$pi_smp),1,cumsum)))
  gibbs_pi_plot = gibbs_pi_plot[-c(1:burnin),]
  means_vec = apply(gibbs_pi_plot,2,mean)
  median_vec = apply(gibbs_pi_plot,2,median)
  
  #construct the line of means, and the point cloud for the posterior distribution at each point of a.vec
  a.vec = mcleod.obj$parameters_list$a.vec
  
  dt_mean_line = data.frame(a.vec = a.vec,means_vec = means_vec,median_vec = median_vec)
  
  gibbs_cloud_size = dim(gibbs_pi_plot)[1] * dim(gibbs_pi_plot)[2]
  dt_gibbs_cloud = data.frame(a.point = rep(NA,gibbs_cloud_size), CDF.value = rep(NA,gibbs_cloud_size))
  pointer = 1
  for(i in 1:(dim(gibbs_pi_plot)[1])){
    for(j in 1:(dim(gibbs_pi_plot)[2])){
      dt_gibbs_cloud[pointer,] = c(a.vec[j],gibbs_pi_plot[i,j]); pointer = pointer + 1
    }
  }
  
  #perform the actual plots. mean line and then point (if needed)
  gg_obj = ggplot(dt_gibbs_cloud)  + ylim(c(0,1)) +
    geom_line(aes(x = a.vec,y = means_vec),colour = 'red',data = dt_mean_line) + xlab('theta')+ylab('CDF')
  
  if(!plot_only_point_estimate){
    gg_obj = gg_obj + geom_point(aes(x = a.point,y = CDF.value),alpha = 0.25,colour = 'gray',data = dt_gibbs_cloud,size = 0.8, shape = 18)
  }
  
  return(gg_obj)

}


#' Obtain MCMC samples. posterior means, and graphs for the slope coefficients, for the mcleod model with covariates
#'
#' The function receives the output of \code{\link{mcleod}}, when trying to estimate the mixing distribution with additional covariates, and performs the following actions: 1) returns the posterior means for the slope coefficients, and the mean acceptance rate for the MCMC sampler; 2) plots a graph of the values of the slope coefficients, by MCMC iteration number (see additional details for the function parameters); 3) plot the suggested vectors of slope coefficients, sampled by the MCMC algorithm at each iteration, and whether the suggested vector was accepted in each iteration (see additional details below); and 4) returns the MCMC samples and the MCMC suggestions for the slope coefficients (toghether with an indicator if the suggested value for \eqn{\vec{\beta}} was accepted in each iteration).
#'
#' @param mcleod.obj Object received from \code{\link{mcleod}}, when covariates are included in the model.
#' @param plot.posterior If set to True (default is True): plot a series of figures showing the values of the slope coefficients by MCMC iteration number. Burnin iterations (excluded from the computation of posterior means) are shown in red.
#' @param plot.MH.proposal.by.iteration If set to True (default is False): plot a series of figures showing the values of the slope coefficients for the slope vectors sampled as suggestions by the MCMC iteration number. Accepted iterations are shown in red.
#' @param aggregate_by the function by which to aggregate by. Pass the function \code{mean} (in R base) to compute the posterior mean of slope coefficients. Pass the function \code{median} (also in R base) to compute the posterior median.
#' @return Returns a list with the following entries:
#' \itemize{
#' \item{posterior.means}{ - a vector with estimates for the posterior means of the slope coefficients (see model description in \code{\link{mcleod}})}
#' \item{acceptance.rate}{ - the ratio of times proposals for slope vectors were accepted in the MCMC sampler.}
#' \item{beta_smp}{ - a matrix of size (nr MCMC samples, nr covariates), giving the MCMC samples for \eqn{\beta}}
#' \item{beta_suggestion}{ - a matrix of size (nr MCMC samples, nr covariates), giving the values sampled as suggestions for \eqn{\beta} in each sample of the algorithm.}
#' \item{proposal_approved}{ - a vector detailing for each MCMC sample if the suggestion sampled for \eqn{\beta} (in the MH step, see vignette and paper) was accepted (value = 1) or rejected (value = 0).}
#' }
#' @export
#'
#' @examples
#' See package vignette
results.covariate.coefficients.posterior = function(mcleod.obj, plot.posterior = T, plot.MH.proposal.by.iteration = F,aggregate_by = mean){

  if(class(mcleod.obj) != CLASS.NAME.MCLEOD){
    stop('input argument for plot.covariate.coefficients.posterior must be an object returned from function mcleod')
  }
  if(mcleod.obj$parameters_list$covariates_given != 1L){
    stop('covariates were not given, no coefficients estimated')
  }
  #extract the values aof covariates, their number, and the number of burnin samples
  covariates = mcleod.obj$parameters_list$covariates
  burnin = mcleod.obj$parameters_list$nr.gibbs.burnin
  nr.covariates = ncol(covariates)
  mean_vec = rep(NA,nr.covariates)
  
  #if we are required, plot for each covariate, in a seperate frame, the distribution of gibbs samples.
  #Gibbs samples in the burinin are marked in red
  
  if(plot.posterior){
    col_vec = rep(1,mcleod.obj$parameters_list$nr.gibbs)
    col_vec[1:burnin] = 2
    par(mfrow=c(nr.covariates,1))  
    for(i in 1:nr.covariates){
      plot(mcleod.obj$additional$original_stat_res$beta_smp[i,],col = col_vec,pch= 20,main = paste0('Coefficient for covariate ',i),xlab = 'Iteration',ylab = 'Coefficient')
    }  
    par(mfrow=c(1,1))
  }
  
  #compute the posterior means, after exclusing the burnin
  posterior_mean_vec = apply(mcleod.obj$additional$original_stat_res$beta_smp[,-c(1:burnin),drop=F],1,aggregate_by)
  
  # plot for each covariate, in a seperate frame, the distribution of suggested proposals. Accepted proposals are marked in red (note that a proposal is for all slope variables simultanuesly)
  if(plot.MH.proposal.by.iteration){
    par(mfrow=c(nr.covariates,1))  
    for(i in 1:nr.covariates){
      #plot, remove proposal for the last iteration
      plot(mcleod.obj$additional$original_stat_res$beta_suggestion[i, -mcleod.obj$parameters_list$nr.gibbs ],
           col = mcleod.obj$additional$original_stat_res$proposal_approved+1,pch= 20,main = paste0('beta',i,' - proposal'),xlab = 'Iteration',ylab = 'Coefficient')    
    }
    par(mfrow=c(1,1))
  }
  
  #return the posterior means and the acceptance rate
  ret = list()
  ret$posterior.means = posterior_mean_vec
  ret$acceptance.rate = mean(mcleod.obj$additional$original_stat_res$proposal_approved)
  ret$beta_smp = mcleod.obj$additional$original_stat_res$beta_smp
  ret$beta_suggestion = mcleod.obj$additional$original_stat_res$beta_suggestion
  ret$proposal_approved = mcleod.obj$additional$original_stat_res$proposal_approved
  return(ret)
}



#' Retrieve the estimated posterior mixing distribution, both density and CDF
#' 
#' The function computes the point-wise posterior mean/median of the mixing distribution, using the MCMC samples computed by \code{\link{mcleod}}. The user may also ask for the mixing distribution (CDF,density and \eqn{\vec{\pi}}) from a specific iteration of the Gibbs sampler.
#' 
#' @details Posterior mean/ median computed in a pointwise manner, across points of the support, in order to obtain an estimate of the mixing distribution (which is a function). Burn-in samples (specified by the function \code{\link{mcleod}}) are removed from the MCMC samples when computing estiamtes.
#' 
#' @param res object resulting from function \code{\link{mcleod}}
#' @param aggregate_by function to aggregate by, such as \code{\link{mean}} or \code{\link{median}} (used for computing the posterior mean / median, respectively, at each point of the support of the mixing distribution). May also be a user supplied function.
#' @param specific_iter Default value is \code{NULL}, meaning the function returns estimates obtained via posterior mean/ median across all gibbs samples. If the user picks a single number, this iteration of the gibbs sample is returned.
#'
#' @return A list with the following fields:
#'  \itemize{
#' \item{CDF}{ - if \code{specific_iter} is set to \code{NULL}, returns an \code{R} \code{\link{approxfun}} function, giving an estimate (pointwise posterior mean/median) of the CDF of the mixing distribution. If a specific iteration is requested, will return the CDF for mixing distribution, found in the \code{specific_iter} iteration of the MCMC sampler.}
#' \item{density}{ - if \code{specific_iter} is set to \code{NULL}, returns an \code{R} \code{\link{approxfun}} function, giving an estimate (pointwise posterior mean/median) of the density of the mixing distribution. If a specific iteration is requested, will return the density for mixing distribution, found in the \code{specific_iter} iteration of the MCMC sampler.}
#' \item{a.vec}{ - Parameter vector \eqn{\vec{a}}, used to define the support of the piece-wise constant mixing distribution. See full definition in the package vignette.}
#' \item{pi_smp}{ - Parameter vector \eqn{\vec{\pi}}, detailing the estiamted probability mass in each segment of \eqn{\vec{a}}. Computed via pointwise posterior mean/median across MCMC samples. See full definition in the package vignette. If a specific iteration is requested, will return the value of \eqn{\vec{\pi}} sampled in that specific iteration.}
#' }
#' @export
#'
#' @examples
#' see package vignette
mcleod.get.posterior.mixing.dist = function(res,aggregate_by = mean,specific_iter = NULL){
  if(class(res) != CLASS.NAME.MCLEOD){
    stop('input argument res for mcleod.get.posterior.mixing.dist must be an object returned from function mcleod')
  }
  
  if(specific_iter<1 | specific_iter> res$parameters_list$nr.gibbs){
    stop(paste0('specific_iter must be a valid MCMC iteration between 1 and ',res$parameters_list$nr.gibbs))
  }
  
  #get the MCMC samples distribution
  pi_smp = t(res$additional$original_stat_res$pi_smp)
  a.vec = res$parameters_list$a.vec
  
  #check if we need a specific iteration or the whole data (and then remove the burnin)
  if(!is.null(specific_iter)){
    pi_smp = pi_smp[specific_iter,]  
  }else{ #drop burnin
    pi_smp = pi_smp[(res$parameters_list$nr.gibbs.burnin+1):res$parameters_list$nr.gibbs,]
    pi_smp = apply(pi_smp,2,aggregate_by)
  }
  
  #estimate density and CDF
  CDF = approxfun(x = a.vec,
                  y = c(0,cumsum(pi_smp)),yleft = 0,yright = 1)
  
  density = approxfun(x = a.vec,y = c(pi_smp,pi_smp[length(pi_smp)]),
                   yleft = 0,yright = 0,method = 'constant')
  
  #return results
  ret = list()
  ret$CDF = CDF
  ret$density = density
  ret$a.vec = a.vec
  ret$pi_smp = pi_smp
}



#' An internal function used for running a random normal intercept regression, in order to initialize the slope coefficients in the MCMC chain.
#' 
#' Function is called from \code{\link{mcleod}}, when paramter \code{beta_init} is set to \code{NULL}. Will attempt to run \code{GLMER} model first, with Binomial/Poisson samples and a random normmal intercept. If fails, will use fixed Binomial/Poisson regression.
#' @param x Observed Counts/ successes for binomial samples
#' @param n Total number of draws for binomial sampeles. For Poisson regression, set this to \code{NULL}
#' @param covariates matrix of covatiates, row are samples, not intercept term needed
#' @param offset_p vector of offset values for the linear predictor. If none are specified by the user, pass a vector of zeroes.
#' @param family either 'binomial' or 'poisson'.
#'
#' @return a vector of length \code{ncol(covariates)}, giving the estimated slope coefficients, from a random normal intercept.
#' @keywords internal
#' @export
#'
#' @examples
init_mcleod_random_intercept_regression = function(x,n,covariates,offset_p,family = 'binomial'){
  #generate data frame for regression.
  model_dt = data.frame(c = x,nc = n-x) #successes and failures
  model_dt = cbind(model_dt,covariates) #add the covariates
  n_Z = ncol(covariates) #number of covariates
  Z_labels = paste0('Z',1:n_Z); #labels for all cols
  names(model_dt) = c('c','nc',  Z_labels);
  colnames(offset_p) = 'offset'
  model_dt$obs = 1:nrow(model_dt); #index for observations
  model_dt$offset_term = offset_p$offset
  
  # construct a formula for the fixed intercept model
  if(family == 'binomial'){
    formula_test_simple = paste0("cbind(c,nc)~ ")
  }else if(family == 'poisson'){
    formula_test_simple = paste0("c~ ");  
  }
  formula_test_simple = paste0(formula_test_simple,paste0('Z',1:n_Z,collapse = ' + '),' +offset(offset_term)');  
  model <- glm(formula = formula_test_simple,family=family,data=model_dt)  
  mcleod_init_basic = model$coefficients[-1]
  
  #under a try block, try a random intercept model
  try({
    library(lme4) #used for normal random intercept regression models
    
    # construct a formula for the model
    
    formula_test = paste0(formula_test_simple," + (1|obs)")

    #fit model
    glmer.model_random_intercept <- glmer(formula_test, 
                                          data=model_dt, family=family);
    # extract fitted slopes (intercept not needed)
    mcleod_init_basic = glmer.model_random_intercept@beta[-1] 
  })
  
  #return result
  return(mcleod_init_basic)
}



#' Internal function, used for performing input checks for input in \code{mcleod.posterior.estimates.random.effect} and \code{mcleod.predictive.CI}
#'
#' 
#' @param X See description in original function
#' @param N See description in original function
#' @param mcleod_res See description in original function
#' @param offset_vec See description in original function
#' @param is_Noise_Poisson See description in original function
#' @param covariates See description in original function
#' @param skip_checks_X_N Should checks be skipped for X and N. This is relecent in \code{mcleod.predictive.CI}, where the value of X is not known yet.
#'
#' @return posterior mean for slope cofficients if available. Else returns the value 0.
#' @keywords internal
#' @export
checks.input.posthoc.analysis = function(X, N, mcleod_res, offset_vec, is_Noise_Poisson,covariates,skip_checks_X_N){
  if(class(mcleod_res) != CLASS.NAME.MCLEOD){
    stop('input argument for mcleod_res must be an object returned from function mcleod')
  }
  
  nr.gibbs.burnin = mcleod_res$parameters_list$nr.gibbs.burnin
  if(!skip_checks_X_N){
    if(length(X)!= length(offset_vec)){
      stop('offset must be same length as X')
    }
    if(!is_Noise_Poisson){
      if(length(N)!=length(X))
        stop('for binomial errors, N must be same length as X')
    }else{
      if(!is.null(N))
        stop('for Poisson errors, N must be set to NULL')
    }  
  }
  if(!xor(mcleod_res$parameters_list$covariates_given,
          is.null(covariates))){
    stop('If mcleod_res trained on data with covariates, covariates must also be given for this data. If mcleod_res trained on data with no covariates, additional covariates cannot be introduced here.')
  }else{
    #check covariates
    if(mcleod_res$parameters_list$covariates_given & !skip_checks_X_N)
      if(nrow(covariates)!= length(X))
        stop('nrow of covariates must be size of data: each observation must have covaraites')
  }
  
  # get the posterior mean of slope coefficients, if available
  if(!is.null(covariates)){
    posterior_mean_vec = apply(mcleod_res$additional$original_stat_res$beta_smp[,-c(1:nr.gibbs.burnin),drop=F],1,mean)  
    if(length(posterior_mean_vec) != ncol(covariates)){
      stop('number of covariates used for training mcleod model, and for covariates argument must be the same')
    }
    if(!(is.null(offset_vec) & is.null(covariates)))
      if(length(offset_vec) != nrow(covariates)){
        stop( 'length of offset_vec and nrow(covariates) must be the same')
      }
  }else{
    posterior_mean_vec = 0
  }
  ret = list()
  ret$posterior_mean_vec = posterior_mean_vec
  return(ret)
}

#' Title
#'
#' @param X 
#' @param mcleod_res 
#' @param covariates 
#' @param method 
#' @param offset_vec 
#'
#' @return
#' @export
#'
#' @examples
mcleod.posterior.estimates.random.effect = function(X,N,mcleod_res,covariates = NULL,method = 'mean',offset_vec = rep(0,length(X))){
  
  if(!(method %in% c('mean','mode'))){
    stop(' method must be either "mean" or "mode" ')
  }
  nr.gibbs.burnin = mcleod_res$parameters_list$nr.gibbs.burnin
  a.vec = mcleod_res$parameters_list$a.vec
  pi_smp = (t(mcleod_res$additional$original_stat_res$pi_smp))[-(1:nr.gibbs.burnin),]
  pi_All = apply(pi_smp,2,mean)
  Noise_Type = mcleod_res$parameters_list$Noise_Type
  is_Noise_Poisson = (Noise_Type == MCLEOD.POISSON.ERRORS) #if false, it is binomial
  
  checks_result = checks.input.posthoc.analysis(X = X,
                                                N = N,
                                                mcleod_res = mcleod_res,
                                                offset_vec = offset_vec,
                                                is_Noise_Poisson = is_Noise_Poisson,
                                                covariates = covariates,
                                                skip_checks_X_N = F)
  
  posterior_mean_vec = checks_result$posterior_mean_vec
  n = length(x)
  Post.k.i = matrix(NA,nrow = n,ncol = length(a.vec)-1)
  
  compute_shift = function(k,posterior_mean_vec){
    if(!is.null(covariates)){
      sample_shift_in_log_odds_scale = sum(covariates[k,] * posterior_mean_vec + offset_vec[k])  
    }else{
      sample_shift_in_log_odds_scale = 0
    }
    return(sample_shift_in_log_odds_scale)
  }
  
  estimate_posterior_prob_vec = function(k,posterior_mean_vec,pi_All){
    sample_shift_in_log_odds_scale = compute_shift(k,posterior_mean_vec)
    if(!is_Noise_Poisson){
      integrand = function(u){
        dbinom(as.numeric(X[k]),size = N[k],
               prob = mcleod::inv.log.odds(u + sample_shift_in_log_odds_scale))
      }
    }else{
      integrand = function(u){
        dpois(as.numeric(X[k]),
              lambda = exp(u + sample_shift_in_log_odds_scale))
      }
    }
   
    Post.k = rep(NA,length(a.vec)-1)
    for(a_index in 1:(length(a.vec)-1)){
      #a_index = 1
      Post.k[a_index] = integrate(integrand,lower = a.vec[a_index],
                                  a.vec[a_index+1])$value*pi_All[a_index]
    }
    Post.k = Post.k/sum(Post.k)
    return(Post.k)
  }
  
  aux_compute_mean = function(k,Post.k.i_vec,posterior_mean_vec){
    sample_shift_in_log_odds_scale = compute_shift(k,posterior_mean_vec)
    ret = 0
    for(i in 1:length(Post.k.i_vec)){
      N_integral = 100
      seq_integral =  seq(from = a.vec[i],to = a.vec[i]+1,  length.out = N_integral)
      if(!is_Noise_Poisson){
        conditional_f_gamma_given_X = dbinom(x = X[k],
                                             size = N[k],
                                             prob = mcleod::inv.log.odds(theta = seq_integral +sample_shift_in_log_odds_scale))        
      }else{
        conditional_f_gamma_given_X = dpois(x = X[k],
                                             lambda = exp(theta = seq_integral +sample_shift_in_log_odds_scale))
      }

      if(sum(conditional_f_gamma_given_X) != 0){
        conditional_f_gamma_given_X = conditional_f_gamma_given_X/sum(conditional_f_gamma_given_X)
        ret  = ret  + Post.k.i_vec[i ] * sum(seq_integral * conditional_f_gamma_given_X)  
      }else{
        ret  = ret  + Post.k.i_vec[i ] * mean(seq_integral)
      }
      
    }
    return(ret)
  }
  
  aux_compute_mode = function(k,posterior_mean_vec){
    sample_shift_in_log_odds_scale = compute_shift(k,posterior_mean_vec)
    ret = 0
    prob = -Inf
    for(i in 1:length(pi_All)){
      u_values = seq(from = a.vec[i],to =a.vec[i+1], by = (a.vec[i+1] - a.vec[i])/10)
      for(u in u_values){
        if(!is_Noise_Poisson){
          temp_prob = dbinom(as.numeric(X[k]),size = N[k],
                             prob = mcleod::inv.log.odds(u + sample_shift_in_log_odds_scale))*pi_All[i]          
        }else{
          temp_prob = dpois(as.numeric(X[k]),
                             lambda = exp(u + sample_shift_in_log_odds_scale))*pi_All[i]          
        }
        
        if(temp_prob>=prob){
          prob = temp_prob
          ret = u
        }
      }
    }
    return(ret)
  }
  
  estimated_random_effects = rep(NA,n)
  for(k in 1:n){
    if(method == 'mean'){
      Post.k.i[k,] = estimate_posterior_prob_vec(k = k,posterior_mean_vec = posterior_mean_vec,pi_All = pi_All)
      estimated_random_effects[k] = aux_compute_mean(k,Post.k.i[k,],posterior_mean_vec) 
    }else{
      estimated_random_effects[k] = aux_compute_mode(k,posterior_mean_vec)
    }
  }
  
  return(estimated_random_effects)
  
}


mcleod.predictive.CI = function(N,mcleod_res,covariates = NULL,method = 'mean', offset_vec = NULL, CI.coverage = 0.95,Nr_Simulated_Values = 1000){
  
  N.gibbs = mcleod_res$parameters_list$nr.gibbs
  N.gibbs.burnin = mcleod_res$parameters_list$nr.gibbs.burnin
  a.vec = mcleod_res$parameters_list$a.vec
  
  
  Noise_Type = mcleod_res$parameters_list$Noise_Type
  is_Noise_Poisson = (Noise_Type == MCLEOD.POISSON.ERRORS) #if false, it is binomial
  
  if(is_Noise_Poisson){
    if(!is.null(covariates)){
      n = nrow(covariates)
    }else{
      n=1
    }
  }else{
    n =length(N)  
  }
  
  if(is.null(offset_vec)){
    offset_vec = rep(0,n)
  }
  
  checks_result = checks.input.posthoc.analysis(X = NULL,N = N,
                                                mcleod_res = mcleod_res,
                                                offset_vec = offset_vec,
                                                is_Noise_Poisson = is_Noise_Poisson,
                                                covariates = covariates,
                                                skip_checks_X_N = T)  
    
  Lower_mcleod = rep(NA,n)
  Upper_mcleod = rep(NA,n)  
  
  for(i_to_sample_for in 1:n){
    sample_gammas = rep(NA,Nr_Simulated_Values)
    beta_sampled_vec = NULL
    
    if(!is.null(covariates)){
      beta_sampled_vec = matrix(data = NA,nrow = Nr_Simulated_Values,ncol = ncol(covariates))
    }
    
    
    for(j in 1:Nr_Simulated_Values){
      ind_post_sample = sample((N.gibbs.burnin+1):N.gibbs,size = 1)
      if(!is.null(covariates))
        beta_sampled_vec[j,] = mcleod_res$additional$original_stat_res$beta_smp[,ind_post_sample]
      pi_vec =  mcleod_res$additional$original_stat_res$pi_smp[,ind_post_sample]
      sample_gammas[j] = sample(1:length(pi_vec),size = 1,prob = pi_vec)
    }
    
    sample_gammas = runif(n = Nr_Simulated_Values, min = a.vec[sample_gammas],max = a.vec[sample_gammas+1])
    
    Linear_Predictor = sample_gammas + offset_vec[i_to_sample_for]
    if(!is.null(covariates)){
      Linear_Predictor = Linear_Predictor + beta_sampled_vec %*% 
        t(covariates[i_to_sample_for,,drop=F])
    }
    
    if(!is_Noise_Poisson){
      sample_prob = mcleod::inv.log.odds( Linear_Predictor )  
    }else{
      sample_prob = exp( Linear_Predictor )
    }
    
    if(!is_Noise_Poisson){
      X_sampled_for_obs = rbinom(n = Nr_Simulated_Values, size = N[i_to_sample_for],
                                 prob = sample_prob)  
    }else{
      X_sampled_for_obs = rpois(n = Nr_Simulated_Values,lambda = sample_prob)
    }
    
    alpha.interval = (1-CI.coverage)
    Q_Lower = quantile(X_sampled_for_obs,probs = alpha.interval/2)
    Q_Upper = quantile(X_sampled_for_obs,probs = 1-alpha.interval/2)
    Lower_mcleod[i_to_sample_for] = Q_Lower
    Upper_mcleod[i_to_sample_for] = Q_Upper
    
  }
  ret = list()
  ret$Lower = Lower_mcleod
  ret$Upper = Upper_mcleod
  return(ret)
}


  
  
  
  