
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
#' @details Details the type of prior for the random effect model. The prior for the distribution of random effects can be either Hierarchical Beta (a Polya Tree, i.e. a full binary tree with L levels, with each leaf corresponding to a segment on the real line, and each internal node associated with a beta random variable) or a two level dirichlet tree. See additional details in the package vignette.
#'
#' @param prior.type The type of prior used. Values are either \code{MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL} or \code{MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET}
#' @param Beta.Heirarchical.Levels The number of levels used (for prior of type MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL). The first level (root of the tree) is 1. A full tree with L=4 levels has 2^(4-1) segments for the prior. 
#' @param Two.Layer.Dirichlet.Intervals When using \code{MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET} as the prior - the total number of leafs in the two-layer dirichet tree.
#' @param Two.Layer.Dirichlet.Nodes.in.First.Layer Number of nodes in the first level of the Dirichlet tree. \code{Two.Layer.Dirichlet.Intervals} must be an integer multiple of \code{Two.Layer.Dirichlet.Nodes.in.First.Layer}
#'
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
#' @param integration_step_size Used for numerical integration when computing the probability for an observed value to originate from a given segment of the underlying density (specifically, from segments of a.vec)
#' @param Fast.Gamma.Used Should a "fast-and-inaccurate" sampling procedure from a gamma density be used.
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


#' An object of type 'mcleod.covariates.estimation.parameters.obj' defining computational and statistical parameters, when estimating a random effect model with covariates.
#' 
#' @details For binomial response, the model with covariates is log(p/(1-p)) = gamma + delta^T Z, where gamma is a random deviate from a general density estimated from the data, delta is a vector of slopes, and Z is a vector of covariates. For poisson errors, the mode equation is log(lambda) = gamma + delta^T Z
#' @param proposal_sd vector of standard deviations for the proposal distribution, for the values of delta. Proposals are generated by adding a random value from N(0,diag(proposal_sd)) to the current MCMC value of delta.
#' @param beta_prior_sd vector of standard deviations for the prior. The prior for delta is assumed to be N(0,diag(beta_prior_sd)). Setting an entry to a negative value (e.g. -1) will use an uninformative prior (constant) for that component of delta.
#' @param beta_init Initial values for beta. default value is zero, however non-cold-start values can be used, e.g. by using values from a binomial, or normal random-intercept binomial regression.
#' @param Manual_Prior_Values When using a predefined prior for delta, a vector of sorted values defining the paritioning points for a piece-wise constant prior for delta.
#' @param Manual_Prior_Probs A vector with length = length(Manual_Prior_Values)-1, giving the probablities for the different segments defined by  \code{Manual_Prior_Values}. Within each segment, it is assumed the prior is uniformly distributed.
#' @param do_P_k_i_hashing A computational workaround for faster computation of bayesian probabilites. See additional details in the appendix describing the computational method, in the dissertation of Brill (2022)
#' @param P_k_i_hashing_resolution_by_theta    Resolution parameter for the computational workaround activated by do_P_k_i_hashing
#'
#' @return an object of type 'mcleod.covariates.estimation.parameters.obj'
#' @export
#'
#' @examples See the package vignette
mcleod.covariates.estimation.parameters = function(proposal_sd = c(0.05),
                                           beta_prior_sd = c(1),
                                           beta_init = c(0),
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
      stop('Manual_Prior_Values must be longer by one entery than nr rows of Manual_Prior_Probs')      
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


#' The one-stop-shop function for density estimates for the mixing distribution, 
#' 
#' The function estimates the mixing distribution of P/lambda for binomial/Poisson samples, respectively.
#' 
#' @details Computes deconvolution estimates for the mixing distribution of the following model: For binomial measurements, X~bin(N,p) and log(p/(1-p)) is random deviate from either a Polya-Tree or Dirichlet Tree distribution, as defined by \code{prior_parameters}. For Poisson errors, the model equation is X~Pois(lambda), where log(lambda) is distributed using a Polya tree or Dirichlet tree distribution.
#' If covariates are included in the model, the model is log(p/(1-p)) = gamma + delta^T Z, where gamma is a random deviate from a general density estimated from the data, delta is a vector of slopes, and Z is a vector of covariates. For poisson errors, the mode equation is log(lambda) = gamma + delta^T Z.
#' 
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
                     offset_vec = rep(0,length(n.smp)),
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

#' Title
#'
#' @param mcleod.obj 
#' @param plot_only_point_estimate 
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param mcleod.obj 
#' @param plot.posterior 
#' @param plot.MH.proposal.by.iteration 
#'
#' @return
#' @export
#'
#' @examples
results.covariate.coefficients.posterior = function(mcleod.obj, plot.posterior = T, plot.MH.proposal.by.iteration = F){

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
  posterior_mean_vec = apply(mcleod.obj$additional$original_stat_res$beta_smp[,-c(1:burnin),drop=F],1,mean)
  
  # plot for each covariate, in a seperate frame, the distribution of suggested proposals. Accepted proposals are marked in red (note that a proposal is for all slope variables simultanuesly)
  if(plot.MH.proposal.by.iteration){
    par(mfrow=c(nr.covariates,1))  
    for(i in 1:nr.covariates){
      #plot, remove proposal for the last iteration
      plot(mcleod.obj$additional$original_stat_res$beta_suggestion[i, -mcleod.obj$parameters_list$nr.gibbs ],col = mcleod.obj$additional$original_stat_res$proposal_approved+1,pch= 20,main = paste0('beta',i,' - proposal'),xlab = 'Iteration',ylab = 'Coefficient')    
    }
    par(mfrow=c(1,1))
  }
  
  #return the posterior means and the acceptance rate
  ret = list()
  ret$posterior.means = posterior_mean_vec
  ret$acceptance.rate = mean(mcleod.obj$additional$original_stat_res$proposal_approved)
  return(ret)
}