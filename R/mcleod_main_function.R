
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#load libraries required:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Constants and definitions:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CLASS.NAME.MCLEOD = 'mcleod.obj'
CLASS.NAME.PRIOR.DEFINITION = 'mcleod.prior.def.obj'
CLASS.NAME.COMPUTATIONAL.PARAMETERS.DEFINITION = 'mcleod.computational.parameters.obj'
CLASS.NAME.COVARIATES.ESTIMATION.PARAMETERS.DEFINITION = 'mcleod.covariates.estimation.parameters.obj'

MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL = 1L
MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET = 0L

MCLEOD.BINOMIAL.ERRORS = 0L
MCLEOD.POISSON.ERRORS = 1L


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Constructors for auxiliary objects:
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Title
#'
#' @param prior.type 
#' @param Beta.Heirarchical.Levels 
#' @param Two.Layer.Dirichlet.Intervals 
#' @param Two.Layer.Dirichlet.Nodes.in.First.Layer
#'
#' @return
#' @export
#'
#' @examples
mcleod.prior.parameters = function(prior.type = MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL,
                                   Beta.Heirarchical.Levels = 6,
                                   Two.Layer.Dirichlet.Intervals = 64,
                                   Two.Layer.Dirichlet.Nodes.in.First.Layer = 8){
  
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
  
  ret = list()
  ret$prior.type = prior.type
  ret$Beta.Heirarchical.Levels = Beta.Heirarchical.Levels
  ret$Two.Layer.Dirichlet.Intervals = Two.Layer.Dirichlet.Intervals
  ret$Two.Layer.Dirichlet.Nodes.in.First.Layer = Two.Layer.Dirichlet.Nodes.in.First.Layer
  class(ret) = CLASS.NAME.PRIOR.DEFINITION
  return(ret)
}

#' Title
#'
#' @param nr.gibbs 
#' @param nr.gibbs.burnin 
#' @param integration_step_size 
#' @param Fast.Gamma.Used 
#' @param Fast.Gamma.Bank.Size 
#'
#' @return
#' @export
#'
#' @examples
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
    warning('warning: nr.gibbs.burnin lower than 100 iterations')
  if(nr.gibbs < 500)
    warning('warning: nr.gibbs lower than 100 iterations')
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


#' Title
#'
#' @param proposal_sd 
#' @param beta_prior_sd 
#' @param beta_init 
#' @param Manual_Prior_Values 
#' @param Manual_Prior_Probs  
#'
#' @return
#' @export
#'
#' @examples
mcleod.covariates.estimation.parameters = function(proposal_sd = c(0.05),
                                           beta_prior_sd = c(5),
                                           beta_init = c(0),
                                           Manual_Prior_Values = NULL,
                                           Manual_Prior_Probs = NULL){
  
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
  class(ret) = CLASS.NAME.COVARIATES.ESTIMATION.PARAMETERS.DEFINITION
  return(ret)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# main function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Title
#'
#' @param x.smp 
#' @param n.smp 
#' @param a.limits 
#' @param Noise_Type 
#' @param covariates 
#' @param prior_parameters 
#' @param computational_parameters 
#' @param covariates_estimation_parameters 
#' @param input_P_k_i 
#'
#' @return
#' @export
#'
#' @examples
mcleod	<- function( x.smp,
                     n.smp,
                     a.limits = c(-4,4),
                     Noise_Type = MCLEOD.BINOMIAL.ERRORS,
                     covariates = NULL,
                     prior_parameters = NULL,
                     computational_parameters = NULL,
                     covariates_estimation_parameters = NULL,
                     input_P_k_i = NULL,
                     exact.numeric.integration = TRUE
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
  
  #%%% Function settings
  
  
  if(is.null(input_P_k_i)){
    P_k_i_is_given = 0L
    P_k_i_precomputed = matrix(1,nrow = 1)
    K			<- length(x.smp)
  }else{
    P_k_i_is_given = 1L
    P_k_i_precomputed = input_P_k_i
    K			<- nrow(P_k_i_precomputed)
  }
  
  #checks
  #Noise type legal
  if(!(Noise_Type %in% c(MCLEOD.BINOMIAL.ERRORS,MCLEOD.POISSON.ERRORS))){
    stop('Noise_Type must be either set to MCLEOD.BINOMIAL.ERRORS or MCLEOD.POISSON.ERRORS')
  }
  if(!P_k_i_is_given & Noise_Type == MCLEOD.BINOMIAL.ERRORS){
   if(length(x.smp) != length(n.smp)){
     stop('For binomial errors, x.smp and n.smp must be of equal length')
   } 
  }
  if(!P_k_i_is_given & Noise_Type == MCLEOD.POISSON.ERRORS){
    if(!is.null(n.smp)){
      stop('For poisson errors, n.smp must be set to NULL')
    }
    n.smp = x.smp
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
  
  I = 2^(L)
  if(Prior_Type == MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET){
    if(!is.null(I_specificy_parameter))
        I = I_specificy_parameter
    if(I/I1 != as.integer(I/I1)){
     stop('I/I1 not an integer!') 
    }
  }
  
  a.vec.used		<- seq(a.min,a.max,length = I+1)  
  
  if(P_k_i_is_given){
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
    # if(P_k_i_is_given & !(is.null(x.smp) & is.null(n.smp))){
    # stop('if P_k_i is given, x.smp and n.smp must be set to NULL')  
    # }
    # x.smp = 1 # replace to numeric value - instead of null type
    # n.smp = 1
  }
  
  #handle manual priors:
  Manual_Prior_Given = c(0L)
  Manual_Prior_Values = c(-4,4)
  Manual_Prior_Probs = c(1)
  if(!is.null(covariates_estimation_parameters$Manual_Prior_Values)){
    Manual_Prior_Given = c(1L)
    Manual_Prior_Values = covariates_estimation_parameters$Manual_Prior_Values
    Manual_Prior_Probs = covariates_estimation_parameters$Manual_Prior_Probs
    if(ncol(covariates) != dim(Manual_Prior_Probs)[2]){
      stop('ncols of  Manual_Prior_Probs must be the same as the number of covariates')
    }
    
  }
    
    
  #%%% call Rcpp wrapper
  res = Wrapper_rcpp_Gibbs(x.smp,
                           n.smp,
                           a.vec.used,
                           nr.gibbs,nr.gibbs.burnin,
                           as.integer(exact.numeric.integration),
                           as.integer(0), #verbose - turned off
                           L,
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
                           Manual_Prior_Probs = Manual_Prior_Probs
                           )
  
  #%%% Wrap results
  
  ret = list()
  class(ret) = CLASS.NAME.MCLEOD
  
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
  
  ret$additional = list(original_stat_res = res)
  return(ret)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions for plotting results
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Title
#'
#' @param mcleod.obj 
#'
#' @return
#' @export
#'
#' @examples
plot.posterior	<- function(mcleod.obj)
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
  gg_obj = ggplot(dt_gibbs_cloud)  + ylim(c(0,1)) +
    geom_line(aes(x = a.vec,y = means_vec),colour = 'red',data = dt_mean_line) +
    geom_point(aes(x = a.point,y = CDF.value),alpha = 0.25,colour = 'gray',data = dt_gibbs_cloud,size = 0.8, shape = 18)+ xlab('theta')+ylab('CDF')
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
  covariates = mcleod.obj$parameters_list$covariates
  burnin = mcleod.obj$parameters_list$nr.gibbs.burnin
  nr.covariates = ncol(covariates)
  mean_vec = rep(NA,nr.covariates)
  
  if(plot.posterior){
    col_vec = rep(1,mcleod.obj$parameters_list$nr.gibbs)
    col_vec[1:burnin] = 2
    par(mfrow=c(nr.covariates,1))  
    for(i in 1:nr.covariates){
      plot(mcleod.obj$additional$original_stat_res$beta_smp[i,],col = col_vec,pch= 20,main = paste0('Coefficient for covariate ',i),xlab = 'Iteration',ylab = 'Coefficient')
    }  
    par(mfrow=c(1,1))
  }
  
  posterior_mean_vec = apply(mcleod.obj$additional$original_stat_res$beta_smp[,-c(1:burnin),drop=F],1,mean)
  
  if(plot.MH.proposal.by.iteration){
    par(mfrow=c(nr.covariates,1))  
    for(i in 1:nr.covariates){
      #plot, remove proposal for the last iteration
      plot(mcleod.obj$additional$original_stat_res$beta_suggestion[i, -mcleod.obj$parameters_list$nr.gibbs ],col = mcleod.obj$additional$original_stat_res$proposal_approved+1,pch= 20,main = paste0('beta',i,' - proposal'),xlab = 'Iteration',ylab = 'Coefficient')    
    }
    par(mfrow=c(1,1))
  }
  
  ret = list()
  ret$posterior.means = posterior_mean_vec
  ret$acceptance.rate = mean(mcleod.obj$additional$original_stat_res$proposal_approved)
  return(ret)
}