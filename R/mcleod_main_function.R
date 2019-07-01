
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
mcleod.prior.parameters = function(prior.type = MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL,
                                   Beta.Heirarchical.Levels = 6,
                                   Two.Layer.Dirichlet.Intervals = 64,
                                   Two.Layer.Dirichlet.Nodes.in.First.Layer = 8){
  
  #checks:
  #Prior type legal
  #intervals legal
  #complete multiple
  
  ret = list()
  ret$prior.type = prior.type
  ret$Beta.Heirarchical.Levels = Beta.Heirarchical.Levels
  ret$Two.Layer.Dirichlet.Intervals = Two.Layer.Dirichlet.Intervals
  ret$Two.Layer.Dirichlet.Nodes.in.First.Layer = Two.Layer.Dirichlet.Nodes.in.First.Layer
  class(ret) = CLASS.NAME.PRIOR.DEFINITION
  return(ret)
}

mcleod.computational.parameters = function(nr.gibbs = 500,
                                                   nr.gibbs.burnin = 250,
                                                   integration_step_size = 0.01,
                                                   Fast.Gamma.Used = F,
                                                   Fast.Gamma.Bank.Size = 1000L){
  #checks:
  #nr.gibbs > nr.gibbs burnin
  #warnings nr.gibbs and burnins
  #warning on integration_step_size
  
  ret = list()
  ret$nr.gibbs = nr.gibbs
  ret$nr.gibbs.burnin = nr.gibbs.burnin
  ret$integration_step_size = integration_step_size
  ret$Fast.Gamma.Used = Fast.Gamma.Used
  ret$Fast.Gamma.Bank.Size = Fast.Gamma.Bank.Size
  class(ret) = CLASS.NAME.COVARIATES.ESTIMATION.PARAMETERS.DEFINITION
  return(ret)
}


mcleod.covariates.estimation.parameters = function(proposal_sd = c(0.05),
                                           beta_prior_sd = c(5),
                                           beta_init = c(0)){
  
  #checks:
  
  ret = list()
  ret$proposal_sd = proposal_sd
  ret$beta_prior_sd = beta_prior_sd
  ret$beta_init = beta_init
  class(ret) = CLASS.NAME.COMPUTATIONAL.PARAMETERS.DEFINITION
  return(ret)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# main function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mcleod	<- function( x.smp,
                     n.smp,
                     a.limits = c(-4,4),
                     Noise_Type = MCLEOD.BINOMIAL.ERRORS,
                     covariates = NULL,
                     prior_parameters = NULL,
                     computational_parameters = NULL,
                     covariates_estimation_parameters = NULL
                     )
{
  
  exact.numeric.integration = TRUE # We force exact numeric integration over dbinom for computation of P_k_i. Normal approximation is not sufficiant.
  
  
  #%%% Retreive prior parameters
  
  
  if(is.null(prior_parameters)){
    prior_parameters = mcleod.prior.parameters()
  }
  #check object type
  L = prior_parameters$Beta.Heirarchical.Levels
  I1 = prior_parameters$Two.Layer.Dirichlet.Nodes.in.First.Layer
  Prior_Type = prior_parameters$prior.type
  I_specificy_parameter = prior_parameters$Two.Layer.Dirichlet.Intervals
  
  #%%% Retreive computational parameters
  if(is.null(computational_parameters)){
    computational_parameters = mcleod.computational.parameters()
  }
  #check object type
  nr.gibbs = computational_parameters$nr.gibbs
  nr.gibbs.burnin = computational_parameters$nr.gibbs.burnin
  integration_step_size = computational_parameters$integration_step_size
  Fast.Gamma.Used = computational_parameters$Fast.Gamma.Used
  Fast.Gamma.Bank.Size = computational_parameters$Fast.Gamma.Bank.Size
  
  #%%% Retreive covariate estimation parameters
  if(is.null(covariates_estimation_parameters)){
    covariates_estimation_parameters = mcleod.covariates.estimation.parameters()
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
  
  
  #checks
  #x.smp, n.smp lengths match
  #a.limits in the correct order
  #Noise type legal
  #for poisson - check if a.limits are in appropriate
  
  
  a.max = a.limits[2]
  a.min = a.limits[1]
  K			<- length(x.smp)
  
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
  
  #check precomputed P_k_i is same size as a.vec - this can be done in Wrapper_rcpp_Gibbs
  #check that if Precomputed P_k_i - there are no covariates given!
    
  
  
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
                           Noise_Type = Noise_Type
                           )
  
  #%%% Wrap results
  
  ret = list()
  class(ret) = CLASS.NAME.MCLEOD
  
  ret$parameters_list = list(
    a.vec = a.vec.used,
    nr.gibbs = nr.gibbs,
    nr.gibbs.burnin = nr.gibbs.burnin,
    L = L,
    Fast.Gamma.Used = Fast.Gamma.Used,
    Fast.Gamma.Bank.Size = Fast.Gamma.Bank.Size
  )
  
  ret$additional = list(original_stat_res = res)
  return(ret)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions for plotting results
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot.posterior	<- function(mcleod.obj)
{
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


