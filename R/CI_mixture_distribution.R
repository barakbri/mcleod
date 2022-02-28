#mcleod package
#- (Tuesday) create functions for estimating rho, for single theha
#- (Wednesday) create function for solving the estimation problem.
#- (Wednesday) create new function for plotting

# Add functions as imports
# Add checks to this file
# Finish package documentation
# Build a vignette

library(hash)
library(doRNG)
library(doParallel)
library(parallel)


CLASS.NAME.MCLEOD.CI = 'mcleod.CI.obj'
CLASS.NAME.MCLEOD.CI.RHO = 'mcleod.CI.obj.rho'
CLASS.NAME.MCLEOD.CI.PARAMETERS = 'mcleod.CI.obj.parameters'
CLASS.NAME.MCLEOD.CI.DECONV.BANK = 'mcleod.CI.obj.deconv.bank'

mcleod.CI.estimation.parameters = function(q_vec = seq(0.1,0.9,0.1),
                                           theta_vec = seq(-3,3,0.25),
                                           a.limits = c(-5,5),
                                           sampling_distribution = 'binomial',
                                           comp_parameters = mcleod.computational.parameters(),
                                           prior_parameters = mcleod.prior.parameters(),
                                           nr.perms = 200,
                                           alpha.CI = 0.95,
                                           rho.estimation.method = NULL,
                                           rho.possible.values = seq(0.1,0.5,0.1),
                                           nr.cores = detectCores() - 1){
  
  library(doRNG)
  library(doParallel)
  library(parallel)
  
  # need to add checks
  ret = list()
  ret$q_vec = q_vec
  ret$n_q = length(q_vec)
  ret$theta_vec = theta_vec
  ret$n_theta = length(theta_vec)
  ret$a.limits = a.limits
  ret$sampling_distribution = sampling_distribution
  ret$comp_parameters = comp_parameters
  ret$prior_parameters = prior_parameters
  ret$rho.estimation.method = rho.estimation.method
  ret$nr.cores = nr.cores
  ret$nr.perms = nr.perms
  ret$alpha.CI = alpha.CI
  class(ret) = CLASS.NAME.MCLEOD.CI.PARAMETERS
  return(ret)
}


mcleod.CI.deconv.bank.constructor = function(N_vec=NULL,CI_param){
  n_theta = CI_param$n_theta
  n_q = CI_param$n_q
  
  median_curve_GE = list()
  median_curve_LE = list()
  
  for(i in 1:n_theta){
    median_curve_GE[[i]] = list()
    median_curve_LE[[i]] = list()
    for(j in 1:n_q){
      median_curve_GE[[i]][[j]] = list()
      median_curve_LE[[i]][[j]] = list()
    }
  }
  ret = list()
  ret$CI_param = CI_param
  ret$median_curve_GE = median_curve_GE
  ret$median_curve_LE = median_curve_LE
  ret$N_vec = N_vec
  class(ret) = CLASS.NAME.MCLEOD.CI.DECONV.BANK
  return(ret)
}


compute_medians_curve = function(mcleod_for_data){
  pi_smp_for_data = (t(mcleod_for_data$additional$original_stat_res$pi_smp))
  pi_smp_for_data = pi_smp_for_data[-(1:mcleod_for_data$parameters_list$nr.gibbs.burnin),]
  cumulative_pi_smp_for_data = t(apply(pi_smp_for_data,1,cumsum))
  median_cumulative_pi_smp_for_data = c(0,apply(cumulative_pi_smp_for_data,2,median))
  return(median_cumulative_pi_smp_for_data)
}


mcleod.CI.deconv.bank.get_median_curves_for_worst_case_hypothesis=function(bank,
                                                                      ind_theta,
                                                                      ind_q,
                                                                      nr.curves,
                                                                      is_GE = T,
                                                                      do_serial = T){
  if(bank$CI_param$sampling_distribution != 'binomial')
    stop('bank only supports binomial')
  
  if(is_GE)
    nr.computed = length(bank$median_curve_GE[[ind_theta]][[ind_q]])
  else
    nr.computed = length(bank$median_curve_LE[[ind_theta]][[ind_q]])
  
  nr.to.compute = nr.curves - nr.computed 
  
  if(nr.to.compute>0){
    
    if(bank$CI_param$sampling_distribution == 'binomial'){
      
      worker_function_Binomial_GE = function(seed){
        set.seed(seed)
        
        current_theta = bank$CI_param$theta_vec[ind_theta]
        current_q = bank$CI_param$q_vec[ind_q]
        N_vec = bank$N_vec
        P_sample = rbinom(n = length(N_vec),size = 1,1-current_q) * inv.log.odds(current_theta)
        X_sampled = rbinom(n = length(N_vec),size = N_vec,prob = P_sample)
        
        library(mcleod)
        
        temp_mcleod = mcleod(x.smp = X_sampled,n.smp = N_vec,
                             a.limits = bank$CI_param$a.limits,
                             exact.numeric.integration = T,
                             computational_parameters = bank$CI_param$comp_parameters,
                             prior_parameters = bank$CI_param$prior_parameters)
        
        return(compute_medians_curve(temp_mcleod))
      }
      if(do_serial){
        res_Binomial_GE = list()
        for(k in 1:nr.to.compute){
          res_Binomial_GE[[k]] = worker_function_Binomial_GE(k)
        }
      }else{
        #this assumes a cluster is registered
        res_Binomial_GE = foreach(seed=(nr.computed+1):nr.curves, .options.RNG=1,
                                  .export = c('ind_theta','ind_q','bank','worker_function_Binomial_GE')) %dorng% {
                                    worker_function_Binomial_GE(seed)
                                    
                                  }
      }
    
      
      for(k in 1:length(res_Binomial_GE)){
        bank$median_curve_GE[[ind_theta]][[ind_q]][[k + nr.computed]] <<- res_Binomial_GE[[k]]
        bank$median_curve_LE[[ind_theta]][[ind_q]][[k + nr.computed]] <<- rev(1- res_Binomial_GE[[k]])
      }
    }
  }
  
  if(is_GE)
    return(bank$median_curve_GE[[ind_theta]][[ind_q]])
  return(bank$median_curve_LE[[ind_theta]][[ind_q]])
  
}



mcleod.CI.deconv.bank.get_median_curves_for_worst_case_hypothesis_at_point = function(bank,
                                                                                      ind_theta,
                                                                                      ind_q,
                                                                                      a_index,
                                                                                      nr.perms,
                                                                                      is_GE = T,
                                                                                      do_serial = F){
  mcleod.CI.deconv.bank.get_median_curves_for_worst_case_hypothesis(bank = bank,
                                                                    ind_theta = ind_theta,
                                                                    ind_q = ind_q,
                                                                    nr.curves = nr.perms,
                                                                    is_GE = is_GE,
                                                                    do_serial = do_serial)
  
  if(is_GE)
    computed_curves_object = bank$median_curve_GE
  else
    computed_curves_object = bank$median_curve_LE
  
  ret = unlist(lapply(computed_curves_object[[ind_theta]][[ind_q]][1:nr.perms],FUN = function(x){x[a_index]}))
  return(ret)    
}