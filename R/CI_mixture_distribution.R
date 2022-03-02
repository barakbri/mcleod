#mcleod package
#- (Tuesday) create functions for estimating rho, for single theha, then for all theths (interpolation)
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
        
        if(is_GE){
          current_theta = bank$CI_param$theta_vec[ind_theta]
          current_q = bank$CI_param$q_vec[ind_q]  
        }else{ # we compute the corresponding GE hypothesis and revert the curves afterwords
          current_theta = bank$CI_param$theta_vec[bank$CI_param$n_theta - ind_theta + 1]
          current_q = bank$CI_param$q_vec[bank$CI_param$n_q - ind_q +1]  
        }
        
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
        res_Binomial = list()
        for(k in 1:nr.to.compute){
          res_Binomial[[k]] = worker_function_Binomial_GE(k)
        }
      }else{
        #this assumes a cluster is registered
        res_Binomial = foreach(seed=(nr.computed+1):nr.curves, .options.RNG=1,
                                  .export = c('ind_theta','ind_q','bank','worker_function_Binomial_GE')) %dorng% {
                                    worker_function_Binomial_GE(seed)
                                    
                                  }
      }
      
      #Collecting results
      ret = list()
      if(is_GE){
        if(nr.computed>0){
          for(k in 1:nr.computed){
            ret[[k]] = bank$median_curve_GE[[ind_theta]][[ind_q]][[k]]
          }
        }
        for(k in 1:length(res_Binomial)){
          ret[[k + nr.computed]] = res_Binomial[[k]]
          bank$median_curve_GE[[ind_theta]][[ind_q]][[k + nr.computed]] <<- res_Binomial[[k]]
          bank$median_curve_LE[[bank$CI_param$n_theta - ind_theta + 1]][[bank$CI_param$n_q - ind_q +1]][[k + nr.computed]] <<- rev(1- res_Binomial[[k]])
        }  
      }else{
        if(nr.computed>0){
          for(k in 1:nr.computed){
            ret[[k]] = bank$median_curve_LE[[ind_theta]][[ind_q]][[k]]
          }
        }
        for(k in 1:length(res_Binomial)){
          ret[[k + nr.computed]] = rev(1- res_Binomial[[k]])
          bank$median_curve_GE[[bank$CI_param$n_theta - ind_theta + 1]][[bank$CI_param$n_q - ind_q +1]][[k + nr.computed]] <<- res_Binomial[[k]]
          bank$median_curve_LE[[ind_theta]][[ind_q]][[k + nr.computed]] <<- rev(1- res_Binomial[[k]])
        }
      }
      
    }# end of case binomial
  }else{ # if there are none to compute, we can just take from bank
    ret = list()
    if(is_GE){
      for(k in 1:nr.curves){
        ret[[k]] = bank$median_curve_GE[[ind_theta]][[ind_q]][[k]]
      }
    }else{
      for(k in 1:nr.curves){
        ret[[k]] = bank$median_curve_LE[[ind_theta]][[ind_q]][[k]]
      }
    }
  }
  
  
  return(ret)
  
}



mcleod.CI.deconv.bank.get_median_curves_for_worst_case_hypothesis_at_point = function(bank,
                                                                                      ind_theta,
                                                                                      ind_q,
                                                                                      a_index,
                                                                                      nr.perms,
                                                                                      is_GE = T,
                                                                                      do_serial = F){
  computed_curves_object = mcleod.CI.deconv.bank.get_median_curves_for_worst_case_hypothesis(bank = bank,
                                                                    ind_theta = ind_theta,
                                                                    ind_q = ind_q,
                                                                    nr.curves = nr.perms,
                                                                    is_GE = is_GE,
                                                                    do_serial = do_serial)
  
  ret = unlist(lapply(computed_curves_object,FUN = function(x){x[a_index]}))
  return(ret)    
}




mcleod.CI.lower_bound_PV_for_worst_case=function(bank,
                                                 ind_q,
                                                 CDF_value,
                                                 is_GE = F){
  
  current_q = bank$CI_param$q_vec[ind_q]
  n = length(bank$N_vec)
  t_value = n*CDF_value
  prob_at_t_value = 0
  if(is_GE & t_value == round(t_value)){
    prob_at_t_value = dbinom(x = t_value,size = n,prob = current_q)
  }
  if(!is_GE & t_value == round(t_value)){
    prob_at_t_value = dbinom(x = n-t_value,size = n,prob = 1-current_q)
  }
  
  if(is_GE)
    return(pbinom(q = t_value,size = n,prob = current_q,lower.tail = F)+prob_at_t_value)
  return(pbinom(q = n-t_value,size = n,prob = 1-current_q,lower.tail = F)+prob_at_t_value)
}


mcleod.CI.PV.at_point = function(bank,
                                 ind_theta,
                                 ind_q,
                                 CDF_value,
                                 a_index,
                                 alpha = 0.05,
                                 nr.perm = 200,
                                 do_check_vs_noiseless_case = T,
                                 do_check_vs_minimum_number_of_required_iterations = T,
                                 is_GE = F,
                                 do_serial = F){
  
  current_q = bank$CI_param$q_vec[ind_q]
  current_theta = bank$CI_param$theta_vec[ind_theta]
  n = length(bank$N_vec)
  
  ret = NA
  if(do_check_vs_noiseless_case){
    ret = mcleod.CI.lower_bound_PV_for_worst_case(bank = bank,ind_q = ind_q,CDF_value = CDF_value,is_GE = is_GE)
    if(ret>=alpha){
      print(paste0('Break early by worst case'))
      return(1)
    }
  }
  
  if(do_check_vs_minimum_number_of_required_iterations){
    minimal_required_number_of_iterations = alpha*nr.perm
    
    null_stat_values_for_quick_test = mcleod.CI.deconv.bank.get_median_curves_for_worst_case_hypothesis_at_point(bank,
                                                                                                                 ind_theta = ind_theta,
                                                                                                                 ind_q = ind_q,
                                                                                                                 a_index = a_index,
                                                                                                                 nr.perms = minimal_required_number_of_iterations,
                                                                                                                 is_GE = is_GE,
                                                                                                                 do_serial = do_serial)
    if(( is_GE & all(null_stat_values_for_quick_test >= CDF_value)) |
       (!is_GE & all(null_stat_values_for_quick_test <= CDF_value))){
      print(paste0('Break early by iterations'))
      return(1)
    }
    
  }
  
  null_stat_values = mcleod.CI.deconv.bank.get_median_curves_for_worst_case_hypothesis_at_point(bank,
                                                                                                ind_theta = ind_theta,
                                                                                                ind_q = ind_q,
                                                                                                a_index = a_index,
                                                                                                nr.perms = nr.perm,
                                                                                                is_GE = is_GE,
                                                                                                do_serial = do_serial)
  
  if(is_GE)
    return(mean(c(null_stat_values,CDF_value) >= CDF_value))
  return(mean(c(null_stat_values,CDF_value) <= CDF_value))
  
}

