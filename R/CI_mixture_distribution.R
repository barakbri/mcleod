#mcleod package
#-  create functions for estimating rho, for single theha, then for all theths (interpolation)
#-  create function for solving the estimation problem.
#-  create new function for plotting

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
                                           rho.estimation.perm = 50,
                                           rho.q_for_calibration = c(0.2,0.4,0.6,0.8),
                                           do_serial = F,
                                           nr.cores = ceiling(detectCores()/2)){
  
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
  ret$do_serial = do_serial
  ret$rho.possible.values = rho.possible.values
  ret$rho.estimation.perm = rho.estimation.perm
  ret$rho.q_for_calibration = rho.q_for_calibration
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
    minimal_required_number_of_iterations = ceiling(alpha*nr.perm)
    
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


mcleod.CI.rho.calibration.constructor = function(
  bank,
  res_mcleod_holdout,
  CDF_holdout,
  alpha.one.sided,
  nr.perm = 200,
  possible_rhos = c(0.1,0.2,0.3,0.4,0.5),
  q_for_rho_optimization = c(0.2,0.4,0.6,0.8),
  verbose = F
){
  optimal_rho_by_theta_for_GE = rep(NA,length(q_for_rho_optimization))
  optimal_rho_by_theta_for_LE = rep(NA,length(q_for_rho_optimization))
  theta_points = rep(NA,length(q_for_rho_optimization))
  for(index_q in 1:length(q_for_rho_optimization)){
    #index_q = 2
    current_q = q_for_rho_optimization[index_q]
    a_ind_for_theta_for_current_q = which.min(abs(CDF_holdout - current_q))
    theta_current_q = res_mcleod_holdout$parameters_list$a.vec[a_ind_for_theta_for_current_q]
    theta_current_q_ind = which.min(abs(bank$CI_param$theta_vec - theta_current_q))
    theta_current_q = bank$CI_param$theta_vec[theta_current_q_ind]
    theta_points[index_q] = theta_current_q
    if(verbose)
      print(paste0('optimizing rho for q=',current_q,', which is equivalent to theta = ',theta_current_q,' in the holdout data'))
    
    q_lower_by_rho = rep(NA,length(possible_rhos))
    point_to_start_testing_GE = qbinom(p = alpha.one.sided,size = length(bank$N_vec),prob = current_q) / length(bank$N_vec)
    points_to_test_GE = which(bank$CI_param$q_vec<= point_to_start_testing_GE)
    points_to_test_GE = tail(points_to_test_GE,n = 3)
    
    if(length(points_to_test_GE)>=2){
      if(verbose)
        print(paste0('selecting rho for GE'))
      for(index_rho in 1:length(possible_rhos)){
        
        current_rho = possible_rhos[index_rho]
        
        a_ind_GE = which(res_mcleod_holdout$parameters_list$a.vec <= theta_current_q - current_rho)
        if(length(a_ind_GE) > 0 ){
          a_ind_GE = max(a_ind_GE)
        }else{
          a_ind_GE = 1
        }
        
        PVs_at_Qs = rep(NA,length(points_to_test_GE))
        for(i in 1:length(points_to_test_GE)){
          PVs_at_Qs[i] = mcleod.CI.PV.at_point(bank = bank,
                                               ind_theta = theta_current_q_ind,
                                               ind_q = points_to_test_GE[i],
                                               CDF_value = CDF_holdout[a_ind_GE],
                                               a_index = a_ind_GE,
                                               alpha = alpha.one.sided,
                                               nr.perm = nr.perm,
                                               do_check_vs_noiseless_case = F,
                                               do_check_vs_minimum_number_of_required_iterations = F,
                                               is_GE = T,
                                               do_serial = bank$CI_param$do_serial)
        }
        
        #create model:
        if(length(unique(PVs_at_Qs))>1){
          PVs_at_Qs = PVs_at_Qs*(nr.perm)/(nr.perm+1) #this is to make sure we dont have 1's. 0's are avoided by adding the test statistic to the perms when computing a PV value
          logits = log.odds(PVs_at_Qs)
          model = lm(logits~bank$CI_param$q_vec[points_to_test_GE])
          b0 = model$coefficients[1]
          b1 = model$coefficients[2]
          q_lower_by_rho[index_rho] = (log.odds(alpha.one.sided) - b0 )/b1
        }else{
          q_lower_by_rho[index_rho] = max(bank$CI_param$q_vec[points_to_test_GE])
        }
      }
      
      optimal_rho_by_theta_for_GE[index_q] =  which.max(q_lower_by_rho)
    }
    
    
    q_upper_by_rho = rep(NA,length(possible_rhos))
    point_to_start_testing_LE = qbinom(p = 1-alpha.one.sided,size = length(bank$N_vec),prob = current_q) / length(bank$N_vec)
    points_to_test_LE = which(bank$CI_param$q_vec>= point_to_start_testing_LE)
    points_to_test_LE = head(points_to_test_LE,n = 3)
    
    if(length(points_to_test_LE)>=2){
      if(verbose)
        print(paste0('selecting rho for LE'))
      for(index_rho in 1:length(possible_rhos)){
        
        current_rho = possible_rhos[index_rho]
        
        a_ind_LE = which(res_mcleod_holdout$parameters_list$a.vec >= theta_current_q + current_rho)
        if(length(a_ind_LE) > 0 ){
          a_ind_LE = min(a_ind_LE)
        }else{
          a_ind_LE = length(res_mcleod_holdout$parameters_list$a.vec)
        }
        
        PVs_at_Qs = rep(NA,length(points_to_test_LE))
        for(i in 1:length(points_to_test_LE)){
          PVs_at_Qs[i] = mcleod.CI.PV.at_point(bank = bank,
                                               ind_theta = theta_current_q_ind,
                                               ind_q = points_to_test_LE[i],
                                               CDF_value = CDF_holdout[a_ind_LE],
                                               a_index = a_ind_LE,
                                               alpha = alpha.one.sided,
                                               nr.perm = nr.perm,
                                               do_check_vs_noiseless_case = F,
                                               do_check_vs_minimum_number_of_required_iterations = F,
                                               is_GE = F,
                                               do_serial = bank$CI_param$do_serial)
        }
        
        #create model:
        PVs_at_Qs = PVs_at_Qs*(nr.perm)/(nr.perm+1) #this is to make sure we dont have 1's. 0's are avoided by adding the test statistic to the perms when computing a PV value
        if(length(unique(PVs_at_Qs))>1){
          logits = log.odds(PVs_at_Qs)
          model = lm(logits~bank$CI_param$q_vec[points_to_test_LE])
          b0 = model$coefficients[1]
          b1 = model$coefficients[2]
          q_upper_by_rho[index_rho] = (log.odds(alpha.one.sided) - b0 )/b1  
        }else{
          q_upper_by_rho[index_rho] = min(bank$CI_param$q_vec[points_to_test_LE])
        }
        
      }
      optimal_rho_by_theta_for_LE[index_q] = which.min(q_upper_by_rho)
    }
    
  }# end of loop over q
  
  for(i in 2:length(optimal_rho_by_theta_for_LE)){
    if(is.na(optimal_rho_by_theta_for_LE[i]))
      optimal_rho_by_theta_for_LE[i] = optimal_rho_by_theta_for_LE[i-1]
  }
  for(i in (length(optimal_rho_by_theta_for_GE)-1):1){
    if(is.na(optimal_rho_by_theta_for_GE[i]))
      optimal_rho_by_theta_for_GE[i] = optimal_rho_by_theta_for_GE[i+1]
  }
  optimal_rho_by_theta_for_LE = c(optimal_rho_by_theta_for_LE[1],
                                  optimal_rho_by_theta_for_LE,
                                  optimal_rho_by_theta_for_LE[length(optimal_rho_by_theta_for_LE)])
  optimal_rho_by_theta_for_GE = c(optimal_rho_by_theta_for_GE[1],
                                  optimal_rho_by_theta_for_GE,
                                  optimal_rho_by_theta_for_GE[length(optimal_rho_by_theta_for_GE)])
  
  optimal_rho_by_theta_for_LE = possible_rhos[optimal_rho_by_theta_for_LE]
  optimal_rho_by_theta_for_GE = possible_rhos[optimal_rho_by_theta_for_GE]
  theta_points = c(min(bank$CI_param$theta_vec),
                   theta_points,
                   max(bank$CI_param$theta_vec))
  
  rho_approx_fun_LE = approxfun(x = theta_points,y = optimal_rho_by_theta_for_LE,method = 'linear',rule = 2)
  rho_approx_fun_GE = approxfun(x = theta_points,y = optimal_rho_by_theta_for_GE,method = 'linear',rule = 2)
  ret = list(rho_approx_fun_LE = rho_approx_fun_LE,rho_approx_fun_GE = rho_approx_fun_GE)
  class(ret) = CLASS.NAME.MCLEOD.CI.RHO
  return(ret)
}



compute_P_values_over_grid_function=function(bank,rho_calibration_obj,res_mcleod_data,median_curve,alpha.one.sided,verbose = F){
  GE.pval.grid = matrix(NA,ncol = bank$CI_param$n_theta,nrow = bank$CI_param$n_q,
                        dimnames = list(paste0('q = ',bank$CI_param$q_vec),paste0('theta = ',bank$CI_param$theta_vec)))
  LE.pval.grid = matrix(NA,ncol = bank$CI_param$n_theta,nrow = bank$CI_param$n_q,
                        dimnames = list(paste0('q = ',bank$CI_param$q_vec),paste0('theta = ',bank$CI_param$theta_vec)))
  
  for(i in 1:bank$CI_param$n_q)
    for(j in 1:bank$CI_param$n_theta){
      if(verbose)
        print(paste0('testing GE at q_ind = ',i,'/',bank$CI_param$n_q,' , theta_ind = ',j,'/',bank$CI_param$n_theta))
      rho_GE = rho_calibration_obj$rho_approx_fun_GE(bank$CI_param$theta_vec[j])
      rho_LE = rho_calibration_obj$rho_approx_fun_LE(bank$CI_param$theta_vec[j])
      a_ind_GE = which(res_mcleod_data$parameters_list$a.vec <= bank$CI_param$theta_vec[j] - rho_GE)
      if(length(a_ind_GE) > 0 ){
        a_ind_GE = max(a_ind_GE)
      }else{
        a_ind_GE = 1
      }
      a_ind_LE = which(res_mcleod_data$parameters_list$a.vec >= bank$CI_param$theta_vec[j] + rho_LE)
      if(length(a_ind_LE) > 0 ){
        a_ind_LE = min(a_ind_LE)
      }else{
        a_ind_LE = length(res_mcleod_data$parameters_list$a.vec)
      }
      GE.pval.grid[i,j] = mcleod.CI.PV.at_point(bank = bank,
                                                ind_theta = j,
                                                ind_q = i,
                                                CDF_value = median_curve[a_ind_GE],
                                                a_index = a_ind_GE,
                                                alpha = alpha.one.sided,
                                                nr.perm = CI_param$nr.perms,
                                                do_check_vs_noiseless_case = T,
                                                do_check_vs_minimum_number_of_required_iterations = T,
                                                is_GE = T,
                                                do_serial = CI_param$do_serial)
      
      if(verbose)
        print(paste0('testing LE at q_ind = ',i,'/',bank$CI_param$n_q,
                     ' , theta_ind = ',j,'/',bank$CI_param$n_theta))
      LE.pval.grid[i,j] = mcleod.CI.PV.at_point(bank = bank,
                                                ind_theta = j,
                                                ind_q = i,
                                                CDF_value = median_curve[a_ind_LE],
                                                a_index = a_ind_LE,
                                                alpha = alpha.one.sided,
                                                nr.perm = CI_param$nr.perms,
                                                do_check_vs_noiseless_case = T,
                                                do_check_vs_minimum_number_of_required_iterations = T,
                                                is_GE = F,
                                                do_serial = CI_param$do_serial)
      
    }
  
  ret = list()
  ret$GE.pval.grid = GE.pval.grid
  ret$LE.pval.grid = LE.pval.grid
  return(ret)
}


compute_CI_curves_function = function(
  bank,
  res_mcleod_data ,
  rho_calibration_obj,
  median_curve,
  alpha.one.sided,
  verbose
){
  maximal_point_for_GE = rep(NA, bank$CI_param$n_theta)
  minimal_point_for_LE = rep(NA, bank$CI_param$n_theta)
  for(i in 1:bank$CI_param$n_theta){
    current_theta = bank$CI_param$theta_vec[i]
    ai_point_GE = min(which(res_mcleod_data$parameters_list$a.vec >= current_theta))
    ai_point_LE = max(which(res_mcleod_data$parameters_list$a.vec <= current_theta))
    maximal_point_for_GE[i] = max(which(bank$CI_param$q_vec<=median_curve[ai_point_GE]),1)
    minimal_point_for_LE[i] = min(which(bank$CI_param$q_vec>=median_curve[ai_point_LE]),bank$CI_param$n_q)
  }
  
  i_star_GE = rep(NA,bank$CI_param$n_theta)
  q_star_GE = rep(NA,bank$CI_param$n_theta)
  for(i in 1:bank$CI_param$n_theta){
    if(verbose){
      print(paste0('Computing GE confidence intervals, log.odds = ',bank$CI_param$theta_vec[i]))
    }
    #ai = 2
    if(i == 1){
      N_21_GE = 1
    }else{
      N_21_GE = min(i_star_GE[i-1] + 1,bank$CI_param$n_q)
    }
    N_22_GE = maximal_point_for_GE[i]
    
    
    rho_GE = rho_calibration_obj$rho_approx_fun_GE(bank$CI_param$theta_vec[i])
    a_ind_GE = which(res_mcleod_data$parameters_list$a.vec <= bank$CI_param$theta_vec[i] - rho_GE)
    if(length(a_ind_GE) > 0 ){
      a_ind_GE = max(a_ind_GE)
    }else{
      a_ind_GE = 1
    }
    
    ind_to_select_from = (max(N_21_GE,1)):N_22_GE
    rejected_ind = rep(NA,length(ind_to_select_from))
    for(u in 1:length(rejected_ind)){
      rejected_ind[u] = mcleod.CI.PV.at_point(bank = bank,
                                              ind_theta = i,
                                              ind_q = ind_to_select_from[u],
                                              CDF_value = median_curve[a_ind_GE],
                                              a_index = a_ind_GE,
                                              alpha = alpha.one.sided,
                                              nr.perm = CI_param$nr.perms,
                                              do_check_vs_noiseless_case = T,
                                              do_check_vs_minimum_number_of_required_iterations = T,
                                              is_GE = T,
                                              do_serial = CI_param$do_serial) <= alpha.one.sided
      if(is.na(rejected_ind[u])) #need to handle bounds and shifts better
        rejected_ind[u] = F  
      if(u==1 & rejected_ind[u] == F ){
        rejected_ind = rep(F,length(ind_to_select_from))
        break
      }
      
    }
    if(i == 1 & rejected_ind[1] == F){
      i_star_GE[i] = 0
    }else if (i == 1){
      i_star_GE[i] = ind_to_select_from[max(which(rejected_ind))]
    }
    
    if(i >1 & rejected_ind[1] == F){
      i_star_GE[i] = i_star_GE[i-1]
    }else if (i >1){
      i_star_GE[i] = ind_to_select_from[max(which(rejected_ind))]
    }
    if(i_star_GE[i] == 0){
      q_star_GE[i] = 0  
    }else{
      q_star_GE[i] = bank$CI_param$q_vec[i_star_GE[i]]   
    }
    
  }
  
  #maximum_point_to_consider_for_LE = max(which(is.finite(minimal_point_for_LE)))
  #maximum_point_to_consider_for_LE = min(maximum_point_to_consider_for_LE,length(a.vec)-1 )
  i_star_LE = rep(NA,bank$CI_param$n_theta)
  q_star_LE = rep(NA,bank$CI_param$n_theta)
  for(i in bank$CI_param$n_theta:1){
    if(verbose){
      print(paste0('Computing LE confidence intervals, log.odds = ',bank$CI_param$theta_vec[i]))
    }
    if(i == bank$CI_param$n_theta){
      N_21_LE = bank$CI_param$n_q
    }else{
      N_21_LE = max(i_star_LE[i+1] - 1,1)
    }
    N_22_LE = minimal_point_for_LE[i]
    
    rho_LE = rho_calibration_obj$rho_approx_fun_LE(bank$CI_param$theta_vec[i])
    a_ind_LE = which(res_mcleod_data$parameters_list$a.vec >= bank$CI_param$theta_vec[i] + rho_LE)
    if(length(a_ind_LE) > 0 ){
      a_ind_LE = min(a_ind_LE)
    }else{
      a_ind_LE = length(res_mcleod_data$parameters_list$a.vec)
    }
    
    ind_to_select_from = N_22_LE:N_21_LE # (max(N_21_LE,1)):N_22_LE
    
    rejected_ind = rep(NA,length(ind_to_select_from))
    
    for(u in length(ind_to_select_from):1){
      rejected_ind[u] = mcleod.CI.PV.at_point(bank = bank,
                                              ind_theta = i,
                                              ind_q = ind_to_select_from[u],
                                              CDF_value = median_curve[a_ind_LE],
                                              a_index = a_ind_LE,
                                              alpha = alpha.one.sided,
                                              nr.perm = CI_param$nr.perms,
                                              do_check_vs_noiseless_case = T,
                                              do_check_vs_minimum_number_of_required_iterations = T,
                                              is_GE = F,
                                              do_serial = CI_param$do_serial) <= alpha.one.sided
      
      if(is.na(rejected_ind[u])) #need to handle bounds and shifts better
        rejected_ind[u] = F
      if(u==length(ind_to_select_from) & rejected_ind[u] == F )
        break
    }
    
    if(i ==  bank$CI_param$n_theta & rejected_ind[length(rejected_ind)] == F){
      i_star_LE[i] = bank$CI_param$n_q+1
    }else if (i == bank$CI_param$n_theta){
      i_star_LE[i] = ind_to_select_from[min(which(rejected_ind))]
    }
    
    if(i < bank$CI_param$n_theta & rejected_ind[length(rejected_ind)] == F){
      i_star_LE[i] = i_star_LE[i+1]
    }else if (i < bank$CI_param$n_theta){
      i_star_LE[i] = ind_to_select_from[min(which(rejected_ind))]
    }
    
    if(i_star_LE[i] == bank$CI_param$n_q+1){
      q_star_LE[i] = 1  
    }else{
      q_star_LE[i] = bank$CI_param$q_vec[i_star_LE[i]]   
    }
    
  }
  ret = list()
  ret$q_star_LE = q_star_LE
  ret$q_star_GE = q_star_GE
  return(ret)
}
