if(F){
  
  CI_param = mcleod.CI.estimation.parameters()
  CI_param$nr.perms
  
  bank = mcleod.CI.deconv.bank.constructor(rep(20,100),CI_param)
  
  cl <- makeCluster(8)
  registerDoParallel(cl)
  
  
  temp = mcleod.CI.deconv.bank.get_median_curves_for_worst_case_hypothesis(bank = bank,
                                                                      ind_theta = 13,
                                                                      ind_q = 5,
                                                                      nr.curves = 1000,
                                                                      is_GE = T,
                                                                      do_serial = F)
  length(bank$median_curve_LE[[13]][[5]])
  plot(bank$median_curve_GE[[13]][[5]][[1]])
  plot(bank$median_curve_LE[[13]][[5]][[1]])
  
  
  system.time({temp = mcleod.CI.deconv.bank.get_median_curves_for_worst_case_hypothesis(bank = bank,
                                                                                   ind_theta = 13,
                                                                                   ind_q = 5,
                                                                                   nr.curves = 1000,
                                                                                   is_GE = T,
                                                                                   do_serial = F)})
  
  mcleod.CI.deconv.bank.get_median_curves_for_worst_case_hypothesis_at_point(bank, ind_theta = 13,
                                                                             ind_q = 5,a_index = 5,
                                                                             nr.perms = 100,
                                                                             is_GE = T,
                                                                             do_serial = F)
  
  mcleod.CI.deconv.bank.get_median_curves_for_worst_case_hypothesis_at_point(bank, ind_theta = 13,
                                                                             ind_q = 5,a_index = 5,
                                                                             nr.perms = 100,
                                                                             is_GE = F,
                                                                             do_serial = F)
  
  
  stopCluster(cl)
}


if(F){
  
  n = 300
  N = rep(20, n)
  set.seed(1)
  X = rbinom(n = n,size = N,prob = runif(n = 100))
  Conf.level = 0.95
  
  alpha.one.sided = (1-Conf.level)/2
  X_rho = X[1:30]; N_rho = N[1:30]
  X_test = X[-(1:30)]; N_test = N[-(1:30)]
  n_test = length(X_rho)
  
  CI_param = mcleod.CI.estimation.parameters(theta_vec = seq(-2,2,0.5))
  
  bank = mcleod.CI.deconv.bank.constructor(N_test,CI_param)
  
  
  res_mcleod_holdout = mcleod(x.smp = X_rho,n.smp =N_rho,
                           a.limits = bank$CI_param$a.limits,
                           computational_parameters = bank$CI_param$comp_parameters,
                           prior_parameters = bank$CI_param$prior_parameters,
                           exact.numeric.integration = T)
  
  CDF_holdout = mcleod:::compute_medians_curve(res_mcleod_holdout)
  
  nr.perm = 200
  possible_rhos = c(0.1,0.2,0.3,0.4,0.5)
  q_for_rho_optimization = c(0.2,0.4,0.6,0.8)
  
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
    print(paste0('optimizing rho for q=',current_q,', which is equivalent to theta = ',theta_current_q,' in the holdout data'))
    
    q_lower_by_rho = rep(NA,length(possible_rhos))
    point_to_start_testing_GE = qbinom(p = alpha.one.sided,size = length(bank$N_vec),prob = current_q) / length(bank$N_vec)
    points_to_test_GE = which(bank$CI_param$q_vec<= point_to_start_testing_GE)
    points_to_test_GE = tail(points_to_test_GE,n = 3)
    
    if(length(points_to_test_GE)>=2){
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
                                               do_serial = T)
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
      print(paste0('selecting rho for LE'))
      for(index_rho in 1:length(possible_rhos)){
        
        current_rho = possible_rhos[index_rho]
        
        a_ind_LE = which(res_mcleod_data$parameters_list$a.vec >= theta_current_q + current_rho)
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
                                               do_serial = T)
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
  optimal_rho_by_theta_for_LE
  optimal_rho_by_theta_for_GE
  
 approxfun()
  
  GE.pval.grid = matrix(NA,ncol = bank$CI_param$n_theta,nrow = bank$CI_param$n_q,
                        dimnames = list(paste0('q = ',bank$CI_param$q_vec),paste0('theta = ',bank$CI_param$theta_vec)))
  LE.pval.grid = matrix(NA,ncol = bank$CI_param$n_theta,nrow = bank$CI_param$n_q,
                        dimnames = list(paste0('q = ',bank$CI_param$q_vec),paste0('theta = ',bank$CI_param$theta_vec)))
  
  
  res_mcleod_data = mcleod(x.smp = X_test,n.smp =N_test,
                           a.limits = bank$CI_param$a.limits,
                           computational_parameters = bank$CI_param$comp_parameters,
                           prior_parameters = bank$CI_param$prior_parameters,
                           exact.numeric.integration = T)
  
  res_mcleod_data$parameters_list$a.vec
  median_curve = mcleod:::compute_medians_curve(res_mcleod_data)
  
  cl <- makeCluster(8)
  registerDoParallel(cl)
  start.time = Sys.time()
  for(i in 1:bank$CI_param$n_q)
    for(j in 1:bank$CI_param$n_theta){
      #i=1;j=5
      print(paste0('testing GE at q_ind = ',i,'/',bank$CI_param$n_q,' , theta_ind = ',j,'/',bank$CI_param$n_theta))
      rho = 0.1
      a_ind_GE = which(res_mcleod_data$parameters_list$a.vec <= bank$CI_param$theta_vec[j] - rho)
      if(length(a_ind_GE) > 0 ){
        a_ind_GE = max(a_ind_GE)
      }else{
        a_ind_GE = 1
      }
      a_ind_LE = which(res_mcleod_data$parameters_list$a.vec >= bank$CI_param$theta_vec[j] + rho)
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
                                           alpha = 0.05,
                                           nr.perm = 100,
                                           do_check_vs_noiseless_case = T,
                                           do_check_vs_minimum_number_of_required_iterations = T,
                                           is_GE = T,
                                           do_serial = F)
      
      print(paste0('testing LE at q_ind = ',i,'/',bank$CI_param$n_q,
                   ' , theta_ind = ',j,'/',bank$CI_param$n_theta))
      LE.pval.grid[i,j] = mcleod.CI.PV.at_point(bank = bank,
                                                ind_theta = j,
                                                ind_q = i,
                                                CDF_value = median_curve[a_ind_LE],
                                                a_index = a_ind_LE,
                                                alpha = 0.05,
                                                nr.perm = 100,
                                                do_check_vs_noiseless_case = T,
                                                do_check_vs_minimum_number_of_required_iterations = T,
                                                is_GE = F,
                                                do_serial = F)
        
    }
  end.time = Sys.time()
  print(end.time-start.time)
  image(t(GE.pval.grid))
  image(t(LE.pval.grid))
  #View(LE.pval.grid)
    
}