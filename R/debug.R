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
  
  CI_param = mcleod.CI.estimation.parameters(theta_vec = seq(-2,2,0.5),
                                             do_serial = F,rho.estimation.perm = 50)
  
  alpha.one.sided = (1-CI_param$alpha.CI)/2
  X_rho = X[1:30]; N_rho = N[1:30]
  X_test = X[-(1:30)]; N_test = N[-(1:30)]
  n_test = length(X_rho)
  
  
  bank = mcleod.CI.deconv.bank.constructor(N_test,CI_param)
  
  cl <- NULL
  if(!CI_param$do_serial){
    cl <- makeCluster(CI_param$nr.cores)
    registerDoParallel(cl)
  }
  
  
  res_mcleod_holdout = mcleod(x.smp = X_rho,n.smp =N_rho,
                           a.limits = bank$CI_param$a.limits,
                           computational_parameters = bank$CI_param$comp_parameters,
                           prior_parameters = bank$CI_param$prior_parameters,
                           exact.numeric.integration = T)
  
  CDF_holdout = mcleod:::compute_medians_curve(res_mcleod_holdout)
  
  start.time = Sys.time()
  rho_calibration_obj = mcleod.CI.rho.calibration.constructor(bank = bank,
                                      res_mcleod_holdout = res_mcleod_holdout,
                                      CDF_holdout = CDF_holdout,
                                      alpha.one.sided = alpha.one.sided,
                                      verbose = T,
                                      nr.perm = CI_param$rho.estimation.perm,
                                      possible_rhos = CI_param$rho.possible.values,
                                      q_for_rho_optimization = CI_param$rho.q_for_calibration)
  end.time = Sys.time()
  print('rho calibration time')
  print(end.time-start.time)
  
  #rho_calibration_obj$rho_approx_fun_LE(seq(-2,2,0.25))
  #rho_calibration_obj$rho_approx_fun_GE(seq(-2,2,0.25))
  
  GE.pval.grid = matrix(NA,ncol = bank$CI_param$n_theta,nrow = bank$CI_param$n_q,
                        dimnames = list(paste0('q = ',bank$CI_param$q_vec),paste0('theta = ',bank$CI_param$theta_vec)))
  LE.pval.grid = matrix(NA,ncol = bank$CI_param$n_theta,nrow = bank$CI_param$n_q,
                        dimnames = list(paste0('q = ',bank$CI_param$q_vec),paste0('theta = ',bank$CI_param$theta_vec)))
  
  res_mcleod_data = mcleod(x.smp = X_test,n.smp =N_test,
                           a.limits = bank$CI_param$a.limits,
                           computational_parameters = bank$CI_param$comp_parameters,
                           prior_parameters = bank$CI_param$prior_parameters,
                           exact.numeric.integration = T)
  
  median_curve = mcleod:::compute_medians_curve(res_mcleod_data)
  
  
  start.time = Sys.time()
  for(i in 1:bank$CI_param$n_q)
    for(j in 1:bank$CI_param$n_theta){
      #i=1;j=5
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
  end.time = Sys.time()
  print('time for computing pvalues over grid of hypotheses')
  print(end.time-start.time)
  
  if(!CI_param$do_serial){
    stopCluster(cl)
  }
  
  image(t(GE.pval.grid))
  image(t(LE.pval.grid))
  #View(LE.pval.grid)
    
}