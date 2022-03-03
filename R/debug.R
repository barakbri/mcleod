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
  compute_P_values_over_grid = T
  compute_CI_curves = F
  verbose = T
  ratio_holdout = 0.1
  
  ret = list()
  
  CI_param = mcleod.CI.estimation.parameters(theta_vec = seq(-2,2,0.5),
                                             do_serial = F,rho.estimation.perm = 50)
  
  alpha.one.sided = (1-CI_param$alpha.CI)/2
  n_holdout = ceiling(ratio_holdout * length(X))
  X_rho = X[1:n_holdout]; N_rho = N[1:n_holdout]
  X_test = X[-(1:n_holdout)]; N_test = N[-(1:n_holdout)]
  
  
  bank <<- mcleod.CI.deconv.bank.constructor(N_test,CI_param)
  
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
  rho_calibration_obj = mcleod.CI.rho.calibration.constructor(bank_original = bank,
                                      res_mcleod_holdout = res_mcleod_holdout,
                                      CDF_holdout = CDF_holdout,
                                      alpha.one.sided = alpha.one.sided,
                                      verbose = T,
                                      nr.perm = CI_param$rho.estimation.perm,
                                      possible_rhos = CI_param$rho.possible.values,
                                      q_for_rho_optimization = CI_param$rho.q_for_calibration)
  end.time = Sys.time()
  rho_calibration_obj$Elapsed_time = end.time-start.time
  if(verbose){
    print('rho calibration time')
    print(end.time-start.time)  
  }
  bank <<- rho_calibration_obj$bank
  rho_calibration_obj$bank <- NULL
  
  
  #rho_calibration_obj$rho_approx_fun_LE(seq(-2,2,0.25))
  #rho_calibration_obj$rho_approx_fun_GE(seq(-2,2,0.25))
  
  res_mcleod_data = mcleod(x.smp = X_test,n.smp =N_test,
                           a.limits = bank$CI_param$a.limits,
                           computational_parameters = bank$CI_param$comp_parameters,
                           prior_parameters = bank$CI_param$prior_parameters,
                           exact.numeric.integration = T)
  
  median_curve = mcleod:::compute_medians_curve(res_mcleod_data)
  
  
  if(compute_P_values_over_grid){
    start.time = Sys.time()
    
    
    pvalues_grid = mcleod:::compute_P_values_over_grid_function(
                               bank_original = bank,
                               rho_calibration_obj = rho_calibration_obj,
                               res_mcleod_data = res_mcleod_data,
                               median_curve = median_curve,
                               alpha.one.sided = alpha.one.sided,
                               verbose = verbose)

    
    
    end.time = Sys.time()
    pvalues_grid$Elapsed_time = end.time-start.time
    ret$pvalues_grid = pvalues_grid
    if(verbose){
      print('time for computing pvalues over grid of hypotheses')
      print(end.time-start.time)  
    }
    bank <<- pvalues_grid$bank
    pvalues_grid$bank <- NULL
    #image(t(pvalues_grid$GE.pval.grid))
    #image(t(pvalues_grid$LE.pval.grid))
  }
  
  
  
  if(compute_CI_curves){
    start.time = Sys.time()
    computed_curves = mcleod:::compute_CI_curves_function(
      bank_original = bank,
      res_mcleod_data = res_mcleod_data,
      rho_calibration_obj = rho_calibration_obj,
      median_curve = median_curve,
      alpha.one.sided = alpha.one.sided,
      verbose = verbose
    )
    end.time = Sys.time()
    computed_curves$Elapsed_time = end.time-start.time
    ret$computed_curves = computed_curves
    if(verbose){
      print('time for computing CI curves ')
      print(end.time-start.time)  
    }
    #bank <<- computed_curves$bank
    #computed_curves$bank <- NULL
  }
  
  ret$bank = bank
  ret$rho_calibration_obj = rho_calibration_obj
  class(ret) = mcleod:::CLASS.NAME.MCLEOD.CI
  
  if(!CI_param$do_serial){
    stopCluster(cl)
  }
  
  return(ret)
}