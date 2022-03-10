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
  N = rep(20,n)
  shape_1 = 2
  shape_2 = 2
  set.seed(1)
  X = rbinom(n = n,size = N,prob = rbeta(n = n,shape1 = shape_1,shape2 = shape_2))
  
  CI_param = mcleod.CI.estimation.parameters(theta_vec = seq(-3,3,0.25),
                                             q_vec = seq(0.02,0.98,0.04),
                                             do_serial = F,
                                             rho.estimation.perm = 50,
                                             nr.perms = 200,alpha.CI = 0.95,
                                             rho.possible.values = seq(0,1.0,0.1),
                                             rho.q_for_calibration = seq(0.1,0.9,0.1))
  
  CI.est.res = mcleod.estimate.CI(X = X,
                                  N = N,
                                  CI_param = CI_param,
                                  ratio_holdout = 0.2,
                                  compute_P_values_over_grid = F,
                                  compute_CI_curves = T,
                                  verbose = T)
  
  

  plot.mcleod.CI(mcleod.CI.obj = CI.est.res,
                 X_axis_as_Prob = T,
                 add_CI_curves_on_top_of_plot = F)

  
  oracle_x = seq(0.01,0.99,0.01)
  lines(oracle_x,pbeta(q = oracle_x,shape1 = shape_1,shape2 = shape_2),col = 'blue',lwd =1.5)

  plot.posterior(CI.est.res$res_mcleod_holdout,plot_only_point_estimate = T) + ggtitle('CDF for holdout data')
  #plot rho:
  x = seq(-3,3,0.1)
  plot(x = x,y= CI.est.res$rho_calibration_obj$rho_approx_fun_LE(x),ylim = c(0,1),
       col = 'red',type = 'l',xlab = 'theta',ylab = 'rho')
  lines(x=x,y=CI.est.res$rho_calibration_obj$rho_approx_fun_GE(x),col = 'blue')
  
  plot(x = x,y= CI.est.res$rho_calibration_obj$rho_approx_fun_LE_non_smoothed(x),ylim = c(0,1),
       col = 'red',type = 'l',xlab = 'theta',ylab = 'rho')
  lines(x=x,y=CI.est.res$rho_calibration_obj$rho_approx_fun_GE_non_smoothed(x),col = 'blue')
  
}

if(F){
  
  n = 500
  N = rep(20,n) 
  shape_1 = 2
  shape_2 = 2
  set.seed(1)
  p = inv.log.odds(rnorm(n,-2,0.5)+3*rbinom(n,1,0.3))
  
  X = rbinom(n = n,size = N,prob = p)
  p_true = inv.log.odds(rnorm(1E6,-2,0.5)+3*rbinom(1E6,1,0.3))
  
  CI_param = mcleod.CI.estimation.parameters(theta_vec = seq(-2.5,2.5,0.25),
                                             q_vec = seq(0.05,0.95,0.05),
                                             do_serial = F,
                                             rho.estimation.perm = 50,
                                             rho.possible.values = seq(0,1.0,0.1),
                                             rho.q_for_calibration = seq(0.1,0.9,0.1))
  
  CI.est.res = mcleod.estimate.CI(X = X,
                                  N = N,
                                  CI_param = CI_param,
                                  ratio_holdout = 0.1,
                                  compute_P_values_over_grid = F,
                                  compute_CI_curves = T,
                                  verbose = T)
  
  
  
  plot.mcleod.CI(mcleod.CI.obj = CI.est.res,
                 X_axis_as_Prob = T,
                 add_CI_curves_on_top_of_plot = F)
  
  
  oracle_x = seq(0.1,0.9,0.01)
  
  lines(oracle_x,(ecdf(p_true))(oracle_x),col = 'blue',lwd =1.5)
  
  #plot rho:
  x = seq(-3,3,0.1)
  plot(x = x,y= CI.est.res$rho_calibration_obj$rho_approx_fun_LE(x),ylim = c(0,1),col = 'red',type = 'l')
  lines(x=x,y=CI.est.res$rho_calibration_obj$rho_approx_fun_GE(x),col = 'blue')
  
  
  
  CI.est.res.fixed.rho.list = list()
  fixed_rho_values = seq(0,1,0.2)
  for(i in 1:length(fixed_rho_values)){
    print(paste0('Solving for fixed rho = ',fixed_rho_values[i]))
    CI_param_rho_fixed = mcleod.CI.estimation.parameters(theta_vec = seq(-2.5,2.5,0.25),
                                               q_vec = seq(0.05,0.95,0.05),
                                               do_serial = F,
                                               rho.estimation.perm = 50,
                                               rho.set.value = fixed_rho_values[i])
    if(i==1){
      perm.object = NULL
    }else{
      perm.object = CI.est.res.fixed.rho.list[[i-1]]
    }
    CI.est.res.fixed.rho.list[[i]] = mcleod.estimate.CI(X = X,
                                                  N = N,
                                                  CI_param = CI_param_rho_fixed,
                                                  ratio_holdout = 0.1,
                                                  compute_P_values_over_grid = F,
                                                  compute_CI_curves = T,
                                                  verbose = F,
                                                  Use_Existing_Permutations_From_Object = CI.est.res)  
  }
  
  
  par(mfrow=c(2,3))
  for(i in 1:length(fixed_rho_values)){
    plot.mcleod.CI(mcleod.CI.obj = CI.est.res,
                   X_axis_as_Prob = T,
                   add_CI_curves_on_top_of_plot = F,title = paste0('rho = ',fixed_rho_values[i]))
    
    
     oracle_x = seq(0.1,0.9,0.01)
    
     lines(oracle_x,(ecdf(p_true))(oracle_x),col = 'blue',lwd =1.5)
    # 
    
    plot.mcleod.CI(mcleod.CI.obj = CI.est.res.fixed.rho.list[[i]],
                   X_axis_as_Prob = T,
                   add_CI_curves_on_top_of_plot = T,CI_curves_color = 'orange')  
  }
  par(mfrow=c(1,1))
  
}


if(F){
  
  n = 200
  N = rep(50,n)
  shape_1 = 2
  shape_2 = 2
  set.seed(1)
  X = rbinom(n = n,size = N,prob = rbeta(n = n,shape1 = shape_1,shape2 = shape_2))
  
  theta_0 = c(-1,1)
  
  CI_param = mcleod.CI.estimation.parameters(theta_vec = theta_0,
                                             q_vec = seq(0.02,0.98,0.04),
                                             do_serial = F,
                                             rho.estimation.perm = 50,
                                             nr.perms = 200,alpha.CI = 0.9,
                                             rho.possible.values = seq(0,1.0,0.1),
                                             rho.q_for_calibration = seq(0.1,0.9,0.1),
                                             rho.theta_for_calibration = theta_0)
  
  
  CI.est.res = mcleod.estimate.CI(X = X,
                                  N = N,
                                  CI_param = CI_param,
                                  ratio_holdout = 0.2,
                                  compute_P_values_over_grid = F,
                                  compute_CI_curves = T,
                                  verbose = T)
  #CI.est.res$pvalues_grid
  c(CI.est.res$computed_curves$q_star_GE,CI.est.res$computed_curves$q_star_LE)
  
  
}




if(F){
  
  n = 200
  N = rep(50,n)
  shape_1 = 2
  shape_2 = 2
  set.seed(1)
  X = rbinom(n = n,size = N,prob = rbeta(n = n,shape1 = shape_1,shape2 = shape_2))
  
  q_0 = 0.5#c(0.3,0.7)
  
  CI_param = mcleod.CI.estimation.parameters(theta_vec = seq(-3,3,0.25),
                                             q_vec = q_0,
                                             do_serial = F,
                                             rho.estimation.perm = 50,
                                             nr.perms = 200,alpha.CI = 0.95,
                                             rho.possible.values = seq(0,1.0,0.1),
                                             rho.q_for_calibration = q_0)
  
  
  CI.est.res = mcleod.estimate.CI(X = X,
                                  N = N,
                                  CI_param = CI_param,
                                  ratio_holdout = 0.2,
                                  compute_P_values_over_grid = F,
                                  compute_CI_curves = T,
                                  verbose = T)
  
  c(CI.est.res$computed_curves$q_star_GE,
    CI.est.res$computed_curves$q_star_LE)
  
  
  # CI.est.res$pvalues_grid$GE.pval.grid
  # CI.est.res$pvalues_grid$LE.pval.grid
  
}