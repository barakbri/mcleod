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
  N = rpois(n = n,lambda = 20)
  shape_1 = 2
  shape_2 = 2
  set.seed(1)
  X = rbinom(n = n,size = N,prob = rbeta(n = n,shape1 = shape_1,shape2 = shape_2))
  CI_param = mcleod.CI.estimation.parameters(theta_vec = seq(-2.5,2.5,0.25),
                                             q_vec = seq(0.05,0.95,0.05),
                                             do_serial = F,
                                             rho.estimation.perm = 50)
  
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
  lines(oracle_x,pbeta(q = oracle_x,shape1 = shape_1,shape2 = shape_2),col = 'blue',lwd =1.5)
  
}