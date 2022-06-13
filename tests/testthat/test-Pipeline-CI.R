test_that("Pipeline-CI", {
  skip(message = "SKIPPED")
  # Generate Data
  n = 500
  N = rep(20,n) 
  set.seed(1)
  p = inv.log.odds(rnorm(n,-2,0.5)+3*rbinom(n,1,0.3))
  
  X = rbinom(n = n,size = N,prob = p)
  p_true = inv.log.odds(rnorm(1E6,-2,0.5)+3*rbinom(1E6,1,0.3))
  
  
  # Estimate CIs
  CI.est.res = mcleod.estimate.CI(X = X, N = N,
                                  CI_param = mcleod.CI.estimation.parameters(rho.set.value = 0.25))

  # Plot
  plot.mcleod.CI(mcleod.CI.obj = CI.est.res) #one command to plots CIs
  
  # Extract CI
  dt_CIs = mcleod.get.CIs.mixing.dist(CI.est.res)

  # Different definitions - construct object  
  CI_param = mcleod.CI.estimation.parameters(
    theta_vec = seq(-4,4,0.1),
    q_vec = seq(0.05,0.95,0.025),
    alpha.CI = 0.9,
    rho.possible.values = seq(0.1,0.3,0.1))
  
  # we pass CI_param as an object, and also set ratio_holdout=0.1 
  CI.est.res = mcleod.estimate.CI(X = X,
                                  N = N,
                                  CI_param = CI_param,
                                  ratio_holdout = 0.1)
  
  
  set.seed(1)
  
  #estimate CI for single q
  CI.res = mcleod.estimate.CI.single.q(X = X,N = N,q = 0.5,
                 CI_param = mcleod.CI.estimation.parameters(theta_vec = seq(-3,3,0.05),
                 rho.set.value = 0.25), 
                 verbose = F)
  
  
  
  
  #estimate CI for single theta
  CI.res = mcleod.estimate.CI.single.theta(X = X, N = N,
                                           theta = 0,
                                           CI_param = mcleod.CI.estimation.parameters(q_vec = seq(0.05,0.95,0.01),
                                                                                      rho.set.value = 0.25) )
  
  #Generate a similar dataset, 200 samples
  n = 200
  N = rep(20,n) 
  set.seed(1)
  p = inv.log.odds(rnorm(n,-2,0.5)+3*rbinom(n,1,0.3))
  X = rbinom(n = n,size = N,prob = p)
  
  # Run with compute_P_values_over_grid = T
  CI.est.res = mcleod.estimate.CI(X = X, N = N,
                                  compute_P_values_over_grid = T,
                                  compute_CI_curves = F,
                                  CI_param = mcleod.CI.estimation.parameters(rho.set.value = 0.25),
                                  verbose = T)
  
  succeed(message = "Success for Pipeline-CI", info = NULL)
  
  
})
