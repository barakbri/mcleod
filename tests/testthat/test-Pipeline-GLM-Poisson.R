test_that("Pipeline-GLM-Poisson", {
  # Generate Data
  
  K = 200 
  set.seed(2)
  
  covariates = matrix(rexp(K,rate = 2),nrow = K) 
  real_beta = 0.5
  u = sample(c(0,1),size = K,replace = T) 
  extrinsic_size = runif(n = K,1,100)
  offset = log(extrinsic_size) 
  
  x = rpois(K,
            lambda = extrinsic_size * exp(rnorm(K,2 + 3*u,0.5) + real_beta* covariates)
  )
  
  #set the number of iterations
  comp_obj = mcleod.computational.parameters(nr.gibbs = 3000,nr.gibbs.burnin = 500)
  
  #fit model
  res = mcleod(x, n.smp = NULL, 
               a.limits = c(-2,8), 
               computational_parameters = comp_obj,
               covariates = covariates, 
               Noise_Type = MCLEOD.POISSON.ERRORS, 
               offset_vec = offset 
  )
  
  # get coefficients
  coeffs = mcleod::results.covariate.coefficients.posterior(res)
  
  # get random effect distribution:
  posterior_dist_gamma = mcleod.get.posterior.mixing.dist(res)
  
  # estimate gamma_i's
  estimated.gamma_is = mcleod.posterior.estimates.random.effect(
    X = x[1:3],
    N = NULL,
    mcleod_res = res,
    covariates = covariates[1:3,,drop=F],
    offset_vec = offset[1:3])
  
  
  #Predictive intervals
  mcleod.predictive.interval(
    N = NULL,
    covariates = matrix(c(0, -1),ncol = 1),
    offset_vec = c(log(1) , log(3)),
    mcleod_res = res,
    Interval.Coverage = 0.95)
  

  succeed(message = "Success for Pipeline-GLM-Poisson", info = NULL)
  
  
})
