test_that("Pipeline-Pointwise-Poisson", {
  # Generate data
  K = 200 
  set.seed(1)
  
  u = sample(c(0,1),size = K,replace = T) 
  x = rpois(K,lambda = exp(rnorm(K,2 + 3*u,0.5)) )

  # Fit model
  res = mcleod(x, n.smp = NULL,a.limits = c(-2,8),Noise_Type = MCLEOD.POISSON.ERRORS)
  
  #Plot
  plot.posterior(res)
  
  #estimate gamma i
  estimated.log.lambda = mcleod.posterior.estimates.random.effect(
    X = c(10,200),
    N = NULL,
    res)
  
  # Predictive interval
  mcleod.predictive.interval(
    N = NULL,
    mcleod_res = res,
    Interval.Coverage = 0.95)
  
  succeed(message = "Success for Pipeline-Pointwise-Poisson", info = NULL)
  
  
})
