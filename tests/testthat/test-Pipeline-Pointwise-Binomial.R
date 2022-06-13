test_that("Pipeline-Pointwise-Binomial", {
  
  #Generate Data
  N = 30
  K = 300
  set.seed(1)
  u = sample(c(0,1),size = K,replace = T)
  x = rbinom(K,size = N,prob =inv.log.odds(rnorm(K,-1+3*u,sd = 0.3)))
  n = rep(N,K)
  
  #Fit model
  res = mcleod(x, n)
  
  
  #Extract results
  posterior_mixing_dist = mcleod.get.posterior.mixing.dist(res)
  
  plot.posterior(res)
  
  #Change parameters 1:
  res = mcleod(x, n, a.limits = c(-5,5))
  
  
  #Change parameters 2:
  prior_obj  = mcleod.prior.parameters( 
    prior.type =MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL, 
    Beta.Heirarchical.Levels = 6
  )
  
  res = mcleod(x, n, prior_parameters = prior_obj) 
  
  #Change parameters 3:
  prior_obj  = mcleod.prior.parameters(
    prior.type =MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET, 
    Two.Layer.Dirichlet.Intervals = 64, 
    Two.Layer.Dirichlet.Nodes.in.First.Layer = 8
  )
  
  res = mcleod(x, n, prior_parameters = prior_obj)
  
  #Change parameters 4:
  comp_obj = mcleod.computational.parameters(nr.gibbs = 500, #define the number of iter.s
                                             nr.gibbs.burnin = 250) 
  
  res = mcleod(x, n, computational_parameters = comp_obj) # pass object as argument
  
  #Estimate gamma_i
  estimated.log.odds = mcleod.posterior.estimates.random.effect(
    X = c(5,35),
    N = c(30,40),
    res)
  
  #Predictive interval
  mcleod.predictive.interval(
    N = c(30,60),
    mcleod_res = res,
    Interval.Coverage = 0.95)
  
  
  succeed(message = "Success for Pipeline-Pointwise-Binomial", info = NULL)
  
})
