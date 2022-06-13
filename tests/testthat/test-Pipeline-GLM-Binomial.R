test_that("Pipeline-GLM-Binomial", {
  # Generate Data:
  N = 30 #Number of draws per binomial observations
  K = 200 #Number of samples
  set.seed(1)
  covariates = matrix(rnorm(K*2,sd = 0.5),nrow = K) #Generate covariates
  colnames(covariates) = c('covariate 1','covariate 2')
  #define slopes:
  real_beta_1 = -1
  real_beta_2 = 1
  #sample
  x = rbinom(K,size = N,
             prob = inv.log.odds(rcauchy(K,location = 0,scale = 0.5) +
                                   real_beta_1*covariates[,1] + real_beta_2*covariates[,2]))
  n = rep(N,K)

  # Fit model
  res = mcleod(x, n, covariates = covariates)
  
  #Extract Coefficients
  coeffs = mcleod::results.covariate.coefficients.posterior(res)
  
  
  #Estimate gamma_i
  estimated.gamma_is = mcleod.posterior.estimates.random.effect(
    X = c(15,40),
    N = c(30,50),
    mcleod_res = res,
    covariates = rbind(c(   0,  0),
                       c(-0.5,0.5))
  )
  
  #Predictive intervals
  mcleod.predictive.interval(
    N = c(30,50),
    mcleod_res = res,
    covariates = rbind(c(   -2,  2),
                       c(  1,-1)),
    Interval.Coverage = 0.95)
  
  
  #Example with different init, slope prior and slope proposal:
  N = 30 # Number of draws per sample
  K = 300 #Number of samples
  set.seed(2)
  covariates = matrix(rexp(K),nrow = K) # exponentially distributed coefficients
  real_beta = -0.5 #the real value of the coefficient
  
  u = sample(c(0,1),size = K,replace = T)
  x = rbinom(K,size = N,
             prob =inv.log.odds(rnorm(K,-1+3*u,sd = 0.3) +
                                  real_beta*covariates))
  n = rep(N,K)
  
  
  #Generate object with parameters for covariate estimation
  coeffs_obj  = mcleod.covariates.estimation.parameters(
    beta_init = 0,
    beta_prior_sd = 1.5,
    proposal_sd = 0.025)
  
  
  comp_params = mcleod.computational.parameters(nr.gibbs = 3000)
  # In addition to covariates, pass coeffs_obj as an argument 
  res = mcleod(x,
               n,
               covariates = covariates,
               covariates_estimation_parameters = coeffs_obj,
               computational_parameters = comp_params
  )
  
  coeffs = mcleod::results.covariate.coefficients.posterior(res,
                                                            plot.posterior = F,
                                                            plot.MH.proposal.by.iteration = T)
  
  
  succeed(message = "Success for Pipeline-GLM-Binomial", info = NULL)
  
  
})
