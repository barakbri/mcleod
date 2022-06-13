test_that("Pipeline-GLM-SettingNonNormalPrior", {
  # Generate data
  N = 30 
  K = 200 
  set.seed(1)
  covariates = matrix(rnorm(K*2,sd = 0.5),nrow = K)
  real_beta_1 = -1 
  real_beta_2 = 1
  x = rbinom(K,size = N,prob = inv.log.odds(rnorm(K,0,sd = 1)
                                            + real_beta_1*covariates[,1] + real_beta_2*covariates[,2]))
  n = rep(N,K)

  
  #Configure prior
  
  beta_prior_points = seq(-5,5,0.01)
  beta_prior_probs = pcauchy(beta_prior_points[-1]) - 
    pcauchy(beta_prior_points[-length(beta_prior_points)])
  beta_prior_probs = beta_prior_probs/ sum(beta_prior_probs)
  
  beta_prior_probs = matrix(c(beta_prior_probs, beta_prior_probs),ncol = 2)
  
  coeffs_obj = mcleod.covariates.estimation.parameters(
    Manual_Prior_Values = beta_prior_points,
    Manual_Prior_Probs = beta_prior_probs)
  
  #Fit model:
  res = mcleod(x, n, covariates = covariates,
               covariates_estimation_parameters = coeffs_obj)
  
  #Obtain results:
  coeffs = mcleod::results.covariate.coefficients.posterior(res,plot.posterior = F)
  
  
  succeed(message = "Success for Pipeline-GLM-SettingNonNormalPrior", info = NULL)
  
})
