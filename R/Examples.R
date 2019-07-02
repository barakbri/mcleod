
#simple example, two covariates
if(F){
  N = 30
  K = 200
  set.seed(1)
  covariates = matrix(rnorm(K*2,sd = 0.5),nrow = K)
  real_beta_1 = -1
  real_beta_2 = 1
  x = rbinom(K,size = N,prob = inv.log.odds(rnorm(K,0,sd = 1) + real_beta_1*covariates[,1] + real_beta_2*covariates[,2]))
  n = rep(N,K)
  model_dt = data.frame(c = x,nc = n-x)
  model_dt = cbind(model_dt,covariates)
  model <- glm(cbind(c,nc) ~.,family=binomial,data=model_dt)
  model$coefficients[-1]
  
  res = mcleod(x, n, prior_parameters = mcleod.prior.parameters(),
                                        a.limits = c(-4,4),
                                        covariates = covariates,
                                        computational_parameters = mcleod.computational.parameters(nr.gibbs = 500),
                                        covariates_estimation_parameters = mcleod.covariates.estimation.parameters(beta_init = model$coefficients[-1])
  )
  res$additional$original_stat_res$elapsed_secs
  
  mcleod::results.covariate.coefficients.posterior(res,plot.MH.proposal.by.iteration = T)
  
  plot.posterior(res)
}




# bimodal over-dispersion, exp covariates
if(F){
  N = 30
  K = 300
  set.seed(2)
  covariates = matrix(rexp(K),nrow = K)
  hist(covariates[,1])
  real_beta = -0.5
  u = sample(c(0,1),size = K,replace = T)
  x = rbinom(K,size = N,prob =inv.log.odds(rnorm(K,-1+3*u,sd = 0.3) + real_beta*covariates))
  n = rep(N,K)
  hist(x/n)
  plot(ecdf(x/n))
  model_dt = data.frame(c = x,nc = n-x)
  model_dt = cbind(model_dt,covariates)
  model <- glm(cbind(c,nc) ~.,family=binomial,data=model_dt)
  model$coefficients[-1]
  
  res = mcleod(x, n, prior_parameters = mcleod.prior.parameters(),
                                        a.limits = c(-4,4),
                                        covariates = covariates,
                                        computational_parameters = mcleod.computational.parameters(nr.gibbs = 500,nr.gibbs.burnin = 250),
                                        covariates_estimation_parameters = mcleod.covariates.estimation.parameters(beta_init = model$coefficients[-1])
  )
  
  res$additional$original_stat_res$elapsed_secs
  
  mcleod::results.covariate.coefficients.posterior(res,plot.MH.proposal.by.iteration = T)
  
  plot.posterior(res)
}

#Poisson example
if(F){
  
  K = 200
  set.seed(2)
  covariates = matrix(rexp(K,rate = 2),nrow = K)
  real_beta = 0.5
  u = sample(c(0,1),size = K,replace = T)
  x = rpois(K,lambda = exp(rnorm(K,2 + 3*u,0.5) + real_beta* covariates) )
  
  res = mcleod(x, n.smp = NULL,
               prior_parameters = mcleod.prior.parameters(prior.type = MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET,
                                                          Two.Layer.Dirichlet.Nodes.in.First.Layer = 16),
                                                                                a.limits = c(-2,8),
               computational_parameters = mcleod.computational.parameters(nr.gibbs = 1000,nr.gibbs.burnin = 500),
                                        covariates = covariates,
                                        Noise_Type = MCLEOD.POISSON.ERRORS
  )
  
  res$additional$original_stat_res$elapsed_secs
  
  mcleod::results.covariate.coefficients.posterior(res,plot.MH.proposal.by.iteration = T)
  
  plot.posterior(res)
   
  
}
