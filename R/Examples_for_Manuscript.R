library(doRNG)

DO_SAVE = T
DO_EFRON_DATA_ANALYSIS = F
DO_BETA_EXAMPLE_ANALYSIS = F

# Create the results dir
RESULTS_DIR = 'E:/CI_Results/'
if(!dir.exists(RESULTS_DIR)){
  dir.create(RESULTS_DIR)
}


#Run the Efron example
if(DO_EFRON_DATA_ANALYSIS){
  # load the data
  memory.limit(20000)
  x.vec = deconvolveR::surg[,2]
  n.vec = deconvolveR::surg[,1]
  
  # run prob based statistic
  set.seed(1)
  CI.res = mcleod.estimate.CI(x.vec = x.vec,
                              n.vec = n.vec,
                              a.max = 4,q_grid = seq(0.025,0.975,0.025),
                              CI.estimation.parameters = mcleod.CI.estimation.parameters(Nr.reps.for.each.n = 1,
                                                                                         nr.cores = detectCores(),
                                                                                         fraction.of.points.computed = 0.2,#0.5,
                                                                                         epsilon.nr.gridpoints = 2),
                              prior_param = mcleod.prior.parameters(prior.type = MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL,Beta.Heirarchical.Levels = 6),comp_param = mcleod.computational.parameters(nr.gibbs = 200,nr.gibbs.burnin = 100),
                              verbose = T
  )
  
  if(DO_SAVE)
    save(CI.res,file = paste0(RESULTS_DIR,'/CI_res.rdata'))
  #load(file = paste0(RESULTS_DIR,'/CI_res.rdata'))
  CI.res$Elapsed_Time_Parallel
  CI.res$Elapsed_Time_Overall
  pdf(file = paste0(RESULTS_DIR,'/Efron_Prob_based_x_is_P.pdf'),width = 10,height = 5)
    plot.mcleod.CI(CI.res,X_axis_as_Prob=T)
  dev.off()
  pdf(file = paste0(RESULTS_DIR,'/Efron_Prob_based_x_is_LogOdds.pdf'),width = 10,height = 5)
    plot.mcleod.CI(CI.res,X_axis_as_Prob=F)
  dev.off()
  
  
  CI.res_median_based = mcleod.estimate.CI.based.on.medians(x.vec = x.vec,
                                                            n.vec = n.vec,
                                                            conf.level = 0.95,
                                                            nr.perm = 200,
                                                            q_grid = seq(0.025,0.975,0.025),
                                                            a.max = 4,
                                                            shift.size.in.log.odds.scale = 0.25,
                                                            verbose = T,
                                                            prior_param = mcleod.prior.parameters(prior.type = MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL,Beta.Heirarchical.Levels = 6),
                                                            comp_param = mcleod.computational.parameters(nr.gibbs = 600,nr.gibbs.burnin = 300))
  
  pdf(file = paste0(RESULTS_DIR,'/Efron_Median_based_x_is_P.pdf'),width = 10,height = 5)
  plot.mcleod.CI(CI.res_median_based,X_axis_as_Prob=T)
  dev.off()
  pdf(file = paste0(RESULTS_DIR,'/Efron_Median_based_x_is_LogOdds.pdf'),width = 10,height = 5)
  plot.mcleod.CI(CI.res_median_based,X_axis_as_Prob=F)
  dev.off()
  CI.res_median_based$Elapsed_Time_Overall
  if(DO_SAVE)
    save(CI.res_median_based,file = paste0(RESULTS_DIR,'/CI_res_median_based.rdata'))
  #load(file = paste0(RESULTS_DIR,'/CI_res_median_based.rdata'))
  
}


if(DO_BETA_EXAMPLE_ANALYSIS){
  memory.limit(20000)

  sample_sizes = c(500,1000,2000)
  Objects_list = list()
  
  for(current_sample_size_index in 1:length(sample_sizes)){
    #current_sample_size_index = 1  
    print(paste0('current_sample_size_index = ',current_sample_size_index))
    objects_item = list()
    set.seed(1)
    K = sample_sizes[current_sample_size_index]
    n.vec = rep(3,K) #rep(2,K)
    p.vec = rbeta(n = K,2,2)
    x.vec = rbinom(K,size = n.vec,p.vec)
    
    
    
    # CI.res_beta_binomial = mcleod.estimate.CI(x.vec = x.vec,
    #                                           n.vec = n.vec,
    #                                           a.max = 4,
    #                                           seq(0.025,0.975,0.025),
    #                                           CI.estimation.parameters = mcleod.CI.estimation.parameters(Nr.reps.for.each.n = 1,
    #                                                                                                      nr.cores = detectCores(),
    #                                                                                                      fraction.of.points.computed = 1,
    #                                                                                                      epsilon.nr.gridpoints = 2),
    #                                           prior_param = mcleod.prior.parameters(prior.type = MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL,Beta.Heirarchical.Levels = 5),
    #                                           comp_param = mcleod.computational.parameters(nr.gibbs = 200,nr.gibbs.burnin = 100),
    #                                           verbose = T)
    # 
    # objects_item$CI.res_beta_binomial = CI.res_beta_binomial
    objects_item$CI.res_beta_binomial = NA
    
    
    CI.res_beta_binomial_based_on_medians = mcleod.estimate.CI.based.on.medians(x.vec = x.vec,
                                                                                n.vec = n.vec,
                                                                                conf.level = 0.95,
                                                                                nr.perm = 200,
                                                                                q_grid = seq(0.025,0.975,0.025),
                                                                                a.max = 4,
                                                                                shift.size.in.log.odds.scale = 0.25,
                                                                                verbose = T,
                                                                                prior_param = mcleod.prior.parameters(prior.type = MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL,Beta.Heirarchical.Levels = 6),
                                                                                comp_param = mcleod.computational.parameters(nr.gibbs = 1000,nr.gibbs.burnin = 500))
    
    objects_item$CI.res_beta_binomial_based_on_medians = CI.res_beta_binomial_based_on_medians
    
    Objects_list[[current_sample_size_index]] = objects_item
  }
  
  
  if(DO_SAVE)
    save(Objects_list,file = paste0(RESULTS_DIR,'/beta_example_obj_list_3.rdata'))
  #load(file = paste0(RESULTS_DIR,'/beta_example_obj_list_2.rdata'))
  
  if(F){
    library(rmutil)
    ind = 1
    CI.res_beta_binomial = Objects_list[[ind ]]$CI.res_beta_binomial
    CI.res_beta_binomial_based_on_medians = Objects_list[[ind ]]$CI.res_beta_binomial_based_on_medians
    # CI.res_beta_binomial$Elapsed_Time_Parallel
    # CI.res_beta_binomial$Elapsed_Time_Overall
    CI.res_beta_binomial_based_on_medians$Elapsed_Time_Overall
    # plot.mcleod.CI(CI.res_beta_binomial)
    # plot.mcleod.CI(CI.res_beta_binomial,X_axis_as_Prob = F)
    
    pdf(file = paste0(RESULTS_DIR,'/BetaBinomial_',sample_sizes[ind],'_median_based_X_is_P.pdf'),width = 10,height = 5)
    plot.mcleod.CI(CI.res_beta_binomial_based_on_medians,X_axis_as_Prob=T)
    abline(v=(0.65),col =  'blue')
    abline(v=(0.75),col =  'blue')
    abline(h=(0.4),col =  'blue')
    beta_CDF_X = inv.log.odds(CI.res_beta_binomial_based_on_medians$a.vec) 
    beta_CDF_Y = pbeta(beta_CDF_X,shape1 = 2,shape2 = 2)
    lines(beta_CDF_X,beta_CDF_Y,col = 'darkgreen',lwd =2)
    dev.off()
    
    pdf(file = paste0(RESULTS_DIR,'/BetaBinomial_',sample_sizes[ind],'_median_based_X_is_LogOdds.pdf'),width = 10,height = 5)
    plot.mcleod.CI(CI.res_beta_binomial_based_on_medians,X_axis_as_Prob=F)
    abline(v=log.odds(0.65),col =  'blue')
    abline(v=log.odds(0.75),col =  'blue')
    abline(h=(0.4),col =  'blue')
    lines(log.odds(beta_CDF_X),beta_CDF_Y,col = 'darkgreen',lwd =2)
    dev.off()
    
  }
  
  
  
}
