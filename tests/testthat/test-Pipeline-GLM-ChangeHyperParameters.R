test_that("Pipeline-GLM-ChangeHyperParameters", {
  # Generate data
  
  set.seed(1)
  
  K = 2000 # Number of samples
  n.vec = rep(20,K) #number of draws for each sample
  
  # function for sampling a beta binomial (2,2)
  generate_sample_Beta_2_2 = function(K){ 
    p.vec = rbeta(n = K,2,2)
    x.vec = rbinom(K,size = n.vec,p.vec)
    return(x.vec)
  }
  
  #Generate date
  x.vec = generate_sample_Beta_2_2(K)
  
  #run sampler:
  mcleod_res = mcleod(x.vec,n.vec,a.limits = c(-4,4),
                      # HBeta prior with 5 levels     
                      prior_param = mcleod.prior.parameters(
                        prior.type = MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL,
                        Beta.Heirarchical.Levels = 5),
                      # Set number of MCMC samples
                      computational_parameters = mcleod.computational.parameters(
                        nr.gibbs = 2000,nr.gibbs.burnin = 1000)
  )
  
  mcleod::plot.posterior(mcleod_res)
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #Example how to change hyper-parameters - Polya tree
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #For, example, how to change all hyper parameters to 2
  # We form two matrices, for $\alpha^{L}$ and $\alpha^{U}$.
  # In each matrix, the [l,i] entry gives the appropriate hyper parameter
  # for the ith node on the lth level.
  
  L = 5 # number of levels
  alpha_L = matrix(NA,nrow = L,ncol = 2^(L-1)) # hyper-parameter matrices 
  alpha_U = matrix(NA,nrow = L,ncol = 2^(L-1))
  
  # Fill all entries to a value of 2
  for(l in 1:L){
    for(k in 1:2^(l-1)){
      alpha_L[l,k] = 2
      alpha_U[l,k] = 2
    }
  }
  
  # we pass the hyper-paramers to the mcleod.prior.parameters(...) function
  mcleod_res = mcleod(x.vec,n.vec,
                      a.limits = c(-4,4),
                      #pass prior hyper parameters
                      prior_param = mcleod.prior.parameters(
                        prior.type = MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL,
                        Beta.Heirarchical.Levels = L,
                        #NOTE: THESE ARE THE MATRICES FOR THE HYPER-PARAMETERS:
                        Prior_Hyper_Parameters_BetaH_L = alpha_L,
                        Prior_Hyper_Parameters_BetaH_U = alpha_U),
                      # set the number of iterations
                      computational_parameters = mcleod.computational.parameters(
                        nr.gibbs = 2000,
                        nr.gibbs.burnin = 1000)
  )
  
  mcleod::plot.posterior(mcleod_res)
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #Example how to change hyper-parameters - 2 layer dirichlet tree
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Set tree structure
  Nodes_in_first_layer = 8
  Total_nr_nodes = 64
  
  #set matrix for hyper parameters.
  Two_Level_Dirichlet_Tree_Hyperparameters = matrix(NA,nrow = 2,ncol = Total_nr_nodes)
  
  #Dirichlet at top of tree is Dirichlet(2,2,2,2,2,2,2,2)
  Two_Level_Dirichlet_Tree_Hyperparameters[1,1:Nodes_in_first_layer] = 2
  
  # At the middle level we have 8 Dirichlet random variables.
  # We set all of them to be Dirichlet(2,2,2,2,2,2,2,2) as well.
  
  Two_Level_Dirichlet_Tree_Hyperparameters[2,1:Total_nr_nodes] = 2
  
  # we pass the hyper-paramers to the mcleod.prior.parameters(...) function
  mcleod_res = mcleod(x.vec,n.vec,
                      a.limits = c(-4,4),
                      exact.numeric.integration = T,
                      #Define prior
                      prior_param = mcleod.prior.parameters(
                        #Define as 2-layer Dirichlet tree
                        prior.type = MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET,
                        #set tree structure
                        Two.Layer.Dirichlet.Intervals = Total_nr_nodes,
                        Two.Layer.Dirichlet.Nodes.in.First.Layer =  Nodes_in_first_layer,
                        #NOTE: HERE WE PASS THE HYPER-PARAMETERS
                        Prior_Hyper_Parameters_2LDT = Two_Level_Dirichlet_Tree_Hyperparameters),
                      
                      #Set the number of iterations
                      computational_parameters = mcleod.computational.parameters(
                        nr.gibbs = 2000,
                        nr.gibbs.burnin = 1000)
  )
  
  mcleod::plot.posterior(mcleod_res)

  succeed(message = "Success for Pipeline-GLM-ChangeHyperParameters", info = NULL)
  
  
})
