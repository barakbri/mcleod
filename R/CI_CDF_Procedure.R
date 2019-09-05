library(hash)
library(doRNG)
library(doParallel)
library(parallel)


CLASS.NAME.MCLEOD.CI = 'mcleod.CI.obj'
CLASS.NAME.MCLEOD.CI.PARAMETERS = 'mcleod.CI.obj.parameters'


#' Title
#'
#' @param Nr.reps.for.each.n 
#' @param nr.cores 
#' @param epsilon.nr.gridpoints 
#' @param fraction.of.points.computed 
#'
#' @return
#' @export
#'
#' @examples
mcleod.CI.estimation.parameters = function(Nr.reps.for.each.n = 1,
                                           nr.cores = detectCores() - 1,
                                           epsilon.nr.gridpoints = 2,
                                           fraction.of.points.computed = 1){
  
  # need to add checks
  ret = list()
  class(ret) = CLASS.NAME.MCLEOD.CI.PARAMETERS
  ret$Nr.reps.for.each.n = Nr.reps.for.each.n
  ret$nr.cores = nr.cores
  ret$epsilon.nr.gridpoints = epsilon.nr.gridpoints
  ret$fraction.of.points.computed = fraction.of.points.computed
  return(ret)
}

#' Title
#'
#' @param x.vec 
#' @param n.vec 
#' @param q_grid 
#' @param a.max 
#' @param comp_param 
#' @param prior_param 
#' @param CI.estimation.parameters 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
mcleod.estimate.CI = function(x.vec,
                              n.vec,
                              q_grid = seq(0.1,0.9,0.1),
                              a.max = c(3),
                              comp_param = mcleod.computational.parameters(nr.gibbs = 200,
                                                                           nr.gibbs.burnin = 100),
                              prior_param = mcleod.prior.parameters(prior.type = MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET,
                                                                    Two.Layer.Dirichlet.Intervals = 64,
                                                                    Two.Layer.Dirichlet.Nodes.in.First.Layer = 8),
                              CI.estimation.parameters = mcleod.CI.estimation.parameters(),
                              verbose = T){
  
  Start.time.overall = Sys.time()
  K = length(n.vec)  
  M=8
  
  if(class(CI.estimation.parameters) != CLASS.NAME.MCLEOD.CI.PARAMETERS)
    stop('CI.estimation.parameters must be the result of mcleod.CI.estimation.parameters(...)')
  
  NR_reps_for_each_n = CI.estimation.parameters$Nr.reps.for.each.n
  NR.CORES = CI.estimation.parameters$nr.cores
  epsilon.nr.gridpoints = CI.estimation.parameters$epsilon.nr.gridpoints
  fraction.of.points.computed = CI.estimation.parameters$fraction.of.points.computed
  
  #Part I: functions for generating P_k_i based on precomputed values:
  
  sorted_n_data = sort(unique(n.vec))
  
  n_to_P_k_i_generator_index = hash()
  for(i in 1:length(sorted_n_data))
    n_to_P_k_i_generator_index[ as.character(sorted_n_data[i]) ] = i
  
  generate_P_k_i_generator_list = function(n_vec = 1:20){
    P_k_i_generator_list = list()
    for(n_i in 1:length(n_vec)){
      n = n_vec[n_i]
      n.smp = rep(n,n+1)
      x.smp = c(0:n)
      
      P_k_i_for_n = mcleod(x.smp = x.smp,n.smp = n.smp,
                           a.limits = c(-a.max,a.max),
                           computational_parameters = mcleod.computational.parameters(nr.gibbs = 2,
                                                                                      nr.gibbs.burnin = 1),
                           prior_parameters = prior_param)
      P_k_i_generator_list[[n_i]] = P_k_i_for_n$additional$original_stat_res$p_k_i
    }
    return(P_k_i_generator_list)
  }
  
  generator_list = suppressWarnings(generate_P_k_i_generator_list(n_vec = sorted_n_data)) 
  
  generate_P_k_i = function(x_to_generate_for){
    P_k_i_generated = matrix(NA,nrow = K,ncol = dim(generator_list[[1]])[2])
    for(k in 1:K){
      P_k_i_generated[k,] = (generator_list[[  n_to_P_k_i_generator_index[[as.character(n.vec[k]) ]] ]])[ x_to_generate_for[k] + 1, ]
    }
    return(P_k_i_generated)
  }
  
  #e.g., you can now call
  #x_to_generate_for = deconvolveR::surg[,2]
  #generated_P_k_i = generate_P_k_i(x_to_generate_for)
  
  
  #Part II: generate deconv model for the original data
  if(verbose){
    cat(paste0(' - Computing deconvolution estimate for original data.\n\r'))
  }
  
  mcleod_for_data = mcleod(x.smp = x.vec,n.smp = n.vec,a.limits = c(-a.max,a.max),
                           computational_parameters = comp_param,
                           prior_parameters = prior_param)
  
  
  #Part III: generete decov models for resampled data at grid points, for LE based shift
  if(verbose){
    cat(paste0(' - Computing deconvolution estimate for permuted, worst-case data generations.\n\r'))
  }
  Start.time = Sys.time()
  a.vec = mcleod_for_data$parameters_list$a.vec
  epsilon  = epsilon.nr.gridpoints*(a.vec[2] - a.vec[1])
  LE_grid = a.vec[-1] - epsilon
  
  mcleods = list()
  
  cl <- makeCluster(NR.CORES)
  registerDoParallel(cl)
  
  for(grid_point_i in 2:(length(a.vec))){ 
    if(verbose)
      cat(paste0('    Computing for grid point ',grid_point_i-1,'\n\r'))
    
    current_a.vec = a.vec[grid_point_i]
    current_a.vec_minus_epsilon = LE_grid[grid_point_i - 1]
    
    worker_function_at_current_a_i = function(n_at_a.vec_grid_point){
      mcleods_at_current_n_and_a_i = list()
      for(rep_i in 1:NR_reps_for_each_n){
        
        
        #compute a single realization of the data, for n items
        theta_sample = sample(
          c(rep(current_a.vec_minus_epsilon,
                n_at_a.vec_grid_point),
            rep(M,
                K-n_at_a.vec_grid_point)
          ),
          size = K
        )
        
        X_sampled = rbinom(n = K,size = n.vec,prob = inv.log.odds(theta = theta_sample))
        generated_P_k_i = generate_P_k_i(X_sampled)
        
        temp_mcleod = mcleod(x.smp = X_sampled,n.smp = n.vec,
                             a.limits = c(-a.max,a.max),
                             input_P_k_i = generated_P_k_i,
                             computational_parameters = comp_param,
                             prior_parameters = prior_param)
        
        n_smp_posterior_mat = t(temp_mcleod$additional$original_stat_res$n_smp)
        mcleods_at_current_n_and_a_i[[rep_i]] = n_smp_posterior_mat[-(1:comp_param$nr.gibbs.burnin),]
      }
      return(mcleods_at_current_n_and_a_i)
    }
    
    parallel_res <- foreach(n_at_a.vec_grid_point=0:K, .options.RNG=1234) %dorng% {
      library(mcleod);
      points_to_compute = seq(from = 0,to = K,by = ceiling(1/fraction.of.points.computed))
      
      if(!(n_at_a.vec_grid_point %in% points_to_compute)){
        return(NA)  
      }
      
      worker_function_at_current_a_i(n_at_a.vec_grid_point)
    }
    
    mcleods_at_currents_a_i = parallel_res
    
    mcleods[[grid_point_i]] = mcleods_at_currents_a_i
    
  }
  stopCluster(cl)
  
  End.time = Sys.time()
  Elapsed_Time_Parallel = End.time - Start.time
  if(verbose){
    cat(paste0(' - Cluster stopped. Total time in parallel: \n\r'))
    print(Elapsed_Time_Parallel)
  }
  #Part IV: compute matrix of P-values
  
  
  pval_LE = matrix(NA,nrow = length(q_grid),ncol = length(a.vec)-1)
  pval_GE = matrix(NA,nrow = length(q_grid),ncol = length(a.vec)-1)
  data_n_smp_posterior_mat = t(mcleod_for_data$additional$original_stat_res$n_smp)
  data_n_smp_posterior_mat = (data_n_smp_posterior_mat)[-(1:comp_param$nr.gibbs.burnin),]
  
  compute_LE_stat = function(q,grid_point_i,n_smp_mat){
    Left_n_smp = n_smp_mat[,1:grid_point_i,drop = F]
    Right_n_smp = n_smp_mat[,-c(1:grid_point_i),drop = F]
    pbeta_vec = pbeta(q,shape1 = 1+ apply(Left_n_smp,1,sum),shape2 = 1+ apply(Right_n_smp,1,sum))
    return(mean(pbeta_vec))
  }
  
  
  if(verbose){
    cat(paste0(' - Computing P-values for determining Upper CI : \n\r'))
  }
  for(ai in 2:(length(a.vec)-1)){
    if(verbose){
      cat(paste0('    Computing at grid point: ',ai,'\n\r'))
    }
    
    grid_point_i = ai 
    
    for(qi in 1:length(q_grid)){
      
      q_0 = q_grid[qi]
      binom_weights = dbinom(0:(K),size = K,prob = q_0)  
      data_based_estimate = compute_LE_stat(q = q_0,
                                            grid_point_i = grid_point_i,
                                            n_smp_mat = data_n_smp_posterior_mat)
      indicators_matrix = matrix(NA,nrow = K+1, ncol = NR_reps_for_each_n)
      for(k in 1:(K+1)){
        for(r in 1:NR_reps_for_each_n){
          if(!is.na((mcleods[[grid_point_i]])[[k]])){ # check that the point wasnt skipped
            current_sampled_posterior = (mcleods[[grid_point_i]])[[k]][[r]]
            
            indicators_matrix[k,r] = compute_LE_stat(q = q_0,
                                                     grid_point_i = grid_point_i,
                                                     n_smp_mat = current_sampled_posterior)  
          }else{
            binom_weights[k] = NA
          }
          
        }  
      }
      
      binom_weights = binom_weights/sum(binom_weights,na.rm = T) # renormalize, in case we had k's skipped
      
      indicators_matrix = 1*(indicators_matrix >= data_based_estimate)
      pval_at_point = sum(apply(indicators_matrix,1,mean,na.rm=T) * binom_weights,na.rm = T)
      pval_LE[qi,ai] = pval_at_point
    }
  }
  
  if(verbose){
    cat(paste0(' - Computing P-values for determining Lower CI : \n\r'))
  }
  for(ai in 1:(length(a.vec)-2)){
    if(verbose){
      cat(paste0('    Computing at grid point: ',ai,'\n\r'))
    }
    grid_point_i = ai 
    
    for(qi in 1:length(q_grid)){
      
      q_0 = q_grid[qi]
      binom_weights = dbinom(K:0,size = K,prob = 1- q_0)  
      data_based_estimate = 1.0 - compute_LE_stat(q = q_0,
                                                  grid_point_i = grid_point_i,
                                                  n_smp_mat = data_n_smp_posterior_mat)
      indicators_matrix = matrix(NA,nrow = K+1, ncol = NR_reps_for_each_n)
      for(k in 1:(K+1)){
        for(r in 1:NR_reps_for_each_n){
          if(!is.na(mcleods[[length(a.vec) - grid_point_i]])[[(K+2) - k]]){ # check that the point wasnt skipped
            current_sampled_posterior = (mcleods[[length(a.vec) - grid_point_i]])[[(K+2) - k]][[r]]
            
            indicators_matrix[k,r] = compute_LE_stat(q = 1-q_0,
                                                     grid_point_i = length(a.vec) - grid_point_i,
                                                     n_smp_mat = current_sampled_posterior)
          }else{
            binom_weights[k] = NA
          }
        }  
      }
      indicators_matrix = 1*(indicators_matrix >= data_based_estimate)
      pval_at_point = sum(apply(indicators_matrix,1,mean) * binom_weights,na.rm = T) #need to debug
      pval_GE[qi,ai] = pval_at_point
    }
  }
  
  #Part V: Correct edge cases, for tests that are too discrete to have power:
  # if(verbose){
  #   cat(paste0(' - Requiring monotonicity of P-values, for grid points with (low-theta, high CDF) or (high-theta, low-CDF)\n\r'))
  # }
  # pval_LE_original = pval_LE
  # pval_GE_original = pval_GE
  # 
  # for(i in 2:ncol(pval_LE)){
  #   for(j in 2:nrow(pval_LE)){
  #     if(pval_LE[j,i] > pval_LE[j-1,i]){
  #       pval_LE[j,i] = pval_LE[j-1,i]
  #     }
  #   }
  # }
  # 
  # for(i in 1:(ncol(pval_GE)-1)){
  #   for(j in (nrow(pval_GE)-1):1){
  #     if(pval_GE[j,i] > pval_GE[j+1,i]){
  #       pval_GE[j,i] = pval_GE[j+1,i]
  #     }
  #   }
  # }
  
  ret = list()
  class(ret) = CLASS.NAME.MCLEOD.CI
  ret$a.vec = a.vec
  ret$q_grid = q_grid
  ret$LE_grid = LE_grid
  ret$pval_LE = pval_LE
  ret$pval_GE = pval_GE
  ret$epsilon = epsilon
  ret$K = K
  ret$Elapsed_Time_Parallel = Elapsed_Time_Parallel
  ret$mcleod_for_data = mcleod_for_data
  # ret$pval_LE_original = pval_LE_original
  # ret$pval_GE_original = pval_GE_original
  ret$Elapsed_Time_Overall = Sys.time() - Start.time.overall
  return(ret)
  
}


#' Title
#'
#' @param mcleod.CI.obj 
#' @param sig.level 
#'
#' @return
#' @export
#'
#' @examples
plot.mcleod.CI=function(mcleod.CI.obj, sig.level = c(0.05,0.1)){
  
  if(class(mcleod.CI.obj) != CLASS.NAME.MCLEOD.CI){
    stop('mcleod.CI.obj must be a result returned from mcleod.estimate.CI')
  }
  
  a.vec = mcleod.CI.obj$a.vec
  q_grid  = mcleod.CI.obj$q_grid
  LE_grid  = mcleod.CI.obj$LE_grid
  pval_LE  = mcleod.CI.obj$pval_LE
  epsilon  = mcleod.CI.obj$epsilon
  pval_GE  = mcleod.CI.obj$pval_GE
  
  dt_plot_LE = matrix(NA,nrow = (length(q_grid) *(length(a.vec)-1)),ncol = 3)
  colnames(dt_plot_LE) = c('theta.LE','q','PV')
  pointer = 1
  for(qi in 1:length(q_grid)){
    for(ai in 2:(length(a.vec)-1)){
      dt_plot_LE[pointer,] = c(LE_grid[ai-1],q_grid[qi],pval_LE[qi,ai]) #c(a.vec[ai+2],q_grid[qi],pval_LE[qi,ai])#
      pointer = pointer +1
    }
  }
  dt_plot_LE = as.data.frame(dt_plot_LE)
  
  dt_plot_GE = matrix(NA,nrow = (length(q_grid) *(length(a.vec)-1)),ncol = 3)
  colnames(dt_plot_GE) = c('theta.GE','q','PV')
  pointer = 1
  for(qi in 1:length(q_grid)){
    for(ai in 2:(length(a.vec)-1)){
      dt_plot_GE[pointer,] = c(LE_grid[ai-1] + 2*epsilon,q_grid[qi],pval_GE[qi,ai]) #c(a.vec[ai+2],q_grid[qi],pval_GE[qi,ai])#
      pointer = pointer +1
    }
  }
  dt_plot_GE = as.data.frame(dt_plot_GE)
  min.a = min(a.vec)
  max.a = max(a.vec)

  nr.iter.burnin = mcleod.CI.obj$mcleod_for_data$parameters_list$nr.gibbs.burnin
  posterior_mean = (t(mcleod.CI.obj$mcleod_for_data$additional$original_stat_res$pi_smp))[-(1:nr.iter.burnin),]
  posterior_mean =  cbind(rep(0,nrow(posterior_mean)),t(apply(posterior_mean,1,cumsum)))
  posterior_mean = apply(posterior_mean,2,mean)
  dt_post_mean = data.frame(theta = a.vec,posterior_mean = posterior_mean)
  
  
  v <- ggplot() + geom_contour(mapping =  aes(theta.LE, q, z = PV),
                               data =  dt_plot_LE,breaks = sig.level) +
    geom_contour(mapping =  aes(theta.GE, q, z = PV),
                 data =  dt_plot_GE,breaks = sig.level) +
    geom_line(mapping = aes(x = theta,y=posterior_mean),col= 'red',data = dt_post_mean)
  print(v  +ylim(c(0,1)) +xlim(c(min.a,max.a)) +xlab('theta'))
}

if(F){
  memory.limit(20000)
  x.vec = deconvolveR::surg[,2]
  n.vec = deconvolveR::surg[,1]

  set.seed(1)
  CI.res = mcleod.estimate.CI(x.vec = x.vec,
                              n.vec = n.vec,
                              a.max = 4,
                              CI.estimation.parameters = mcleod.CI.estimation.parameters(Nr.reps.for.each.n = 1,
                                                                                         nr.cores = detectCores(),
                                                                                         fraction.of.points.computed = 0.33,
                                                                                         epsilon.nr.gridpoints = 2),
                              prior_param = mcleod.prior.parameters(
                                prior.type = MCLEOD.PRIOR.TYPE.TWO.LAYER.DIRICHLET,
                                Two.Layer.Dirichlet.Intervals = 48,
                                Two.Layer.Dirichlet.Nodes.in.First.Layer = 4),
                              verbose = T
                              )

  CI.res$Elapsed_Time_Parallel
  CI.res$Elapsed_Time_Overall
  plot.mcleod.CI(CI.res)
}






