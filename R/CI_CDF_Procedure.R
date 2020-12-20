library(hash)
library(doRNG)
library(doParallel)
library(parallel)


CLASS.NAME.MCLEOD.CI = 'mcleod.CI.obj'
CLASS.NAME.MCLEOD.CI.MEDIAN.BASED.STATISTIC = 'mcleod.CI.obj.median.based.stat'
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
  
  library(doRNG)
  library(doParallel)
  library(parallel)
  
  # need to add checks
  ret = list()
  class(ret) = CLASS.NAME.MCLEOD.CI.PARAMETERS
  ret$Nr.reps.for.each.n = Nr.reps.for.each.n
  ret$nr.cores = nr.cores
  ret$epsilon.nr.gridpoints = epsilon.nr.gridpoints
  ret$fraction.of.points.computed = fraction.of.points.computed
  return(ret)
}


generate_P_k_i_matrix_cache = function(n.vec,a.max,prior_param,comp_param){
  
  sorted_n_data = sort(unique(n.vec))
  
  n_to_P_k_i_generator_index = hash()
  for(i in 1:length(sorted_n_data))
    n_to_P_k_i_generator_index[ as.character(sorted_n_data[i]) ] = i
  
  generate_P_k_i_generator_list = function(n_vec = 1:20){
    P_k_i_generator_list = hash() #list()
    for(n_i in 1:length(n_vec)){
      n = n_vec[n_i]
      n.smp = rep(n,n+1)
      x.smp = c(0:n)
      
      P_k_i_for_n = mcleod(x.smp = x.smp,n.smp = n.smp,
                           a.limits = c(-a.max,a.max),
                           computational_parameters = mcleod.computational.parameters(nr.gibbs = 2,
                                                                                      nr.gibbs.burnin = 1,
                                                                                      integration_step_size = comp_param$integration_step_size,
                                                                                      Fast.Gamma.Used = comp_param$Fast.Gamma.Used,
                                                                                      Fast.Gamma.Bank.Size = comp_param$Fast.Gamma.Bank.Size),
                           prior_parameters = prior_param)
      P_k_i_generator_list[[as.character(n_i)]] = P_k_i_for_n$additional$original_stat_res$p_k_i
    }
    return(P_k_i_generator_list)
  }
  
  
  generator_list = suppressWarnings(generate_P_k_i_generator_list(n_vec = sorted_n_data)) 
  
  generate_P_k_i = function(x_to_generate_for,n.dim,n.vec){
    K=length(n.vec)
    P_k_i_generated = matrix(NA,nrow = K,ncol = n.dim)
    for(k in 1:K){
      P_k_i_generated[k,] = (generator_list[[ 
        as.character(
          n_to_P_k_i_generator_index[[   as.character(n.vec[k])  ]]
        )
        ]])[ x_to_generate_for[k] + 1, ]
    }
    return(P_k_i_generated)
  }
  
  ret = list()
  ret$sorted_n_data = sorted_n_data
  ret$n_to_P_k_i_generator_index = n_to_P_k_i_generator_index
  ret$generate_P_k_i_generator_list = generate_P_k_i_generator_list
  ret$generator_list = generator_list
  ret$generate_P_k_i = generate_P_k_i
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
  gen_object = generate_P_k_i_matrix_cache(n.vec,a.max,prior_param,comp_param)
  sorted_n_data = gen_object$sorted_n_data
  n_to_P_k_i_generator_index = gen_object$n_to_P_k_i_generator_index
  generate_P_k_i_generator_list = gen_object$generate_P_k_i_generator_list
  generator_list = gen_object$generator_list
  generate_P_k_i = gen_object$generate_P_k_i
  
  #e.g., you can now call:
  #x_to_generate_for = deconvolveR::surg[,2]
  #generated_P_k_i = generate_P_k_i(x_to_generate_for,length(a.vec)-1,n.vec)
  
  
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
    #points per n * length N.vec * fraction of points computed * length(a.vec) / cores * time of data run *2
    
    estiamted_time_to_run_in_seconds = mcleod_for_data$additional$original_stat_res$elapsed_secs *
      CI.estimation.parameters$fraction.of.points.computed *
      length(n.vec) *
      CI.estimation.parameters$Nr.reps.for.each.n *
      (length(mcleod_for_data$parameters_list$a.vec)-1) / NR.CORES *2 *1.1
    #The 1.1 is due to some overhead...
    cat(paste0('    Estimated time to run worst-case data generations : ',estiamted_time_to_run_in_seconds,' seconds.\n\r'))
      
  }
  Start.time = Sys.time()
  a.vec = mcleod_for_data$parameters_list$a.vec
  epsilon  = epsilon.nr.gridpoints*(a.vec[2] - a.vec[1])
  LE_grid = a.vec[-1] - epsilon
  
  mcleods = list()
  
  compute_LE_stat = function(q,grid_point_i,n_smp_mat){
    Left_n_smp = n_smp_mat[,1:grid_point_i,drop = F]
    Right_n_smp = n_smp_mat[,-c(1:grid_point_i),drop = F]
    pbeta_vec = pbeta(q,shape1 = 1+ apply(Left_n_smp,1,sum),shape2 = 1+ apply(Right_n_smp,1,sum))
    return(mean(pbeta_vec))
  }
  
  cl <- makeCluster(NR.CORES)
  registerDoParallel(cl)
  
  worker_configurations = expand.grid(
                                    grid_point_i = 2:(length(a.vec)-1),
                                    n_at_a.vec_grid_point=0:K
                                      )
  worker_function_at_current_a_i = function(n_at_a.vec_grid_point,grid_point_i){

      current_a.vec = a.vec[grid_point_i]
      current_a.vec_minus_epsilon = LE_grid[grid_point_i - 1]

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
        generated_P_k_i = generate_P_k_i(X_sampled,length(a.vec)-1,n.vec)

        temp_mcleod = mcleod(x.smp = X_sampled,n.smp = n.vec,
                             a.limits = c(-a.max,a.max),
                             input_P_k_i = generated_P_k_i,
                             computational_parameters = comp_param,
                             prior_parameters = prior_param)

        n_smp_posterior_mat = t(temp_mcleod$additional$original_stat_res$n_smp)
        n_smp_posterior_mat = n_smp_posterior_mat[-(1:comp_param$nr.gibbs.burnin),]
        mcleods_at_current_n_and_a_i[[rep_i]] = rep(NA,length(q_grid))
        for(qi in 1:length(q_grid)){
          
          q_0 = q_grid[qi]
          
          mcleods_at_current_n_and_a_i[[rep_i]][qi] = compute_LE_stat(q = q_0,
                                                      grid_point_i = grid_point_i,
                                                      n_smp_mat = n_smp_posterior_mat)
        }
      }
      return(mcleods_at_current_n_and_a_i)
    }
  
  # job_id = 64
  # n_at_a.vec_grid_point = worker_configurations$n_at_a.vec_grid_point[job_id]
  # grid_point_i = worker_configurations$grid_point_i[job_id]
  # points_to_compute = seq(from = 0,to = K,by = ceiling(1/fraction.of.points.computed))
  # if(!(n_at_a.vec_grid_point %in% points_to_compute)){
  #   return(NA)
  # }
  # worker_function_at_current_a_i(n_at_a.vec_grid_point,grid_point_i)
  # 
  parallel_res <- foreach(job_id = 1:nrow(worker_configurations), .options.RNG=1234) %dorng% {
    library(mcleod);
    n_at_a.vec_grid_point = worker_configurations$n_at_a.vec_grid_point[job_id]
    grid_point_i = worker_configurations$grid_point_i[job_id]
    
    points_to_compute = seq(from = 0,to = K,by = ceiling(1/fraction.of.points.computed))
    
    if(!(n_at_a.vec_grid_point %in% points_to_compute)){
      return(NA)
    }

    worker_function_at_current_a_i(n_at_a.vec_grid_point,grid_point_i)
  }
  
  #collect the results
  for(grid_point_i in 2:(length(a.vec))){
    mcleods[[grid_point_i]] = list()  
  }
  for(job_id in 1:length(parallel_res)){
    n_at_a.vec_grid_point = worker_configurations$n_at_a.vec_grid_point[job_id]
    grid_point_i = worker_configurations$grid_point_i[job_id]
     
    mcleods[[grid_point_i]][[n_at_a.vec_grid_point + 1]] = parallel_res[[job_id]]
  }
  
  # for(grid_point_i in 2:(length(a.vec))){ 
  #   
  #   worker_function_at_current_a_i = function(n_at_a.vec_grid_point,grid_point_i){
  #     
  #     current_a.vec = a.vec[grid_point_i]
  #     current_a.vec_minus_epsilon = LE_grid[grid_point_i - 1]
  #     
  #     mcleods_at_current_n_and_a_i = list()
  #     for(rep_i in 1:NR_reps_for_each_n){
  #       
  #       
  #       #compute a single realization of the data, for n items
  #       theta_sample = sample(
  #         c(rep(current_a.vec_minus_epsilon,
  #               n_at_a.vec_grid_point),
  #           rep(M,
  #               K-n_at_a.vec_grid_point)
  #         ),
  #         size = K
  #       )
  #       
  #       X_sampled = rbinom(n = K,size = n.vec,prob = inv.log.odds(theta = theta_sample))
  #       generated_P_k_i = generate_P_k_i(X_sampled,length(a.vec)-1,n.vec)
  #       
  #       temp_mcleod = mcleod(x.smp = X_sampled,n.smp = n.vec,
  #                            a.limits = c(-a.max,a.max),
  #                            input_P_k_i = generated_P_k_i,
  #                            computational_parameters = comp_param,
  #                            prior_parameters = prior_param)
  #       
  #       n_smp_posterior_mat = t(temp_mcleod$additional$original_stat_res$n_smp)
  #       mcleods_at_current_n_and_a_i[[rep_i]] = n_smp_posterior_mat[-(1:comp_param$nr.gibbs.burnin),]
  #     }
  #     return(mcleods_at_current_n_and_a_i)
  #   }
  #   
  #   parallel_res <- foreach(n_at_a.vec_grid_point=0:K, .options.RNG=1234) %dorng% {
  #     library(mcleod);
  #     points_to_compute = seq(from = 0,to = K,by = ceiling(1/fraction.of.points.computed))
  #     
  #     if(!(n_at_a.vec_grid_point %in% points_to_compute)){
  #       return(NA)  
  #     }
  #     
  #     worker_function_at_current_a_i(n_at_a.vec_grid_point)
  #   }
  #   
  #   mcleods_at_currents_a_i = parallel_res
  #   
  #   mcleods[[grid_point_i]] = mcleods_at_currents_a_i
  #   
  # }
  #stopCluster(cl)
  
  End.time = Sys.time()
  Elapsed_Time_Parallel = End.time - Start.time
  if(verbose){
    cat(paste0(' - Cluster stopped. Total time in parallel for worst case: \n\r'))
    print(Elapsed_Time_Parallel)
  }
  #Part IV: compute matrix of P-values
  
  
  pval_LE = matrix(NA,nrow = length(q_grid),ncol = length(a.vec)-1)
  pval_GE = matrix(NA,nrow = length(q_grid),ncol = length(a.vec)-1)
  data_n_smp_posterior_mat = t(mcleod_for_data$additional$original_stat_res$n_smp)
  data_n_smp_posterior_mat = (data_n_smp_posterior_mat)[-(1:comp_param$nr.gibbs.burnin),]
  
  
  
  if(verbose){
    cat(paste0(' - Computing P-values for determining Upper CI \n\r'))
  }
  
  
  worker_PV_computer_LE = function(grid_point_i){
    ret_vec  =rep(NA,length(q_grid))
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
            #current_sampled_posterior = (mcleods[[grid_point_i]])[[k]][[r]]
            
            #indicators_matrix[k,r] = compute_LE_stat(q = q_0,
            #                                         grid_point_i = grid_point_i,
            #                                         n_smp_mat = current_sampled_posterior)  
            
            indicators_matrix[k,r] = (mcleods[[grid_point_i]])[[k]][[r]][qi]
          }else{
            binom_weights[k] = NA
          }
          
        }  
      }
      
      binom_weights = binom_weights/sum(binom_weights,na.rm = T) # renormalize, in case we had k's skipped
      
      indicators_matrix = 1*(indicators_matrix >= data_based_estimate)
      pval_at_point = sum(apply(indicators_matrix,1,mean,na.rm=T) * binom_weights,na.rm = T)
      ret_vec[qi] = pval_at_point
      #pval_LE[qi,ai] = pval_at_point
    }
    return(ret_vec)
  }
  
  
  #cl <- makeCluster(NR.CORES)
  #registerDoParallel(cl)
  parallel_res_LE <- foreach(ai = 2:(length(a.vec)-1), .options.RNG=1234) %dorng% {
    library(mcleod);
    grid_point_i = ai 
    return(list(ai=ai,PV = worker_PV_computer_LE(grid_point_i)))
  }
  
  #stopCluster(cl) 
  
  for(j in 1:length(parallel_res_LE)){
    pval_LE[,parallel_res_LE[[j]]$ai] = parallel_res_LE[[j]]$PV
  }
  
  if(verbose){
    cat(paste0(' - Computing P-values for determining Lower CI \n\r'))
  }
  
  worker_PV_computer_GE = function(grid_point_i){
    ret_vec  =rep(NA,length(q_grid))
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
            # current_sampled_posterior = (mcleods[[length(a.vec) - grid_point_i]])[[(K+2) - k]][[r]]
            
            # indicators_matrix[k,r] = compute_LE_stat(q = 1-q_0,
            #                                          grid_point_i = length(a.vec) - grid_point_i,
            #                                          n_smp_mat = current_sampled_posterior)
            
            indicators_matrix[k,r] = (mcleods[[length(a.vec) - grid_point_i]])[[(K+2) - k]][[r]][length(q_grid)-qi+1]
          }else{
            binom_weights[k] = NA
          }
        }  
      }
      indicators_matrix = 1*(indicators_matrix >= data_based_estimate)
      pval_at_point = sum(apply(indicators_matrix,1,mean) * binom_weights,na.rm = T) #need to debug
      #pval_GE[qi,ai] = pval_at_point
      ret_vec[qi] = pval_at_point
    }
    return(ret_vec)
  }
  
  #cl <- makeCluster(NR.CORES)
  #registerDoParallel(cl)
  parallel_res_GE <- foreach(ai = 1:(length(a.vec)-2), .options.RNG=1234) %dorng% {
    library(mcleod);
    grid_point_i = ai 
    return(list(ai=ai,PV = worker_PV_computer_GE(grid_point_i)))
  }
  
  stopCluster(cl) 
  
  for(j in 1:length(parallel_res_GE)){
    pval_GE[,parallel_res_GE[[j]]$ai] = parallel_res_GE[[j]]$PV
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
#' @param conf.level 
#' @param X_axis_as_Prob 
#'
#' @return
#' @export
#'
#' @examples
plot.mcleod.CI=function(mcleod.CI.obj, conf.level = c(0.95),X_axis_as_Prob = T){
  
  OBJECT_IS_CLASS.NAME.MCLEOD.CI = (class(mcleod.CI.obj) == CLASS.NAME.MCLEOD.CI)
  OBJECT_IS_CLASS.NAME.MCLEOD.CI.MEDIAN.BASED.STATISTIC = (class(mcleod.CI.obj) == CLASS.NAME.MCLEOD.CI.MEDIAN.BASED.STATISTIC)
  
  if(!OBJECT_IS_CLASS.NAME.MCLEOD.CI & !OBJECT_IS_CLASS.NAME.MCLEOD.CI.MEDIAN.BASED.STATISTIC){
    stop('mcleod.CI.obj must be a result returned from mcleod.estimate.CI or mcleod.estimate.CI.based.on.medians')
  }
  
  shift.for.draw = NA
  if(OBJECT_IS_CLASS.NAME.MCLEOD.CI.MEDIAN.BASED.STATISTIC)
    shift.for.draw = mcleod.CI.obj$shift.size.in.log.odds.scale
  if(OBJECT_IS_CLASS.NAME.MCLEOD.CI)
    shift.for.draw = mcleod.CI.obj$epsilon
  
  function_for_Pval_LE = function(q_ind,ai){
    return(mcleod.CI.obj$pval_LE[q_ind,ai])
  }
  
  function_for_Pval_GE = function(q_ind,ai){
    return(mcleod.CI.obj$pval_GE[q_ind,ai])
  }
  
  if(OBJECT_IS_CLASS.NAME.MCLEOD.CI){
    nr_points_trim_GE = sum(apply(is.na(mcleod.CI.obj$pval_GE),2,sum)>0)
    nr_points_trim_LE = sum(apply(is.na(mcleod.CI.obj$pval_LE),2,sum)>0)  
  }
  
  
  
  
  for(conf.level.ind in 1:length(conf.level)){
    #conf.level.ind = 1
    if(OBJECT_IS_CLASS.NAME.MCLEOD.CI){
      curve_obj = compute.mcleod.CI.curve(mcleod.CI.obj$mcleod_for_data,
                                          mcleod.CI.obj$q_grid,
                                          function_for_Pval_LE,
                                          function_for_Pval_GE,
                                          nr_points_trim_LE,
                                          nr_points_trim_GE,
                                          conf.level[conf.level.ind])  
    }
    if(OBJECT_IS_CLASS.NAME.MCLEOD.CI.MEDIAN.BASED.STATISTIC){
      curve_obj = mcleod.CI.obj$curve_obj
    }
    
    x_axis = curve_obj$a.vec
    x_axis_GE = x_axis + shift.for.draw
    x_axis_LE = x_axis - shift.for.draw
    x_axis_label = 'theta'
    if(X_axis_as_Prob){
      x_axis = inv.log.odds(x_axis)
      x_axis_GE = inv.log.odds(x_axis_GE)
      x_axis_LE = inv.log.odds(x_axis_LE)
      x_axis_label = 'P'
    }
    if(conf.level.ind == 1)
      plot(x_axis,curve_obj$median_cumulative_pi_smp_for_data, col =  'red',type = 'b',pch = 20,xlab = x_axis_label,ylab = 'CDF')
    lines(x_axis,curve_obj$q_star_LE,col = 'black')
    lines(x_axis,curve_obj$q_star_GE,col = 'black')
  }
  
  
}


compute_medians_curve = function(mcleod_for_data){
  pi_smp_for_data = (t(mcleod_for_data$additional$original_stat_res$pi_smp))
  pi_smp_for_data = pi_smp_for_data[-(1:mcleod_for_data$parameters_list$nr.gibbs.burnin),]
  cumulative_pi_smp_for_data = t(apply(pi_smp_for_data,1,cumsum))
  median_cumulative_pi_smp_for_data = c(0,apply(cumulative_pi_smp_for_data,2,median))
  return(median_cumulative_pi_smp_for_data)
}



compute.mcleod.CI.curve=function(mcleod_for_data,
                                 q_grid,
                                 function_for_Pval_LE,
                                 function_for_Pval_GE,
                                 nr_points_trim_LE,
                                 nr_points_trim_GE,
                                 conf.level = c(0.95),verbose = F){
  #mcleod.CI.obj = CI.res
  # if(class(mcleod.CI.obj) != CLASS.NAME.MCLEOD.CI){
  #   stop('mcleod.CI.obj must be a result returned from mcleod.estimate.CI')
  # }
  # 
  a.vec = mcleod_for_data$parameters_list$a.vec
  #q_grid  = mcleod.CI.obj$q_grid
  #LE_grid  = mcleod.CI.obj$LE_grid
  #pval_LE  = mcleod.CI.obj$pval_LE
  #epsilon  = mcleod.CI.obj$epsilon
  #pval_GE  = mcleod.CI.obj$pval_GE
  
  
  median_cumulative_pi_smp_for_data = compute_medians_curve(mcleod_for_data)
  
  
  
  maximal_point_for_GE = rep(NA, length(a.vec))
  minimal_point_for_LE = rep(NA, length(a.vec))
  for(ai in 2:(length(a.vec)-1)){
    maximal_point_for_GE[ai] = max(which(q_grid<=median_cumulative_pi_smp_for_data[ai]),1)
    minimal_point_for_LE[ai] = min(which(q_grid>=median_cumulative_pi_smp_for_data[ai]),length(q_grid))
  }
  
  
  i_star_GE = rep(NA,length(a.vec))
  q_star_GE = rep(NA,length(a.vec))
  for(ai in 2:(length(a.vec)-1-nr_points_trim_GE)){
    if(verbose){
      print(paste0('Computing GE confidence intervals, log.odds = ',a.vec[ai]))
    }
    #ai = 2
    if(ai == 2){
      N_21_GE = 1
    }else{
      N_21_GE = min(i_star_GE[ai-1] + 1,length(q_grid))
    }
    N_22_GE = maximal_point_for_GE[ai]
    
    ind_to_select_from = (max(N_21_GE,1)):N_22_GE
    rejected_ind = rep(NA,length(ind_to_select_from))
    for(u in 1:length(rejected_ind)){
      rejected_ind[u] = function_for_Pval_GE(ind_to_select_from[u],ai)<= (1-conf.level)/2  
      if(is.na(rejected_ind[u])) #need to handle bounds and shifts better
        rejected_ind[u] = F  
      if(u==1 & rejected_ind[u] == F ){
        rejected_ind = rep(F,length(ind_to_select_from))
        break
      }
        
    }
    if(ai == 2 & rejected_ind[1] == F){
      i_star_GE[ai] = 0
    }else if (ai == 2){
      i_star_GE[ai] = ind_to_select_from[max(which(rejected_ind))]
    }
    
    if(ai >2 & rejected_ind[1] == F){
      i_star_GE[ai] = i_star_GE[ai-1]
    }else if (ai >2){
      i_star_GE[ai] = ind_to_select_from[max(which(rejected_ind))]
    }
    
    if(i_star_GE[ai] == 0){
      q_star_GE[ai] = 0  
    }else{
      q_star_GE[ai] = q_grid[i_star_GE[ai]]   
    }
    
  }
  
  #maximum_point_to_consider_for_LE = max(which(is.finite(minimal_point_for_LE)))
  #maximum_point_to_consider_for_LE = min(maximum_point_to_consider_for_LE,length(a.vec)-1 )
  i_star_LE = rep(NA,length(a.vec))
  q_star_LE = rep(NA,length(a.vec))
  for(ai in (length(a.vec)-1):(nr_points_trim_LE +1)){
    if(verbose){
      print(paste0('Computing LE confidence intervals, log.odds = ',a.vec[ai]))
    }
    #ai = nr_points_trim_LE +1
    #ai = 31#32
    if(ai == length(a.vec)-1){
      N_21_LE = length(q_grid)
    }else{
      N_21_LE = max(i_star_LE[ai+1] - 1,1)
    }
    N_22_LE = minimal_point_for_LE[ai]
    
    ind_to_select_from = N_22_LE:N_21_LE # (max(N_21_LE,1)):N_22_LE
    
    rejected_ind = rep(NA,length(ind_to_select_from))
    
    for(u in length(ind_to_select_from):1){
      rejected_ind[u] = function_for_Pval_LE(ind_to_select_from[u],ai)<= (1-conf.level)/2
      if(is.na(rejected_ind[u])) #need to handle bounds and shifts better
        rejected_ind[u] = F
      if(u==length(ind_to_select_from) & rejected_ind[u] == F )
        break
    }
    
    if(ai == length(a.vec)-1 & rejected_ind[length(rejected_ind)] == F){
      i_star_LE[ai] = length(q_grid)+1
    }else if (ai == length(a.vec)-1){
      i_star_LE[ai] = ind_to_select_from[min(which(rejected_ind))]
    }
    
    if(ai <length(a.vec)-1 & rejected_ind[length(rejected_ind)] == F){
      i_star_LE[ai] = i_star_LE[ai+1]
    }else if (ai <length(a.vec)-1){
      i_star_LE[ai] = ind_to_select_from[min(which(rejected_ind))]
    }
    
    if(i_star_LE[ai] == length(q_grid)+1){
      q_star_LE[ai] = 1  
    }else{
      q_star_LE[ai] = q_grid[i_star_LE[ai]]   
    }
    
  }
  
  ret = list(a.vec = a.vec,q_star_LE = q_star_LE,q_star_GE = q_star_GE,median_cumulative_pi_smp_for_data = median_cumulative_pi_smp_for_data)
  return(ret)
}




mcleod.estimate.CI.based.on.medians = function(x.vec,
                              n.vec,
                              conf.level = 0.95,nr.perm = 200,
                              q_grid = seq(0.1,0.9,0.1),
                              a.max = c(3),
                              comp_param = mcleod.computational.parameters(nr.gibbs = 200,nr.gibbs.burnin = 100),
                              prior_param = mcleod.prior.parameters(prior.type = MCLEOD.PRIOR.TYPE.BETA.HEIRARCHICAL,Beta.Heirarchical.Levels = 6),
                              shift.size.in.log.odds.scale = NULL,
                              nr.cores = detectCores() - 1,
                              verbose = T){
  
  Start.time.overall = Sys.time()
  K = length(n.vec)  
  M=8
  
  #Part I: generate deconv model for the original data
  if(verbose){
    cat(paste0(' - Computing deconvolution estimate for original data.\n\r'))
  }
  
  mcleod_for_data = mcleod(x.smp = x.vec,n.smp = n.vec,a.limits = c(-a.max,a.max),
                           computational_parameters = comp_param,
                           prior_parameters = prior_param)
  
  median_curve_for_data = compute_medians_curve(mcleod_for_data)
  
  a.vec = mcleod_for_data$parameters_list$a.vec
  
  if(is.null(shift.size.in.log.odds.scale)){
    shift.size.in.log.odds.scale = 2*(a.vec[2]-a.vec[1])
  }
  
  
  nr_points_trim_GE = 1
  nr_points_trim_LE = 1
  
  function_for_Worst_Case_PV_computation_parallelized = function(q_ind,ai,is_LE){
    
    NR_OVER_TO_QUIT = nr.perm*(1-conf.level)/2
    
    current_temp_dir = paste0(tempdir(),'/parallel_worker_',q_ind,'_',ai,'_',is_LE,'/')
    if(!dir.exists(current_temp_dir)){
      dir.create(current_temp_dir)
    }else{
      for(u in list.files(current_temp_dir,full.names = T)){
        file.remove(u)
      }  
    }
    parallel_worker = function(worker_id){
      nr_perm_for_worker = ceiling(nr.perm/nr.cores)
      worker_nr_done = 0
      worker_nr_more_extreme = 0
      for(bi in 1:nr_perm_for_worker){
        if(is_LE){
          # assign q_ind prob to location ai - epsilon 
          # and assign (1-q_ind) to location inf
          theta_sample = sample(x = c(a.vec[ai]- shift.size.in.log.odds.scale,a.vec[ai]+M),
                                size = K,
                                replace = T,
                                prob = c(q_grid[q_ind],1-q_grid[q_ind]))
          
        }else{
          # assign q_ind prob to location -inf 
          # and assign (1-q_ind) to location ai + epsilon  
          theta_sample = sample(x = c(a.vec[ai]-M, a.vec[ai]+ shift.size.in.log.odds.scale),
                                size = K,
                                replace = T,
                                prob = c(q_grid[q_ind],1-q_grid[q_ind]))
        }
        
        
        sampled_x_vec = rbinom(K,size = n.vec,prob = inv.log.odds(theta_sample))
        
        
        mcleod_for_resampled_data = mcleod(x.smp = sampled_x_vec,n.smp = n.vec,a.limits = c(-a.max,a.max),
                                           computational_parameters = comp_param,
                                           prior_parameters = prior_param)
        
        median_curve_for_resampled_data = compute_medians_curve(mcleod_for_resampled_data)
        
        worker_nr_done = worker_nr_done +1
        if((is_LE & median_curve_for_resampled_data[ai]<=median_curve_for_data[ai]) |
           (!is_LE & median_curve_for_resampled_data[ai]>=median_curve_for_data[ai]) ){
          worker_nr_more_extreme = worker_nr_more_extreme +1
          ghost_object = list(a=1) # not letting me save a blank list with no name for some reason
          save(ghost_object,file = paste0(current_temp_dir,'signal_file_worker_',
                                    worker_id,'_sig_nr_',worker_nr_more_extreme,'.rdata'))
        }
        if(length(list.files(current_temp_dir,full.names = T)) > NR_OVER_TO_QUIT) #no need to to more permutations, we wont reject
          return(list(worker_nr_more_extreme = worker_nr_more_extreme,worker_nr_done = worker_nr_done))
      }
      
      return(list(worker_nr_more_extreme = worker_nr_more_extreme,worker_nr_done = worker_nr_done))
    }
    
    
    variables_to_pass_as_vector = c('mcleod_for_data',
                                    'median_curve_for_data',
                                    'a.vec',
                                    'K',
                                    'M',
                                    'x.vec',
                                    'n.vec',
                                    'conf.level',
                                    'nr.perm',
                                    'q_grid',
                                    'a.max',
                                    'comp_param',
                                    'prior_param',
                                    'shift.size.in.log.odds.scale',
                                    'nr.cores')
    
    perm_worker_results <- foreach(wi = 1:nr.cores, .options.RNG=1234,.export = variables_to_pass_as_vector) %dorng% {
      library(mcleod);
      return(parallel_worker(wi))
    }
    
    total_nr_more_extreme = 0
    total_nr_computed = 0
    for(wi in 1:length(perm_worker_results)){
      total_nr_more_extreme = total_nr_more_extreme + perm_worker_results[[wi]]$worker_nr_more_extreme
      total_nr_computed = total_nr_computed + perm_worker_results[[wi]]$worker_nr_done
    }
    
    if(dir.exists(current_temp_dir))
      unlist(current_temp_dir)
    
    return((total_nr_more_extreme+1)/(total_nr_computed+1))
  }
  
  function_for_Worst_Case_PV_computation = function(q_ind,ai,is_LE){
    
    NR_OVER_TO_QUIT = nr.perm*(1-conf.level)/2
    nr_more_extreme = 0
    for(bi in 1:nr.perm){
      if(is_LE){
        # assign q_ind prob to location ai - epsilon 
        # and assign (1-q_ind) to location inf
        theta_sample = sample(x = c(a.vec[ai]- shift.size.in.log.odds.scale,a.vec[ai]+M),
                              size = K,
                              replace = T,
                              prob = c(q_grid[q_ind],1-q_grid[q_ind]))
        
      }else{
        # assign q_ind prob to location -inf 
        # and assign (1-q_ind) to location ai + epsilon  
        theta_sample = sample(x = c(a.vec[ai]-M, a.vec[ai]+ shift.size.in.log.odds.scale),
                              size = K,
                              replace = T,
                              prob = c(q_grid[q_ind],1-q_grid[q_ind]))
      }
      
      
      sampled_x_vec = rbinom(K,size = n.vec,prob = inv.log.odds(theta_sample))
      
      
      mcleod_for_resampled_data = mcleod(x.smp = sampled_x_vec,n.smp = n.vec,a.limits = c(-a.max,a.max),
                                         computational_parameters = comp_param,
                                         prior_parameters = prior_param)
      
      median_curve_for_resampled_data = compute_medians_curve(mcleod_for_resampled_data)
      
      if((is_LE & median_curve_for_resampled_data[ai]<=median_curve_for_data[ai]) |
         (!is_LE & median_curve_for_resampled_data[ai]>=median_curve_for_data[ai]) ){
        nr_more_extreme = nr_more_extreme +1
      }
      if(nr_more_extreme > NR_OVER_TO_QUIT) #no need to to more permutations, we wont reject
        return(1)
      
    }
    
    return((nr_more_extreme +1 )/(nr.perm+1))
  }
  
  
  function_for_Pval_LE = function(q_ind,ai){
    #return(function_for_Worst_Case_PV_computation(q_ind, ai, is_LE = T))
    return(function_for_Worst_Case_PV_computation_parallelized(q_ind, ai, is_LE = T))
    
  }
  
  function_for_Pval_GE = function(q_ind,ai){
    #return(function_for_Worst_Case_PV_computation(q_ind, ai, is_LE = F))
    return(function_for_Worst_Case_PV_computation_parallelized(q_ind, ai, is_LE = F))
  }
  
  #function_for_Pval_LE(39,33)
  
  cl <- makeCluster(nr.cores)
  registerDoParallel(cl)
  
  curve_obj = compute.mcleod.CI.curve(mcleod_for_data,
                                      q_grid,
                                      function_for_Pval_LE,
                                      function_for_Pval_GE,
                                      nr_points_trim_LE,
                                      nr_points_trim_GE,
                                      conf.level,verbose = verbose)
  
  stopCluster(cl) 
  
  ret = list()
  class(ret) = CLASS.NAME.MCLEOD.CI.MEDIAN.BASED.STATISTIC
  ret$a.vec = a.vec
  ret$q_grid = q_grid
  #ret$pval_LE = pval_LE
  #ret$pval_GE = pval_GE
  ret$shift.size.in.log.odds.scale = shift.size.in.log.odds.scale
  ret$K = K
  ret$mcleod_for_data = mcleod_for_data
  ret$curve_obj = curve_obj
  ret$Elapsed_Time_Overall = Sys.time() - Start.time.overall
  return(ret)
  
}
