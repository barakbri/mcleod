Estimate_Latent_With_Covarites	<- function(
                                       x.smp,
                                       n.smp,
                                       theta.vec = seq(-4,4,1),
                                       a.max = 4,
                                       L = 4,
                                       I1 = 4,
                                       nr.gibbs = 200,
                                       nr.gibbs.burnin = min(nr.gibbs / 2 , 100),
                                       Fast.Gamma.Used = F,
                                       Fast.Gamma.Bank.Size = 1000L,
                                       VERBOSE = F,
                                       Prior_Type = 0,
                                       I_specificy_parameter = NULL
                                       )
{
  exact.numeric.integration = TRUE # We force exact numeric integration over dbinom for computation of P_k_i. Normal approximation is not sufficiant.
  NR.DIGITS.FOR.LABELS = 4 # This sets the resolution of plot at the moment. setting this to 4 means you cannot go over L=7 (L=7 is maximum)
  COMPUTE_FOR_ANOTHER_PARAMETER = 3 #number of extra grid points to compute over in each column, after Confidence band is reached, for plotting smooth splines (must be at least 2?)
  
  #	Function settings
  a.max = abs(a.max)
  a.min = -1* (a.max)
  K			<- length(x.smp)
  
  theta.vec = sort(unique(theta.vec,-1*theta.vec)) 
  N.theta		<- length(theta.vec)
  Fast.Gamma.Bank = matrix(1,nrow = 1)
  Fast.Gamma.Used.p = 0
  if(Fast.Gamma.Used){
    Fast.Gamma.Used.p = 1
    Fast.Gamma.Bank = rcpp_Generate_Fast_Gamma_Bank(Fast.Gamma.Bank.Size)
  }
  
  I = 2^(L)
  if(Prior_Type == 0){
    if(!is.null(I_specificy_parameter))
        I = I_specificy_parameter
    if(I/I1 != as.integer(I/I1)){
     stop('I/I1 not an integer!') 
    }
  }
  
  a.vec.used		<- seq(a.min,a.max,length = I+1)
  
  res = Wrapper_rcpp_Gibbs(x.smp,
                           n.smp,
                           a.vec.used,
                           nr.gibbs,nr.gibbs.burnin,
                           as.integer(exact.numeric.integration),
                           as.integer(0), #verbose - turned off
                           L,
                           #1, #should computations be returned - this is turned on, since we need the P_k_i's
                           Fast.Gamma.Used.p,
                           Fast.Gamma.Bank,
                           PriorType = Prior_Type,
                           I1 = I1)  
  

  
  ret = list()
  class(ret) = 'NPCI.Object'
  
  ret$parameters_list = list(
    theta.vec = theta.vec,
    a.vec = a.vec.used,
    nr.gibbs = nr.gibbs,
    nr.gibbs.burnin = nr.gibbs.burnin,
    L = L,
    Fast.Gamma.Used = Fast.Gamma.Used,
    Fast.Gamma.Bank.Size = Fast.Gamma.Bank.Size
  )
  
  ret$additional = list(original_stat_res = res)
  return(ret)
}


plot.posterior	<- function(NPCI.obj)
{
 
  library(ggplot2)
  
  burnin = NPCI.obj$parameters_list$nr.gibbs.burnin
  gibbs_pi_plot = cbind(0,t(apply(t(NPCI.obj$additional$original_stat_res$pi_smp),1,cumsum)))
  gibbs_pi_plot = gibbs_pi_plot[-c(1:burnin),]
  means_vec = apply(gibbs_pi_plot,2,mean)
  median_vec = apply(gibbs_pi_plot,2,median)
  
  a.vec = NPCI.obj$parameters_list$a.vec
  
  dt_mean_line = data.frame(a.vec = a.vec,means_vec = means_vec,median_vec = median_vec)
  
  gibbs_cloud_size = dim(gibbs_pi_plot)[1] * dim(gibbs_pi_plot)[2]
  dt_gibbs_cloud = data.frame(a.point = rep(NA,gibbs_cloud_size), CDF.value = rep(NA,gibbs_cloud_size))
  pointer = 1
  for(i in 1:(dim(gibbs_pi_plot)[1])){
    for(j in 1:(dim(gibbs_pi_plot)[2])){
      dt_gibbs_cloud[pointer,] = c(a.vec[j],gibbs_pi_plot[i,j]); pointer = pointer + 1
    }
  }
  gg_obj = ggplot(dt_gibbs_cloud)  + ylim(c(0,1)) +
    geom_line(aes(x = a.vec,y = means_vec),colour = 'red',data = dt_mean_line) +
    geom_point(aes(x = a.point,y = CDF.value),alpha = 0.25,colour = 'gray',data = dt_gibbs_cloud,size = 0.8, shape = 18)
 return(gg_obj)

}



if(F){
  N = 200
  K = 200
  set.seed(1)
  x = rbinom(K,size = N,prob = NPCI:::inv.log.odds(rnorm(K,0,sd = 1)))
  n = rep(N,K)
  
  res = Estimate_Latent_With_Covarites(x, n, L = 6,
                                                     I1=8,
                                                     VERBOSE = T,
                                                     a.max = 4,
                                                     theta.vec = seq(-2,2,0.125), Prior_Type = 1
                                                     )
  plot.posterior(res)
  res$additional$original_stat_res$beta_smp
  
}

