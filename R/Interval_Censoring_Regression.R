CLASS.NAME.MCLEOD.INTERVAL.CENSORING = 'mcleod.ice.obj'


#' Title
#'
#' @param lu.mat 
#' @param L 
#' @param a.rng 
#' @param n.gbbs 
#' @param n.gbbs.burnin 
#'
#' @return
#' @export
#'
#' @examples
mcleod.interval.censoring.density.estimation = function(lu.mat,L = 8, a.rng = c(0,1100) , n.gbbs = 200, n.gbbs.burnin=100){
  
  I			<- 2^L

  tree.parent.lst		<- as.list(1:L)
  for(l in 1:L)
  {
    if(l == 1)	tree.parent.lst[[l]]	<- as.list(0)
    if(l > 1)	tree.parent.lst[[l]]	<- as.list(rep(1:(2^(l-2)),each = 2))
    names(tree.parent.lst[[l]])	<- paste("Level",l," -- Node",1:(2^(l-1)))
  }
  
  
  tree.child.node.lst		<- as.list(1:L)
  for(l in 1:L)
  {
    tree.child.node.lst[[l]]	<- as.list(1:(2^(l-1)))
    names(tree.child.node.lst[[l]])		<- paste("Level",l," --  Node",1:(2^(l-1)))
    for(node in 1:(2^(l-1)))
    {
      tree.child.node.lst[[l]][[node]]	<- c(2*(node-1)+1,2*(node-1)+2)
    }
  }
  
  tree.a.ind.lst		<- as.list(1:L)
  for(l in 1:L)
  {
    tree.a.ind.lst[[l]]	<- as.list(1:(2^(l-1)))
    names(tree.a.ind.lst[[l]])		<- paste("Level",l," -- Node",1:(2^(l-1)))
    for(node in 1:(2^(l-1)))
    {
      p.ind<- (1:(2^(L-l))) + (node-1)*2^(L-l+1)
      q.ind<- p.ind + 2^(L-l) 
      tree.a.ind.lst[[l]][[node]]	<- list(p.ind = p.ind,q.ind = q.ind)
    }
  }
  
  
  #	This a matrix of indices that is used to efficiently sample distributions 
  
  prod.ind.mat	<- array(dim=c(L,I))
  for(l in 1:L)
  {
    p.ind			<- (2^(l-1)):((2^l)-1)
    q.ind			<- p.ind + 2^L-1
    p.q				<- c(rbind(p.ind,q.ind))
    prod.ind.mat[l,]<- rep(p.q,each = 2^(L-l))
  }
  
  a.vec			<- seq(a.rng[1] ,a.rng[2] ,length = I + 1)
  a.med			<- (a.vec[-1] + a.vec[-(I+1)])/2  
  ind.mat			<- outer(lu.mat[,1],a.med,"<") * outer(lu.mat[,2],a.med,">")
  
  K			<- dim(lu.mat)[1]
  pi.gbbs		<- array(dim=c(n.gbbs,I))
  delta.gbbs	<- array(dim=c(n.gbbs,K))
  
  for(g in 1:n.gbbs)
  {
    
     		
    
    #	Initialize  pi.gbbs
    if(g==1)	pi.gbbs[g,]		<- 1/I
    
    
    #	Gibbs step 1: for k = 1 ... K   sample  theta.gbbs[g,k]] conditionally pi.gbbs
    
    nvec		<- rep(0,I)
    
    for(k in 1:K)
    {	
      delta.gbbs[g,k]		<- sample(1:I,1,prob = pi.gbbs[g,]*ind.mat[k,] / sum(pi.gbbs[g,]*ind.mat[k,]))
      nvec[delta.gbbs[g,k]]	<- nvec[delta.gbbs[g,k]] + 1
    }
    
    
    
    #	Gibbs step 2: for each node independently sample node prob conditionally on n.vec 
    
    if(g < n.gbbs)
    {
      
      pp.vec		<- rep(NA,2^L-1)
      node.ind	<- 1
      
      for(l in 1:L)
        for(node in 1:(2^(l-1)))
        {
          n.p					<- sum(nvec[tree.a.ind.lst[[l]][[node]]$p.ind])
          n.q					<- sum(nvec[tree.a.ind.lst[[l]][[node]]$q.ind])
          pp.vec[node.ind]	<- rbeta(1,n.p + 1,n.q + 1)
          node.ind			<- node.ind + 1	
        }
      
      pq.vec				<- c(pp.vec,1-pp.vec)
      pq.mat				<- array(pq.vec[prod.ind.mat],dim=c(L,I))
      pi.gbbs[g+1,]		<- apply(pq.mat,2,prod)
    }
  }		
  
  
  ret = list(a.vec = a.vec,
             pi.gbbs = pi.gbbs,
             delta.gbbs = delta.gbbs,
             a.med = a.med,
             L = L,
             a.rng = a.rng,
             n.gbbs = n.gbbs,
             n.gbbs.burnin=n.gbbs.burnin,
             lu.mat = lu.mat)
  
  class(ret) = CLASS.NAME.MCLEOD.INTERVAL.CENSORING
  return(ret)
  
}



#' Title
#'
#' @param mcleod.ice.obj 
#'
#' @return
#' @export
#'
#' @examples
plot.posterior.interval.censoring = function(mcleod.ice.obj){
  if(class(mcleod.ice.obj)!= CLASS.NAME.MCLEOD.INTERVAL.CENSORING){
    stop('Input for plot.posterior.interval.censoring(...) must be object returned from mcleod.interval.censoring.density.estimation(...)')
  }
  n.gbbs = mcleod.ice.obj$n.gbbs
  n.gbbs.burnin = mcleod.ice.obj$n.gbbs.burnin
  a.vec = mcleod.ice.obj$a.vec
  pi.gbbs = mcleod.ice.obj$pi.gbbs
  vec.for.limits = mcleod.ice.obj$lu.mat[,2]
  vec.for.limits = vec.for.limits[is.finite(vec.for.limits)]
  limits_censoring = c(min(vec.for.limits), max(vec.for.limits))
  ind_to_plot = which(a.vec >= min(limits_censoring) & a.vec <= max(limits_censoring) )
  
  for(g in (n.gbbs.burnin+1):n.gbbs)
  {
    
    if(g == n.gbbs.burnin + 1) plot(a.vec[ind_to_plot],1 - c(cumsum(pi.gbbs[g,]))[ind_to_plot], type = "l", lwd = 0.2, col = 'gray',xlim = limits_censoring,ylim = c(0,1),xlab = 'Time Units',ylab = 'Survival')
    if(g >  n.gbbs.burnin + 1) lines(a.vec[ind_to_plot],1 - c(cumsum(pi.gbbs[g,]))[ind_to_plot], lwd = 0.2, col = 'gray')
  }
  
  lines(a.vec[ind_to_plot],1 - c(cumsum(apply(pi.gbbs[(n.gbbs.burnin+1):n.gbbs,],2,mean)))[ind_to_plot],col = 3,lwd = 2)
  
}

if(F){
  
  #install.packages('icenReg')
  
  ############################################################################
  
  #	3. Apply hBeta sampler to icenReg package mice data example
  
  ############################################################################
  
  library(icenReg)
  data(miceData)
  
  
  miceData.mat	<- cbind(miceData$l,miceData$u)
  fit.all <- ic_np(miceData.mat)
  
  
  a	<- mcleod.interval.censoring.density.estimation(miceData.mat,a.rng = c(0,1200),n.gbbs = 200)
  
  plot.posterior.interval.censoring(a)
  
  lines(fit.all,col=2,lwd = 2)
  
  ############################################################################
  
  #	4. Compare hBeta to NPMLE survival estimate
  
  ############################################################################
  
  r.surv.time	<- function(n)
  {
    u.min		<- c(0,700,800)
    u.max		<- c(700,800,1200)
    i			<- sample(1:3,n,replace = T,prob = c(0.3,0.5,0.2))
    return(round(u.min[i] + runif(n) * (u.max[i] - u.min[i])))
  }
  
  plot(ecdf(r.surv.time(1500)))
  
  
  
  n <- 20		# 20,  50, 150, 200, 400
  surv.time			<- r.surv.time(n)
  
  surv.interval.mat	<- array(dim=c(n,4))
  for(i in 1:n)		surv.interval.mat[i,]	<- round(sort(c(0,runif(2),1)) * 1200)
  
  rnk.vec	<- apply(cbind(surv.time+0.5,surv.interval.mat[,2:3]),1,rank)[1,]
  l.vec	<- t(surv.interval.mat)[rnk.vec + (0:(n-1))*4]
  u.vec	<- t(surv.interval.mat)[rnk.vec + 1 + (0:(n-1))*4]
  u.vec	<- ifelse(u.vec == 1200,Inf,u.vec)
  
  # cbind(surv.interval.mat,surv.time,rnk.vec,l.vec,u.vec)
  
  a	<- mcleod.interval.censoring.density.estimation(cbind(l.vec,u.vec),a.rng = c(0,1300),n.gbbs = 200)
  
  plot.posterior.interval.censoring(a)
  
  lines(ic_np(cbind(l.vec,u.vec)),col=4,lwd = 2)
  
  lines(sort(r.surv.time(10^4)),seq(1,0,length  = 10^4),col=2,lwd = 2)

}





