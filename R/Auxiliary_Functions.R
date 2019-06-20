# Auxiliary Functions and Wrapper:


log.odds			<- function(p) 		log(p / (1-p))
inv.log.odds		<- function(theta)	exp(theta) / (1 +exp(theta))
revert.location = function(i,n){rev(1:n)[i]}

#wrapper for the CPP function
Wrapper_rcpp_Gibbs = function(x.vec,
                              n.vec,
                              a.vec,
                              CDF,
                              theta,
                              nr.gibbs,
                              nr.gibbs.burnin,
                              numeric_integration,
                              Verbose,
                              L,
                              ReturnComputations,
                              FastGammaUsed=0,
                              FastGammaBank=matrix(1,nrow = 1),
                              P_k_i_is_given = 0L,
                              P_k_i_precomputed = matrix(1,nrow = 1),
                              PriorType = 0,
                              I1 = length(a.vec)-1){ 
  
    return(
      rcpp_Gibbs_Prob_Results_BetaH(x.vec,
                                    n.vec,
                                    a.vec,
                                    CDF,
                                    theta,
                                    nr.gibbs,
                                    nr.gibbs.burnin,
                                    as.integer(numeric_integration),
                                    as.integer(Verbose),
                                    L,
                                    as.integer(ReturnComputations),
                                    as.integer(0),
                                    rep(1/(length(a.vec)-1),(length(a.vec)-1)),
                                    FastGammaUsed,
                                    FastGammaBank,
                                    P_k_i_is_given,
                                    P_k_i_precomputed,
                                    0.01, #Integration step size
                                    PriorType,
                                    I1) 
    )
}



