# Auxiliary Functions and Wrapper:


log.odds			<- function(p) 		log(p / (1-p))
inv.log.odds		<- function(theta)	exp(theta) / (1 +exp(theta))
revert.location = function(i,n){rev(1:n)[i]}

#wrapper for the CPP function
Wrapper_rcpp_Gibbs = function(x.vec,
                              n.vec,
                              a.vec,
                              nr.gibbs,
                              nr.gibbs.burnin,
                              numeric_integration,
                              Verbose,
                              L,
                              FastGammaUsed=0,
                              FastGammaBank=matrix(1,nrow = 1),
                              P_k_i_is_given = 0L,
                              P_k_i_precomputed = matrix(1,nrow = 1),
                              PriorType = 0,
                              I1 = length(a.vec)-1,
                              covariates_given = 0,
                              covariates = matrix(c(1),nrow = 1),
                              proposal_sd = c(1),
                              beta_prior_sd = c(1),
                              beta_init = c(1),
                              integtation_step_size = 0.01){ 
  
    return(
      rcpp_Gibbs_Prob_Results(x.vec,
                                    n.vec,
                                    a.vec,
                                    nr.gibbs,
                                    nr.gibbs.burnin,
                                    as.integer(numeric_integration),
                                    as.integer(Verbose),
                                    L,
                                    as.integer(0),
                                    rep(1/(length(a.vec)-1),(length(a.vec)-1)),
                                    FastGammaUsed,
                                    FastGammaBank,
                                    P_k_i_is_given,
                                    P_k_i_precomputed,
                                    integtation_step_size, #Integration step size
                                    PriorType,
                                    I1,
                                    covariates_given,
                                    covariates,
                                    proposal_sd,
                                    beta_prior_sd,
                                    beta_init) 
    )
}



