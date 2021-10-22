# Auxiliary Functions and Wrapper:


#' Title
#'
#' @param p 
#'
#' @return
#' @export
#'
#' @examples
log.odds			<- function(p){
  log(p / (1-p))
} 		

#' Title
#'
#' @param theta 
#'
#' @return
#' @export
#'
#' @examples
inv.log.odds		<- function(theta){
  exp(theta) / (1 +exp(theta))
}	

#wrapper for the CPP function
Wrapper_rcpp_Gibbs = function(x.vec,
                              n.vec,
                              a.vec,
                              nr.gibbs,
                              nr.gibbs.burnin,
                              numeric_integration,
                              Verbose,
                              L,
                              Prior_Hyper_Parameters_BetaH_L,
                              Prior_Hyper_Parameters_BetaH_U,
                              Prior_Hyper_Parameters_2LDT,
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
                              integration_step_size = 0.01,
                              Noise_Type = c(0L),
                              Manual_Prior_Given = c(0L),
                              Manual_Prior_Values = c(-4,4),
                              Manual_Prior_Probs = c(1),
                              do_P_k_i_hashing = c(0L),
                              P_k_i_hashing_resolution = 0.0001,
                              offset_vec = rep(0,length(x.vec))
                              ){ 
  
    return(
      rcpp_Gibbs_Prob_Results(x.vec,
                                    n.vec,
                                    a.vec,
                                    nr.gibbs,
                                    nr.gibbs.burnin,
                                    as.integer(numeric_integration),
                                    as.integer(Verbose),
                                    L,
                                    Prior_Hyper_Parameters_BetaH_L,
                                    Prior_Hyper_Parameters_BetaH_U,
                                    Prior_Hyper_Parameters_2LDT,
                                    as.integer(0),
                                    rep(1/(length(a.vec)-1),(length(a.vec)-1)),
                                    FastGammaUsed,
                                    FastGammaBank,
                                    P_k_i_is_given,
                                    P_k_i_precomputed,
                                    integration_step_size, 
                                    PriorType,
                                    I1,
                                    covariates_given,
                                    covariates,
                                    proposal_sd,
                                    beta_prior_sd,
                                    beta_init,
                                    Noise_Type,
                                    Manual_Prior_Given,
                                    Manual_Prior_Values,
                                    Manual_Prior_Probs,
                                    do_P_k_i_hashing,
                                    P_k_i_hashing_resolution,
                                    offset_vec
                              ) 
    )
}



