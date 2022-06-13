// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_Gibbs_Prob_Results
List rcpp_Gibbs_Prob_Results(NumericVector x_vec, NumericVector n_vec, NumericVector a_vec, IntegerVector n_gibbs, IntegerVector n_gibbs_burnin, IntegerVector IsExact, IntegerVector Verbose, IntegerVector L, NumericMatrix Prior_Hyper_Parameters_BetaH_L, NumericMatrix Prior_Hyper_Parameters_BetaH_U, NumericMatrix Prior_Hyper_Parameters_2LDT, IntegerVector InitGiven, NumericVector Init, IntegerVector Sample_Gamma_From_Bank, NumericMatrix Bank, IntegerVector P_k_i_is_given, NumericMatrix P_k_i_precomputed, NumericVector Pki_Integration_Stepsize, IntegerVector Prior_Type, IntegerVector Two_Layer_Dirichlet_I1, IntegerVector covariates_given, NumericMatrix covariates, NumericVector proposal_sd, NumericVector beta_prior_sd, NumericVector beta_init, IntegerVector Noise_Type, IntegerVector manual_beta_dist_given, NumericVector manual_beta_dist_values, NumericMatrix manual_beta_dist_Probs, IntegerVector do_P_k_i_hashing, NumericVector P_k_i_hashing_resolution, NumericVector offset_vec);
RcppExport SEXP _mcleod_rcpp_Gibbs_Prob_Results(SEXP x_vecSEXP, SEXP n_vecSEXP, SEXP a_vecSEXP, SEXP n_gibbsSEXP, SEXP n_gibbs_burninSEXP, SEXP IsExactSEXP, SEXP VerboseSEXP, SEXP LSEXP, SEXP Prior_Hyper_Parameters_BetaH_LSEXP, SEXP Prior_Hyper_Parameters_BetaH_USEXP, SEXP Prior_Hyper_Parameters_2LDTSEXP, SEXP InitGivenSEXP, SEXP InitSEXP, SEXP Sample_Gamma_From_BankSEXP, SEXP BankSEXP, SEXP P_k_i_is_givenSEXP, SEXP P_k_i_precomputedSEXP, SEXP Pki_Integration_StepsizeSEXP, SEXP Prior_TypeSEXP, SEXP Two_Layer_Dirichlet_I1SEXP, SEXP covariates_givenSEXP, SEXP covariatesSEXP, SEXP proposal_sdSEXP, SEXP beta_prior_sdSEXP, SEXP beta_initSEXP, SEXP Noise_TypeSEXP, SEXP manual_beta_dist_givenSEXP, SEXP manual_beta_dist_valuesSEXP, SEXP manual_beta_dist_ProbsSEXP, SEXP do_P_k_i_hashingSEXP, SEXP P_k_i_hashing_resolutionSEXP, SEXP offset_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x_vec(x_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_vec(n_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a_vec(a_vecSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_gibbs(n_gibbsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_gibbs_burnin(n_gibbs_burninSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type IsExact(IsExactSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Verbose(VerboseSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Prior_Hyper_Parameters_BetaH_L(Prior_Hyper_Parameters_BetaH_LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Prior_Hyper_Parameters_BetaH_U(Prior_Hyper_Parameters_BetaH_USEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Prior_Hyper_Parameters_2LDT(Prior_Hyper_Parameters_2LDTSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type InitGiven(InitGivenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Init(InitSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Sample_Gamma_From_Bank(Sample_Gamma_From_BankSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Bank(BankSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type P_k_i_is_given(P_k_i_is_givenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P_k_i_precomputed(P_k_i_precomputedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Pki_Integration_Stepsize(Pki_Integration_StepsizeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Prior_Type(Prior_TypeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Two_Layer_Dirichlet_I1(Two_Layer_Dirichlet_I1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type covariates_given(covariates_givenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type proposal_sd(proposal_sdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_prior_sd(beta_prior_sdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Noise_Type(Noise_TypeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type manual_beta_dist_given(manual_beta_dist_givenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type manual_beta_dist_values(manual_beta_dist_valuesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type manual_beta_dist_Probs(manual_beta_dist_ProbsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type do_P_k_i_hashing(do_P_k_i_hashingSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P_k_i_hashing_resolution(P_k_i_hashing_resolutionSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset_vec(offset_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_Gibbs_Prob_Results(x_vec, n_vec, a_vec, n_gibbs, n_gibbs_burnin, IsExact, Verbose, L, Prior_Hyper_Parameters_BetaH_L, Prior_Hyper_Parameters_BetaH_U, Prior_Hyper_Parameters_2LDT, InitGiven, Init, Sample_Gamma_From_Bank, Bank, P_k_i_is_given, P_k_i_precomputed, Pki_Integration_Stepsize, Prior_Type, Two_Layer_Dirichlet_I1, covariates_given, covariates, proposal_sd, beta_prior_sd, beta_init, Noise_Type, manual_beta_dist_given, manual_beta_dist_values, manual_beta_dist_Probs, do_P_k_i_hashing, P_k_i_hashing_resolution, offset_vec));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_Generate_Fast_Gamma_Bank
NumericMatrix rcpp_Generate_Fast_Gamma_Bank(IntegerVector Size);
RcppExport SEXP _mcleod_rcpp_Generate_Fast_Gamma_Bank(SEXP SizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Size(SizeSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_Generate_Fast_Gamma_Bank(Size));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_Generate_Gamma_from_Fast_Gamma_Bank
NumericVector rcpp_Generate_Gamma_from_Fast_Gamma_Bank(NumericVector x, NumericMatrix Bank);
RcppExport SEXP _mcleod_rcpp_Generate_Gamma_from_Fast_Gamma_Bank(SEXP xSEXP, SEXP BankSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Bank(BankSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_Generate_Gamma_from_Fast_Gamma_Bank(x, Bank));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_Gibbs_Prob_Results_Multiple
List rcpp_Gibbs_Prob_Results_Multiple(IntegerVector threadpool_size, List x_vec_list, NumericVector n_vec, NumericVector a_vec, IntegerVector n_gibbs, IntegerVector n_gibbs_burnin, IntegerVector IsExact, IntegerVector Verbose, IntegerVector L, NumericMatrix Prior_Hyper_Parameters_BetaH_L, NumericMatrix Prior_Hyper_Parameters_BetaH_U, NumericMatrix Prior_Hyper_Parameters_2LDT, IntegerVector InitGiven, NumericVector Init, IntegerVector Sample_Gamma_From_Bank, NumericMatrix Bank, IntegerVector P_k_i_is_given, List P_k_i_precomputed_list, NumericVector Pki_Integration_Stepsize, IntegerVector Prior_Type, IntegerVector Two_Layer_Dirichlet_I1, IntegerVector covariates_given, NumericMatrix covariates, NumericVector proposal_sd, NumericVector beta_prior_sd, NumericVector beta_init, IntegerVector Noise_Type, IntegerVector manual_beta_dist_given, NumericVector manual_beta_dist_values, NumericMatrix manual_beta_dist_Probs, IntegerVector do_P_k_i_hashing, NumericVector P_k_i_hashing_resolution, NumericVector offset_vec);
RcppExport SEXP _mcleod_rcpp_Gibbs_Prob_Results_Multiple(SEXP threadpool_sizeSEXP, SEXP x_vec_listSEXP, SEXP n_vecSEXP, SEXP a_vecSEXP, SEXP n_gibbsSEXP, SEXP n_gibbs_burninSEXP, SEXP IsExactSEXP, SEXP VerboseSEXP, SEXP LSEXP, SEXP Prior_Hyper_Parameters_BetaH_LSEXP, SEXP Prior_Hyper_Parameters_BetaH_USEXP, SEXP Prior_Hyper_Parameters_2LDTSEXP, SEXP InitGivenSEXP, SEXP InitSEXP, SEXP Sample_Gamma_From_BankSEXP, SEXP BankSEXP, SEXP P_k_i_is_givenSEXP, SEXP P_k_i_precomputed_listSEXP, SEXP Pki_Integration_StepsizeSEXP, SEXP Prior_TypeSEXP, SEXP Two_Layer_Dirichlet_I1SEXP, SEXP covariates_givenSEXP, SEXP covariatesSEXP, SEXP proposal_sdSEXP, SEXP beta_prior_sdSEXP, SEXP beta_initSEXP, SEXP Noise_TypeSEXP, SEXP manual_beta_dist_givenSEXP, SEXP manual_beta_dist_valuesSEXP, SEXP manual_beta_dist_ProbsSEXP, SEXP do_P_k_i_hashingSEXP, SEXP P_k_i_hashing_resolutionSEXP, SEXP offset_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type threadpool_size(threadpool_sizeSEXP);
    Rcpp::traits::input_parameter< List >::type x_vec_list(x_vec_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n_vec(n_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a_vec(a_vecSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_gibbs(n_gibbsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type n_gibbs_burnin(n_gibbs_burninSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type IsExact(IsExactSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Verbose(VerboseSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Prior_Hyper_Parameters_BetaH_L(Prior_Hyper_Parameters_BetaH_LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Prior_Hyper_Parameters_BetaH_U(Prior_Hyper_Parameters_BetaH_USEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Prior_Hyper_Parameters_2LDT(Prior_Hyper_Parameters_2LDTSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type InitGiven(InitGivenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Init(InitSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Sample_Gamma_From_Bank(Sample_Gamma_From_BankSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Bank(BankSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type P_k_i_is_given(P_k_i_is_givenSEXP);
    Rcpp::traits::input_parameter< List >::type P_k_i_precomputed_list(P_k_i_precomputed_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Pki_Integration_Stepsize(Pki_Integration_StepsizeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Prior_Type(Prior_TypeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Two_Layer_Dirichlet_I1(Two_Layer_Dirichlet_I1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type covariates_given(covariates_givenSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type proposal_sd(proposal_sdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_prior_sd(beta_prior_sdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Noise_Type(Noise_TypeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type manual_beta_dist_given(manual_beta_dist_givenSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type manual_beta_dist_values(manual_beta_dist_valuesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type manual_beta_dist_Probs(manual_beta_dist_ProbsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type do_P_k_i_hashing(do_P_k_i_hashingSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type P_k_i_hashing_resolution(P_k_i_hashing_resolutionSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset_vec(offset_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_Gibbs_Prob_Results_Multiple(threadpool_size, x_vec_list, n_vec, a_vec, n_gibbs, n_gibbs_burnin, IsExact, Verbose, L, Prior_Hyper_Parameters_BetaH_L, Prior_Hyper_Parameters_BetaH_U, Prior_Hyper_Parameters_2LDT, InitGiven, Init, Sample_Gamma_From_Bank, Bank, P_k_i_is_given, P_k_i_precomputed_list, Pki_Integration_Stepsize, Prior_Type, Two_Layer_Dirichlet_I1, covariates_given, covariates, proposal_sd, beta_prior_sd, beta_init, Noise_Type, manual_beta_dist_given, manual_beta_dist_values, manual_beta_dist_Probs, do_P_k_i_hashing, P_k_i_hashing_resolution, offset_vec));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mcleod_rcpp_Gibbs_Prob_Results", (DL_FUNC) &_mcleod_rcpp_Gibbs_Prob_Results, 32},
    {"_mcleod_rcpp_Generate_Fast_Gamma_Bank", (DL_FUNC) &_mcleod_rcpp_Generate_Fast_Gamma_Bank, 1},
    {"_mcleod_rcpp_Generate_Gamma_from_Fast_Gamma_Bank", (DL_FUNC) &_mcleod_rcpp_Generate_Gamma_from_Fast_Gamma_Bank, 2},
    {"_mcleod_rcpp_Gibbs_Prob_Results_Multiple", (DL_FUNC) &_mcleod_rcpp_Gibbs_Prob_Results_Multiple, 33},
    {NULL, NULL, 0}
};

RcppExport void R_init_mcleod(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
