// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppNumerical.h>
#include <math.h>
#include <Rmath.h>
#include <ctime>



using namespace Rcpp;
using namespace Numer;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]


/**
 * Class is used to sample faster from a gamma random variable. \
 * This is done by pre-sampling gamma random variables with different rates and recombining them when needed.
 * The bank itself is represented as a matrix, so it could be wrapped and parsed back to R.
 */
class Fast_Gamma_Sampler{
  
  //Constants for the indices of different shapes in the bank
  static const int INDEX_0_1  = 0;
  static const int INDEX_1    = 1;
  static const int INDEX_10   = 2;
  
  //Total number of bank columns
  static const int BANK_COLS  = 3;
  
  /**
   * Bank constructor, samples each of the columns. Receives as parameter the required bank size.
   */
  public:
  NumericMatrix generate_Bank(IntegerVector BankSize){
    int _sample_size = BankSize(0);
    NumericMatrix Bank(_sample_size , BANK_COLS);
    Bank(_,INDEX_0_1)  = rgamma(_sample_size,0.1,1.0);
    Bank(_,INDEX_1)    = rgamma(_sample_size,1,1.0);
    Bank(_,INDEX_10)   = rgamma(_sample_size,10,1.0);
    return(Bank);
  }
    
  /**
   * Generate a required Gamma random variable with shape given by the argument x using the bank.
   */
  double sample_Gamma_from_Bank(double x,NumericMatrix Bank){
    
    double VALUES_OF_COLS[BANK_COLS] = {0.1,1,10}; 
    double _x = x;
    double _cummulator = 0.0;
    int _bank_size = Bank.nrow();
    double _bank_size_d = (double)_bank_size;
    int _sampled_pointer = 0;
    
    // go over the different coloumns, from largest to smallest, and sample values from each col.
    // While sampling, aggregating the hape value for each column to a cummulator. 
    // Stop when a given resolution is acheived.
    
    for(int i=BANK_COLS - 1 ; i>=0 ; i--){
      while(_x >= VALUES_OF_COLS[i]){
        _sampled_pointer = (int)((_bank_size_d)* unif_rand());
        _cummulator += Bank(_sampled_pointer, i);
        _x -= VALUES_OF_COLS[i];
      }
    }
    
    return _cummulator;
  }
  
};
 

// ###################################################################################################
//    Class for Gibbs sampling based estimation of ECDF, for noisy random variables.
//    by Barak Brill, based on an earlier R implementation by Pallavi Basu and Daniel Yekutieli.
//    ---------------------------------------------------------------------------------------------
//    This class is built and called from R. The test statistic is the posterior probability of a beta binomial tree, over a set of terminal nodes.
//    See documentation in R and paper (link in R package) for full details about the method.
// ###################################################################################################

class Gibbs_Sampler{

  // Should exact integration be used?
  bool ExactIntegration;
  
  // Should the code be verbose, and show messages.
  bool Verbose;
  
  // Number of intervals in theta space (partition size of a)
  int I;
  
  // Number of samples
  int K;
  
  // Number of Gibbs iterations
  int n_gibbs;
  
  // Number of Gibbs iterations to be removed, when computing test statistic
  int n_gibbs_burn_in;
  
  // sets the number of levels for the beta tree, irrelevant on dirichlet prior
  int L; 
  
  // marks the type of prior to be used- 0= 2 Layer dirichlet, 1= L level beta tree
  int Prior_Type;
  
  // Parameters I1,I2 for Two Layer Dirichlet
  int TwoLayerDirichlet_I1;
  int TwoLayerDirichlet_I2;
  
  
  // Theta value for which to compute the Gibbs sampler statistic
  NumericVector theta;
  
  // vector of observed counts (over observations)
  NumericVector x_vec;
  
  // vector for number of tries (for each observation)
  NumericVector n_vec;
  
  // vector of interval edges (array of length I + 1) of theta space
  NumericVector a_vec;
  
  // p_k_i of posterior probabilities
  NumericMatrix p_k_i;
  
  // Gibbs sampler state array, for pi (prior distribution for p's)
  NumericMatrix pi_smp;
  
  // Gibbs sampler state array, for n (number of observations in each a interval)
  NumericMatrix n_smp;
  
  // Variables used for setting the initial conditions for the Gibbs sampler
  bool InitGiven;
  NumericVector PredefinedInit;
  
  // variables used for time testing
  clock_t begin;
  clock_t end;
  
  // Is the mechanism for fast-sampling gamma variables used
  bool Fast_Sample_Gamma;
  
  // Gamma Bank for fast sampling, represented as a matrix
  NumericMatrix Fast_Sample_Gamma_Bank;
  
  // Object for fast gamma sampling
  Fast_Gamma_Sampler FGS;
  
  // stepsize for numeric intergration of P_k,i
  double dbinom_integration_stepsize;
  
  bool covariates_given;
  
  NumericMatrix covariates;
  
  int Nr_covariates;
  
  NumericVector proposal_sd;
  
  NumericVector beta_prior_sd;
  
  NumericMatrix beta_smp;
  
  NumericMatrix beta_suggestion;
  
  NumericVector proposal_approved;
  
  NumericVector beta;
  
  NumericMatrix p_k_i_suggestion;
  
  NumericVector a_v_plus_theta;
  
  NumericVector p_k_i_vec_computation_result;
  
  double elapsed_secs;
  
  public:  
    
  /*
   * Constructor for the Gibbs Sampler Class
   */
  Gibbs_Sampler(NumericVector x_vec_p,                       // x vector
                     NumericVector n_vec_p,                  // n vector
                     NumericVector a_vec_p,                  // Interval partition
                     IntegerVector n_gibbs_p,                // total number of gibs iterations
                     IntegerVector n_gibbs_burnin_p,         // number of gibbs iterations discarded - "burn in"
                     IntegerVector IsExact_p,                // should integration be exact - ?
                     IntegerVector Verbose_p,                // should be verbose ?
                     IntegerVector L_,                       // Number of tree levels
                     IntegerVector Init_Given,               // Is a initialization point for the Gibbs sampler given
                     NumericVector Init,                     // The given initialization 
                     IntegerVector Fast_Gamma_Sample_p,      // Is a fast sampling gamma bank given
                     NumericMatrix Fast_Gamma_Sample_Bank_p, // FGS bank as matrix, to be used
                     IntegerVector P_k_i_is_given,           // Is a matrix of P_k_i values been precomputed?
                     NumericMatrix P_k_i_precomputed,        // The precomputed matrix of P_k_i values
                     NumericVector integration_stepsize_p,   // Integration step size
                     IntegerVector Prior_Type_p,             // Type of prior used, 0 = 2 Layer dirichlet, 1 = L level binomial tree.
                     IntegerVector Two_Layer_Dirichlet_I1_p, // Parameter I1 for 2-Layer_Dirichlet
                     IntegerVector covariates_given_p,
                     NumericMatrix covariates_p,
                     NumericVector proposal_sd_p,
                     NumericVector beta_prior_sd_p,
                     NumericVector beta_init
                     ){
    
    GetRNGstate(); // Take the seed
    
    //Should exact integration be used?
    ExactIntegration = (IsExact_p[0] == 1) ? true: false;
    
    // Should the code be verbose (e.g., present times)
    Verbose = (Verbose_p[0] == 1) ? true: false;
    
    InitGiven = (Init_Given[0] == 1) ? true: false;
    PredefinedInit = Init;
    
    // Initialize class variables by input
    K = x_vec_p.length();
    n_gibbs = n_gibbs_p(0);
    n_gibbs_burn_in = n_gibbs_burnin_p(0);
    
    theta = NumericVector(K);
    L = L_(0);
    I = std::pow(2,L); // removed +1
    
    Prior_Type = Prior_Type_p(0);
    
    if(Prior_Type == 0){
      TwoLayerDirichlet_I1 = Two_Layer_Dirichlet_I1_p(0);
      TwoLayerDirichlet_I2 = (a_vec_p.length()-1) / TwoLayerDirichlet_I1;
    }
      
    
    dbinom_integration_stepsize = integration_stepsize_p(0);
    
    //handling of fast gamma sample for BetaH
    Fast_Sample_Gamma =  (Fast_Gamma_Sample_p[0] == 1) ? true: false;
    Fast_Sample_Gamma_Bank = Fast_Gamma_Sample_Bank_p;
    
    // Data and CDF function grid
    x_vec = x_vec_p;
    n_vec = n_vec_p;
    a_vec = a_vec_p;
    a_v_plus_theta = NumericVector(a_vec.length()); // preallocated, used for computing the rows of p_k_i
    covariates_given = (covariates_given_p[0] == 1) ? true: false;
    covariates = covariates_p;
    proposal_sd = proposal_sd_p;
    Nr_covariates = covariates.ncol();
    beta_prior_sd = beta_prior_sd_p;
    beta_smp = NumericMatrix(Nr_covariates,n_gibbs);
    if(covariates_given){
      beta = beta_init;
      compute_theta_vec(beta,theta);
    }else{
      beta = NumericVector(Nr_covariates);
    }
      
    
    beta_suggestion = NumericMatrix(Nr_covariates,n_gibbs);
    proposal_approved = NumericVector(n_gibbs);
    p_k_i_vec_computation_result = NumericVector(I);
    p_k_i_suggestion = NumericMatrix(K,I);  
    begin = clock();  
    
    ///////////////////compute p_k_i
    
    //handling precomputation of P_k_i
    if(P_k_i_is_given[0] == 1){
      p_k_i = NumericMatrix(K,I);  
      p_k_i = P_k_i_precomputed;
    }else{
      p_k_i = NumericMatrix(K,I);  
      compute_p_k_i(x_vec,n_vec,theta,a_vec,p_k_i);  
    }
    
    
    //call the Gibbs sampler:
    if(Prior_Type == 0)
      perform_Gibbs_Sampling_2_Layer_Dirichlet();
    
    if(Prior_Type == 1)
      perform_Gibbs_Sampling_BetaHeirarchical();
    
    end = clock();
    
    //measure times:
    elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    
    if(Verbose)
      Rprintf("\n\r Elapsed time for Gibbs sampler (s) : %lf \n\r",elapsed_secs);
    
    // Put the seed bank in place
    PutRNGstate();
  }
    
  /*
   * Function returns the p_k_i matrix computed for this class
   */
  NumericMatrix get_p_k_i(){return (p_k_i);}
  
  /*
   * Function returns the matrix of sampled observations - by the gibbs sampler
   */
  NumericMatrix get_n_smp(){return(n_smp);}  
  
  /*
   * Function returns the matrix of sampled pi's - by the gibbs sampler
   */  
  NumericMatrix get_pi_smp(){return(pi_smp);}  
    
  /*
   * Function returns the matrix of sampled beta's - by the gibbs sampler
   */  
  NumericMatrix get_beta_smp(){return(beta_smp);}  
    
  NumericMatrix get_beta_suggestion(){return(beta_suggestion);}  
    
  NumericVector get_proposal_approved(){return(proposal_approved);}  
    
  double get_elapsed_secs(){return(elapsed_secs);}
    
    
  //***************************
  // Function for handling covariates
  //***************************
  
  //function for folding theta from covariates
  void compute_theta_vec(NumericVector beta_candidate, NumericVector theta_output){
    double _cummulator = 0;
    int i;
    int j;
    for(i=0;i<theta_output.length();i++){
      _cummulator = 0;
      for(j=0;j<beta_candidate.length();j++){
        _cummulator += beta_candidate(j) * covariates(i,j);
      }
      theta_output(i) = _cummulator;
    }
  }
    
  //function for beta suggestion
  void generate_Beta_suggestion(NumericVector current_beta, NumericVector beta_suggestion_out){
    for(int i=0;i<current_beta.length();i++){
      beta_suggestion_out(i) = current_beta(i) + rnorm(1,0,proposal_sd(i))(0);
    }
  }
  
  
  //function for computing log likelihood from for a solution
  double compute_loglikelihood(NumericMatrix current_P_k_i, NumericVector current_pi, NumericVector current_beta, double _ll_term_added_for_pi){
    double _ret_ll = 0;
    double _current_sample_prob = 0;
    
    //iterate over observations:
    for(int i=0;i<K;i++){
      _current_sample_prob = 0;
      
      for(int j=0;j<current_pi.length();j++){
        _current_sample_prob += current_P_k_i(i,j)*current_pi(j);  
      }
      _ret_ll += log(_current_sample_prob);
    }
    
    
    //handle term added for pi:
    for(int j=0;j<current_beta.length();j++){
      _ret_ll += log(
                    R::dnorm( current_beta(j), 0.0, beta_prior_sd(j),0) //should change to 1 and check - maybe their log is better...
                    );
    }
    
    // handle prior for beta:
    _ret_ll += _ll_term_added_for_pi;
    
    return(_ret_ll);
  }
  
  //***************************
  // Auxilary functions for the Beta Heirarchical case
  //***************************
  inline int left_descendant(int l,int i, int l2){return i * pow(2, l2 - l);} // the index in layer l2, of the left most descendent of the ith node on the lth layer
  inline int right_descendant(int l,int i, int l2){return (i + 1 ) * pow(2, l2 - l) - 1;} // the index in layer l2, of the right most descendent of the ith node on the lth layer
  inline int left_bottom_descendant(int l,int i){return left_descendant(l,i,L);} // the index in layer L, of the left most descendent of the ith node on the lth layer
  inline int right_bottom_descendant(int l,int i){return right_descendant(l,i,L);} // the index in layer L, of the right most descendent of the ith node on the lth layer
  inline int index_of_left_child_next_level(int i){return 2*i;} // the index of the left son, on the l+1 layer, of the ith node on the ith layer
  inline int index_of_right_child_next_level(int i){return (2*i + 1);} // the index of the right son, on the l+1 layer, of the ith node on the ith layer
  
  
  /**
   * Specific verison of cumsum, that adds a value of zero before.
   */
  void CumSum_ZeroPrefix(NumericVector x,NumericVector res){ 
    res(0) = 0.0;
    for(int i=1;i<res.length();i++){
      res(i) = res(i-1) + x(i-1);
    }
  }
  
  
  /**
   * Function for summing over a vector - from starting index, to ending index, inclusive
   */
  double SumOverVector(NumericVector x, int starting_ind, int ending_ind_inclusive){
    double ret =0;
    for(int i = starting_ind ; i <= ending_ind_inclusive  ; i++){
      ret += x[i];
    }
    return(ret);
  }  
  
  /**
   * Inline function for sampling a beta random variable, using the FGS object
   */
  inline double my_rbeta(double p,double q){
    double gamma_p = FGS.sample_Gamma_from_Bank(p,Fast_Sample_Gamma_Bank);
    double gamma_q = FGS.sample_Gamma_from_Bank(q,Fast_Sample_Gamma_Bank);
    return(gamma_p/(gamma_p + gamma_q));
  }
    
  /***********************************************************************
   * Function performs the actual Gibbs sampling - for a Beta Heirarchical tree generating prior
   ***********************************************************************/
  void perform_Gibbs_Sampling_BetaHeirarchical(){
    bool SHOW_DEBUG_MSGS = false;
    
    pi_smp = NumericMatrix(I,n_gibbs); // gibbs samples of the pis generated for each node of the tree.
    n_smp  = NumericMatrix(I,n_gibbs); // gibbs samples of the number of P's sampled, in each iteration, for each node of the tree.
    
    NumericVector _pst_vec(I); // used for the posterior distribuion, for each of the K observations
    NumericVector _pst_vec_cumsum(I+1);// cumulative distribution of the above variable
    
    NumericVector _n_vec_cum_sum(I+1);
    
    //for beta sampling
    double _ll_current_solution;
    double _ll_suggestion;
    NumericVector _beta_candidate_placeholder(beta.length());
    NumericVector _theta_candidate_placeholder(K);
    //NumericMatrix _pki_candidate_placeholder(K,I);
    
    
    // We initialize the state. This is actually done using a dirichlet prior
    NumericVector _init_state = sum_matrix_over_rows(p_k_i);
    _init_state = 10.0/((double) K) * _init_state + NumericVector(I,1.0);
    pi_smp(_,0) = rdirichlet(_init_state);
    
    // If an init location for the sampler is given, we used it here.
    
    if(InitGiven){
      pi_smp(_,0) = PredefinedInit;  
    }
    
    
    // used to iterate over the tree, when sampling a node from a heirarchical beta
    int node_in_l = 0; 
    int node_in_next_l = 0;
    double p_sum;
    double q_sum;
    double temp_beta;
    int _place_l_s, _place_l_e, _place_r_s, _place_r_e;
    double _u; // used for checking MH approval
    
    for(int gibbs_itr = 0 ; gibbs_itr < n_gibbs ; gibbs_itr++){
      
      if(SHOW_DEBUG_MSGS && Verbose)
        Rprintf("Starting Gibbs Iter %d \n\r",gibbs_itr);
      
      if(covariates_given){
        // record the current beta
        for(int q=0;q < beta.length();q++)
          beta_smp(q,gibbs_itr) = beta(q);    
      }
      
      // Gibbs step 1: for k = 1 ... K and l = 1 ... L, sample  delta.k[l] conditionally on delta.k[l], x.vec[k], pi.gbbs
      
      for(int k=0 ; k<K ; k++){
        
        
        if(SHOW_DEBUG_MSGS && Verbose  && k==0)
          Rprintf("Computing tree note for k = %d \n\r",k);
        
        node_in_l = 0;
        node_in_next_l = 0;
        
        _pst_vec = p_k_i(k,_) * pi_smp(_,gibbs_itr);
        
        if(SHOW_DEBUG_MSGS && Verbose && k==0){
          Rprintf("Posterior computed, length: %d \n\r",_pst_vec.length());
          for( int i=0;i< _pst_vec.length();i++)
            Rprintf("%lf ",_pst_vec(i));
          Rprintf("\n\r");
        }
        
        
        CumSum_ZeroPrefix(_pst_vec,_pst_vec_cumsum);
        
        if(SHOW_DEBUG_MSGS && Verbose  && k==0){
          Rprintf("CumSum Posterior computed, length: %d \n\r",_pst_vec_cumsum.length());
          for( int i=0;i< _pst_vec_cumsum.length();i++)
            Rprintf("%lf ",_pst_vec_cumsum(i));
          Rprintf("\n\r");
        }
        
        for(int l = 1 ; l <= L; l++){
          
          
          if(SHOW_DEBUG_MSGS && Verbose  && k==0){
            Rprintf("l: %d , node: %d, ", l, node_in_l);
          }
          
          _place_l_e = right_bottom_descendant( l , index_of_left_child_next_level(node_in_l) ); 
          _place_l_s = left_bottom_descendant(  l ,index_of_left_child_next_level(node_in_l) );
          
          
          if(SHOW_DEBUG_MSGS && Verbose  && k==0){
            Rprintf("ind left child: %d, _place_l_s: %d ,  _place_l_e: %d, ",  index_of_left_child_next_level(node_in_l),  _place_l_s, _place_l_e);
          }
          
          p_sum = _pst_vec_cumsum( _place_l_e  + 1 ) - _pst_vec_cumsum( _place_l_s);
          
          _place_r_e = right_bottom_descendant(  l , index_of_right_child_next_level(node_in_l) );
          _place_r_s = left_bottom_descendant(   l ,index_of_right_child_next_level(node_in_l) );
          
          if(SHOW_DEBUG_MSGS && Verbose  && k==0){
            Rprintf("ind right child: %d, _place_r_s: %d , _place_r_e:%d, ", index_of_right_child_next_level(node_in_l),  _place_r_s, _place_r_e );
          }
          
          q_sum = _pst_vec_cumsum( _place_r_e  + 1) -  _pst_vec_cumsum( _place_r_s);
          
          
          node_in_next_l = 2 * node_in_l;
          if(unif_rand() <= q_sum/(q_sum + p_sum))
            node_in_next_l++;  
          
          if(SHOW_DEBUG_MSGS && Verbose  && k==0){
            Rprintf(" node_in_next_l: %d \n\r", node_in_next_l);
          }
          
          node_in_l = node_in_next_l;
        }
        
        if(SHOW_DEBUG_MSGS && Verbose  && k==0){
          Rprintf(" Selected node: %d \n\r", node_in_l);
        }
        
        n_smp(node_in_l, gibbs_itr)++;
      }
      
      
      
      // Gibbs step 2: for each node independently sample node-pi conditionally on n.vec and then hierarchically construct pi.smp
      
      
      if(gibbs_itr < n_gibbs - 1)
      {
        
        //_n_vec_cum_sum = CumSum_ZeroPrefix(n_smp(_,gibbs_itr) );
        CumSum_ZeroPrefix(n_smp(_,gibbs_itr),_n_vec_cum_sum );
        
        for(int i=0;i<I;i++)
          pi_smp(i,gibbs_itr + 1) = 1.0;
        
        for(int l=1;l<=L;l++){
          for(int i=0;i< pow(2,l-1);i++){ 
            _place_l_e = right_bottom_descendant( l , index_of_left_child_next_level(i) ); 
            _place_l_s = left_bottom_descendant(  l ,index_of_left_child_next_level(i) ); 
            
            p_sum = _n_vec_cum_sum( _place_l_e + 1) - _n_vec_cum_sum( _place_l_s );
            
            _place_r_e = right_bottom_descendant( l , index_of_right_child_next_level(i) ); 
            _place_r_s = left_bottom_descendant(  l , index_of_right_child_next_level(i) ); 
            
            q_sum = _n_vec_cum_sum( _place_r_e + 1) - _n_vec_cum_sum( _place_r_s );
            
            if(Fast_Sample_Gamma){
              temp_beta = my_rbeta(p_sum + 1, q_sum + 1);
            }else{
              temp_beta = rbeta(1,p_sum + 1, q_sum + 1)[0];  
            }
            
            // place the results
            for(int p = _place_l_s ; p <= _place_l_e ; p++)
              pi_smp(p, gibbs_itr + 1) = pi_smp(p, gibbs_itr + 1) * temp_beta;
            
            for(int p = _place_r_s ; p <= _place_r_e ; p++)
              pi_smp(p, gibbs_itr + 1) = pi_smp(p, gibbs_itr + 1) * (1.0 - temp_beta) ;
            
            if((SHOW_DEBUG_MSGS && Verbose)){
              Rprintf("Recompiling pi, p_sum = %lf, q_aum = %lf, beta = %lf \n\r",p_sum, q_sum, temp_beta);
              for(int p = _place_l_s ; p <= _place_r_e ; p++)
                Rprintf( "%lf ", pi_smp(p,gibbs_itr + 1) );
              
              Rprintf("\n\r");
            }
            
          }
        }
        
      } // end of check on final gibbs iteration
      
      if(covariates_given && gibbs_itr < n_gibbs - 1){
        
        //compute likelihood of current solution:
        _ll_current_solution =  compute_loglikelihood(p_k_i, pi_smp(_,gibbs_itr+1), beta, 0.0); // it is assumed we dont need to compute the ll term for the prior, because we will be using the MH rule
        
        // generate a new proposal
        generate_Beta_suggestion(beta, _beta_candidate_placeholder);
        beta_suggestion(_,gibbs_itr) = _beta_candidate_placeholder;
        
        // generate theta for the new beta candidate:
        compute_theta_vec(_beta_candidate_placeholder, _theta_candidate_placeholder);
        
        //generate Pki for solution:
        compute_p_k_i(x_vec,n_vec,_theta_candidate_placeholder,a_vec, p_k_i_suggestion);  
        
        // compute likelihood of candidate
        _ll_suggestion = compute_loglikelihood(p_k_i_suggestion, pi_smp(_,gibbs_itr+1), _beta_candidate_placeholder, 0.0); // it is assumed we dont need to compute the ll term for the prior, because we will be using the MH rule
        
        // check for approval
        _u = Rf_runif(0, 1);
        // if approved
        if(_u<= exp(_ll_suggestion - _ll_current_solution)){
          // - replace beta
          std::copy( _beta_candidate_placeholder.begin(), _beta_candidate_placeholder.end(), beta.begin() ) ;
          
          // - replace theta
          std::copy( _theta_candidate_placeholder.begin(), _theta_candidate_placeholder.end(), theta.begin() ) ;
          
          // - replace pki - will be used in the next iteration
          p_k_i = p_k_i_suggestion;
            
          // record approval
          proposal_approved(gibbs_itr) = 1;
          
        }else{
          
          // record disapprobal 
          proposal_approved(gibbs_itr) = 0;
          
        }
        
      }
      
    } // end of for loop on gibbs iter
    
    return;
  }
    
    
  inline int two_layer_dirichlet_left_descendant_by_I1_index(int i1){return ((i1-1)*(TwoLayerDirichlet_I2));} 
  inline int two_layer_dirichlet_right_descendant_by_I1_index(int i1){return ((i1)*(TwoLayerDirichlet_I2) - 1);} 
    
  /***********************************************************************
   * Function performs the actual Gibbs sampling - for a 2-Layer dirichlet tree generating prior
   ***********************************************************************/
  void perform_Gibbs_Sampling_2_Layer_Dirichlet(){
    bool SHOW_DEBUG_MSGS = false;
    
    if(SHOW_DEBUG_MSGS){
      Rprintf("Entered  perform_Gibbs_Sampling_2_Layer_Dirichlet \n\r");
      Rprintf("I1 %d  I2 %d\n\r",TwoLayerDirichlet_I1,TwoLayerDirichlet_I2);  
    }
    
    
    
    pi_smp = NumericMatrix(I,n_gibbs); // gibbs samples of the pis generated for each node of the tree.
    n_smp  = NumericMatrix(I,n_gibbs); // gibbs samples of the number of P's sampled, in each iteration, for each node of the tree.
    
    NumericVector _pst_vec(I); // used for the posterior distribuion, for each of the K observations
    NumericVector _pst_vec_cumsum(I+1);// cumulative distribution of the above variable
    
    NumericVector _n_vec_cum_sum(I+1);
    
    // We initialize the state. This is actually done using a dirichlet prior
    NumericVector _init_state = sum_matrix_over_rows(p_k_i);
    _init_state = 10.0/((double) K) * _init_state + NumericVector(I,1.0);
    pi_smp(_,0) = rdirichlet(_init_state);
    
    // If an init location for the sampler is given, we used it here.
    if(InitGiven){
      pi_smp(_,0) = PredefinedInit;  
    }
    
    // used to iterate over the tree, when sampling a node from a heirarchical beta
    NumericVector _sum_of_pst_over_first_layer(TwoLayerDirichlet_I1);
    int _sampled_i1;
    int _sampled_segment;
    NumericVector _pst_sub_vector(TwoLayerDirichlet_I2);
    NumericVector _dirichlet_parameters_first_layer(TwoLayerDirichlet_I1);
    NumericVector _dirichlet_parameters_second_layer(TwoLayerDirichlet_I2);
    NumericVector _dirichlet_sampled_first_layer(TwoLayerDirichlet_I1);
    NumericVector _dirichlet_sampled_second_layer(TwoLayerDirichlet_I2);
    
    for(int gibbs_itr = 0 ; gibbs_itr < n_gibbs ; gibbs_itr++){
      
      if(SHOW_DEBUG_MSGS && Verbose )
        Rprintf("Starting Gibbs Iter %d \n\r",gibbs_itr);
      
      // Gibbs step 1: for k = 1 ... K and l = 1 ... L, sample  delta.k[l] conditionally on delta.k[l], x.vec[k], pi.gbbs
      
      for(int k=0 ; k<K ; k++){
        
        
        if(SHOW_DEBUG_MSGS && Verbose  && k==0)
          Rprintf("Computing gibbs sample for k = %d \n\r",k);
        
        
        _pst_vec = p_k_i(k,_) * pi_smp(_,gibbs_itr);
        
        CumSum_ZeroPrefix(_pst_vec,_pst_vec_cumsum);
        
        
        // compute sums of posteriors on I1 dirichlet items in the top level
        for(int i1=1;i1<=TwoLayerDirichlet_I1;i1++){
          if(false && gibbs_itr ==0 && k==0)
            Rprintf("DEBUG i1 %d access: [%d,%d] ",i1, two_layer_dirichlet_left_descendant_by_I1_index(i1), two_layer_dirichlet_right_descendant_by_I1_index(i1) + 1);
          
          _sum_of_pst_over_first_layer[i1 - 1] =  _pst_vec_cumsum( two_layer_dirichlet_right_descendant_by_I1_index(i1) + 1 ) -
                                          _pst_vec_cumsum( two_layer_dirichlet_left_descendant_by_I1_index(i1));
          if(false && gibbs_itr ==0 && k==0)
            Rprintf(" value =  %lf \n\r",_sum_of_pst_over_first_layer[i1 - 1]);
        }
        
        _sum_of_pst_over_first_layer = _sum_of_pst_over_first_layer/ sum(_sum_of_pst_over_first_layer);
        
        // sample an an index i1
        _sampled_i1 = sample_ind_with_weights(_sum_of_pst_over_first_layer); // note that this number is zero based
          
        // get posteriror for the relevant sub-vector on level 2, with I2 items.
        for(int i2=1;i2<=TwoLayerDirichlet_I2;i2++){
          _pst_sub_vector[i2-1] = _pst_vec( (_sampled_i1) * TwoLayerDirichlet_I2 + i2 - 1);  
        }
        _pst_sub_vector = _pst_sub_vector/sum(_pst_sub_vector);
          
        // sample an index in the subvect, find the relevant segment in the total vector:.
        _sampled_segment = _sampled_i1 * TwoLayerDirichlet_I2 + sample_ind_with_weights(_pst_sub_vector);
        
        // update N_smp
        
        n_smp(_sampled_segment, gibbs_itr) = n_smp(_sampled_segment, gibbs_itr) + 1.0;
        
      }
      
      // Gibbs step 2: for each node independently sample node-pi conditionally on n.vec and then hierarchically construct pi.smp
     
      if(gibbs_itr < n_gibbs - 1)
      {
        
        CumSum_ZeroPrefix(n_smp(_,gibbs_itr),_n_vec_cum_sum);
        
        for(int i=0;i<I;i++)
          pi_smp(i,gibbs_itr + 1) = 1.0;
        
        //sample dirichlet vector for level 1, size is I1, based on the sums of n_smp in the relevant subvectors.
        //    compute sums:
        for(int i1=1;i1<=TwoLayerDirichlet_I1;i1++){
          _dirichlet_parameters_first_layer(i1-1) = _n_vec_cum_sum( two_layer_dirichlet_right_descendant_by_I1_index(i1) + 1) -
                                                  _n_vec_cum_sum( two_layer_dirichlet_left_descendant_by_I1_index(i1)) +
                                                  1.0;
        }

        //    sample:
        rdirichlet_no_generation(_dirichlet_parameters_first_layer, _dirichlet_sampled_first_layer);
        
        //iterate over I1 subvectors
        for(int i1=1;i1<=TwoLayerDirichlet_I1;i1++){
          //for each subvector of size i1, sample a subvector based on the entries of n_smp in the relevant subvector.
          //   collect parameters:
          for(int i2 = 1;i2<= TwoLayerDirichlet_I2;i2++){
            _dirichlet_parameters_second_layer(i2 - 1) = n_smp(two_layer_dirichlet_left_descendant_by_I1_index(i1) + i2 - 1 , gibbs_itr) + 1.0;
          }
          
          //    sample:
          rdirichlet_no_generation(_dirichlet_parameters_second_layer,_dirichlet_sampled_second_layer);
          
          //multiply by the relevant entry in I1, and place the result in pi_smp
          for(int i2 = 1;i2<= TwoLayerDirichlet_I2;i2++){
            pi_smp(two_layer_dirichlet_left_descendant_by_I1_index(i1) + i2 - 1, gibbs_itr + 1) = _dirichlet_sampled_first_layer(i1-1) * _dirichlet_sampled_second_layer(i2-1);
          }
        }//end of iterate over I1 subvectors
      
        
      } // end of check on final gibbs iteration
    } // end of for loop on gibbs iter
    return;
  }
  
  
  /*
  * Testing Function - Used For Debugging, mainly for the auxilary functions and parameters
  */
  void DebugPrint(){
    if(ExactIntegration)
      Rprintf("Performing Exact Integration \n\r");
    else
      Rprintf("Performing Appoximate Integration \n\r");
    
    Rprintf("log odds 0.5: \n\r");
    Rprintf("%lf \n\r",log_odds(0.5));
    Rprintf("checking inv log odds:\n\r");
    Rprintf("%lf \n\r",inv_log_odds(log_odds(0.5)));
    
    Rprintf("checking dbinom:\n\r");
    Rprintf("%lf \n\r",integrand(log_odds(0.5),1,1));
    
    NumericVector a = NumericVector::create(1.0,1.0,1.0);
    NumericVector n = NumericVector::create(1.0,1.0,1.0);
    
    Rprintf("checking log_d_Dir_Mult:\n\r");
    Rprintf("%lf \n\r",log_d_Dir_Mult(n,a)); //  to be checked against R...
    
  }
    
  /*
   * compute.p.k.i.vec - computes the posterior probability of a measurement to 
   * be received from a value of p, coming from each of the segments of a_v, assuming uniform distribution in each section.
   * the value of ExactIntegration for the class sets whether exact numeric integration over dbinom is done,
   * instead of a normal approximation.
   * 
   * results returned in object p_k_i_vec_computation_result
   * 
   */
  void compute_p_k_i_vec(double x, double n, NumericVector a_v){
    
    double p_hat	= ( (double)(x + (double)0.5) ) / ( (double)( n + (double)1.0 ));
    
    double theta_hat = log_odds(p_hat);
    
    double theta_se	= sqrt(((double)1.0) / ((double)(n * p_hat * (1.0 - p_hat) )) );
    
    double _temp_p_nrm_ul, _temp_p_nrm_ll;
    double _temp;
    double _integral = 0.0;
    double _integral_h = dbinom_integration_stepsize;
    double _integral_p = 0.0;
    double _integral_sum = 0.0;
    
    double _this_density = -1;
    double _last_density = -1;
    for(int i=0;i < I; i++){
      //Check if Exact Numerical Integration is needed
      if(ExactIntegration &&
         -4.0 <= ((a_v(i) - theta_hat) / theta_se) &&
         ((a_v(i+1) - theta_hat) / theta_se) <= 4.0){
         
         //do integration
         _integral = 0.0;
         _integral_p = a_v(i);
         while(_integral_p +_integral_h <= a_v(i+1)){
           if(_last_density == -1){ // this rule is activated only on the first iteration, since we don't have a density computed from the last step
             _last_density = Rf_dbinom(x,n,inv_log_odds(_integral_p),0);
           }
           _this_density = Rf_dbinom(x,n,inv_log_odds(_integral_p +_integral_h),0);
           _integral   +=   (_last_density + _this_density) * _integral_h / 2.0;
           //_integral   +=   Rf_dbinom(x,n,inv_log_odds(_integral_p),0) * _integral_h;
           _integral_p += _integral_h;
           _last_density = _this_density;
         }
         _temp = _integral;
      }else{
        _temp_p_nrm_ul = Rf_pnorm5(a_v(i+1) , theta_hat, theta_se, 1, 0); 
        _temp_p_nrm_ll = Rf_pnorm5(a_v(i)   , theta_hat, theta_se, 1, 0); 
        _temp = (_temp_p_nrm_ul - _temp_p_nrm_ll) / ( (a_v(i+1)) - a_v(i) );
      }
      
      p_k_i_vec_computation_result(i) = _temp;
      _integral_sum += _temp;
    }
    
    
    //normalize:
    for(int i=0; i< I;i++){
      p_k_i_vec_computation_result(i) = p_k_i_vec_computation_result(i)/_integral_sum;
    }
  }
    
  
  /*
   * compute P_k_i matrix, row by row (row = observation).
   * P_k_i gives the posterior probability of the kth obervation to be received from a p coming from the ith interval of the a grid
   * 
   * results returned in p_mat - must be of dimensions K X I
   */
  void compute_p_k_i(NumericVector x_v, NumericVector n_v,NumericVector theta, NumericVector a_v,NumericMatrix p_mat){
    
    //NumericVector a_v_plus_theta(a_v.length());
    for(int k = 0; k < K ; k++)
    {
      for(int j=0;j<a_v_plus_theta.length();j++){
        a_v_plus_theta(j) = a_v(j) +theta(k);
      }
      
      compute_p_k_i_vec( x_v(k), n_v(k), a_v_plus_theta );
      
      for(int j=0;j<I;j++){ // need to benchmark against memcopy - however this requires a row object - not sure how it will work...
        p_mat(k,j) =  p_k_i_vec_computation_result(j);
      }
    }
  }
    
  /*
   * ################################
   * Auxilary Functions:
   * ##############################
   */
  
  
  /*
   * Function for log odds, given by the original R definition of: log.odds			<- function(p) 		log(p / (1-p))
   */
  inline double log_odds(double p)
  {
    return (log(p/(1.0-p)));
  };
  
  
  /*
   * Function for inverse log odds , originally given by inv.log.odds		<- function(theta)	exp(theta) / (1 +exp(theta))
   */  
  inline double inv_log_odds(double theta)
  {
    double _exp_theta = exp(theta);
    return (_exp_theta/(1.0 + _exp_theta));
  };
  
  /*
   * Function for integrand for binomial distribution
   * Given Originally by : integrand 			<-  function(theta,x=2,n=10)		dbinom(x,n,inv.log.odds(theta))
   */
  inline double integrand(double theta,int x, int n)
  {
    return (Rf_dbinom((double)x,(double)n,inv_log_odds(theta),0));
  };
  
  
  /*
   * Function for computing the logarithm of the Dirichlet Multinomial distribution
   * Given Originally by:
   * log.d.Dir.Mult			<- function(n.vec,alpha.vec)	
   * {	
   *  (lgamma(sum(n.vec)+1) + lgamma(sum(alpha.vec)) - lgamma(sum(n.vec+alpha.vec))) + 
   *  (sum( lgamma(n.vec + alpha.vec) - lgamma(n.vec+1) - lgamma(alpha.vec)))
   * }
   */
  
  double log_d_Dir_Mult(NumericVector n_vec, NumericVector alpha_vec){
    double _sum_n_vec = sum(n_vec);
    double _sum_alpha_vec = sum(alpha_vec);
    double _sum_n_vec_plus_alpha_vec = sum(alpha_vec +n_vec);
    
    double _term1 = lgamma(_sum_n_vec + 1.0) + lgamma(_sum_alpha_vec) - lgamma(_sum_n_vec_plus_alpha_vec);
    double _term2 = sum(lgamma(n_vec + alpha_vec) - lgamma(n_vec+1.0) - lgamma(alpha_vec));
    
    return(_term1 + _term2);
  }

  
  /*
   * Function used for producing a single sample, of a dirichlet distribution, with alpha vector x
   */
  NumericVector rdirichlet(NumericVector x){
    NumericVector _ret = NumericVector(x.length());
    rdirichlet_no_generation(x,_ret);
    return(_ret);
  }
  
  /*
   * Function used for producing a single sample, of a dirichlet distribution, with alpha vector x, using preallocated space 
   */
  void rdirichlet_no_generation(NumericVector x,NumericVector ret){
    for(int i=0;i<x.length();i++){
      ret(i) = (rgamma(1,x(i),1))(0);
    } 
    double _sum_ret = sum(ret);
    for(int i=0;i<x.length();i++){
      ret(i) = ret(i)/_sum_ret;
    } 
  }
  
  /**
   * Function returns the summation over rows, for each of the columns of a matrix
   */
  NumericVector sum_matrix_over_rows(NumericMatrix x){
    NumericVector ret(x.ncol());
    for(int i=0;i<x.ncol();i++)
      ret(i) = sum(x(_,i));
    return(ret);
  }
  
  /**
   * Function returns the summation over rows, for each of the columns of a matrix
   */
  NumericVector sum_matrix_over_cols(NumericMatrix x){
    NumericVector ret(x.nrow());
    for(int i=0;i<x.nrow();i++)
      ret(i) = sum(x(i,_));
    return(ret);
  }
  
  /**
   * Function returns the index of the interval, zero based, to be sampled with weights w.
   * It assumes sum(w) == 1.
   */
  int sample_ind_with_weights(NumericVector w){
    
    // time it over 100 - the elapsed time inside this function is zero when timed with the clock
    
    double u = Rf_runif(0, 1);
    double _cummulator = 0;
    for(int i=0;i<w.length();i++){
      _cummulator += w(i);
      if(u <= _cummulator)
        return(i);
    }
    
    return(w.length() - 1); // we should not be reaching this phrase in any case, need to check that.
  }
  
  
};

// #################################################################
// The Next section gives the entry points to the CPP code via RCPP:
// #################################################################

// [[Rcpp::export]]
List rcpp_Gibbs_Prob_Results(NumericVector x_vec,
                                   NumericVector n_vec,
                                   NumericVector a_vec,
                                   IntegerVector n_gibbs,
                                   IntegerVector n_gibbs_burnin,
                                   IntegerVector IsExact,
                                   IntegerVector Verbose,
                                   IntegerVector L,
                                   IntegerVector InitGiven,
                                   NumericVector Init,
                                   IntegerVector Sample_Gamma_From_Bank,
                                   NumericMatrix Bank,
                                   IntegerVector P_k_i_is_given,
                                   NumericMatrix P_k_i_precomputed,
                                   NumericVector Pki_Integration_Stepsize,
                                   IntegerVector Prior_Type,
                                   IntegerVector Two_Layer_Dirichlet_I1,
                                   IntegerVector covariates_given,
                                   NumericMatrix covariates,
                                   NumericVector proposal_sd,
                                   NumericVector beta_prior_sd,
                                   NumericVector beta_init){
  
  Gibbs_Sampler _gibbs(x_vec, n_vec, a_vec, n_gibbs, n_gibbs_burnin, IsExact, Verbose, L, InitGiven, Init, Sample_Gamma_From_Bank, Bank, P_k_i_is_given, P_k_i_precomputed,Pki_Integration_Stepsize,Prior_Type, Two_Layer_Dirichlet_I1,covariates_given,covariates,proposal_sd,beta_prior_sd,beta_init);
  
  List ret;
  ret["p_k_i"]             = _gibbs.get_p_k_i();
  ret["n_smp"]             = _gibbs.get_n_smp();
  ret["pi_smp"]            = _gibbs.get_pi_smp();  
  ret["beta_smp"]          = _gibbs.get_beta_smp();
  ret["beta_suggestion"]   = _gibbs.get_beta_suggestion();
  ret["proposal_approved"] = _gibbs.get_proposal_approved();
  ret["elapsed_secs"]      = _gibbs.get_elapsed_secs();
  return(ret);
}


// [[Rcpp::export]]
NumericMatrix rcpp_Generate_Fast_Gamma_Bank(IntegerVector Size){
  Fast_Gamma_Sampler _fgs;
  NumericMatrix _res =_fgs.generate_Bank(Size);
  return(_res);
}


// [[Rcpp::export]]
NumericVector rcpp_Generate_Gamma_from_Fast_Gamma_Bank(NumericVector x,NumericMatrix Bank){
  Fast_Gamma_Sampler _fgs;
  NumericVector _res(x.length());
  for(int i=0;i<x.length();i++){
    _res(i) = _fgs.sample_Gamma_from_Bank(x(i),Bank);
  }
  return(_res);
}

// [[Rcpp::export]]
List rcpp_hello_world() {

    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;

    return z ;
}


//this section is an example for the quadrature based integration - need to incorporate into the code...
// P(0.3 < X < 0.8), X ~ Beta(a, b)
class BetaPDF: public Func
{
private:
  double a;
  double b;
public:
  BetaPDF(double a_, double b_) : a(a_), b(b_) {}
  
  double operator()(const double& x) const
  {
    return R::dbeta(x, a, b, 0);
  }
};

// [[Rcpp::export]]
Rcpp::List integrate_test()
{
  const double a = 3, b = 10;
  const double lower = 0.3, upper = 0.8;
  const double true_val = R::pbeta(upper, a, b, 1, 0) -
    R::pbeta(lower, a, b, 1, 0);
  
  BetaPDF f(a, b);
  double err_est;
  int err_code;
  const double res = integrate(f, lower, upper, err_est, err_code);
  return Rcpp::List::create(
    Rcpp::Named("true") = true_val,
    Rcpp::Named("approximate") = res,
    Rcpp::Named("error_estimate") = err_est,
    Rcpp::Named("error_code") = err_code
  );
}
