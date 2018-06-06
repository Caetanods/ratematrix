#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// This will have a function to run the MCMC and all the cpp functions that this depends on. The idea is that R functions will prepare the prior parameters and the objects needed to run this function. The output of the function will be the files mcmc and log.

// #######################################################
// ###### Support Functions
// #######################################################

// #######################################################
// ###### The logLikelihood functions
// #######################################################

double logLikNode_C(arma::vec ss, arma::mat sigma_len, arma::mat sigma_len_inv, int k) {
  double val;
  double signal;
  arma::log_det(val, signal, sigma_len); // val is logdeterminant and signal should be 0.
  return -0.5 * ( k * log(2 * arma::datum::pi) + val + as_scalar(trans(ss) * sigma_len_inv * ss));
}

// [[Rcpp::export]]
double logLikPrunningMCMC_C(arma::mat X, int k, int p, arma::uvec nodes, arma::uvec des, arma::uvec anc, arma::uvec names_anc, arma::mat mapped_edge, arma::cube R, arma::vec mu) {
  // X is the data matrix with information for the tips.
  // k is the number of traits in the data.
  // p is the number of regimes in the data.
  // nodes is a vector with the indexes for the nodes. STARTING IN 0!!!
  // des is a vector with the descendents of the nodes. STARTING IN 0!!!
  // anc is a vector with the ancestrals of the nodes. STARTING IN 0!!!
  // mapped_edge is a matrix with the stochastic map.
  // R is a cube with the evolutionary rate matrices.
  // mu is a vector with the root values.

  // All vector of indexes are straight from R. So need to take care when calling values from containers, since R indexing starts from 1.

  // Initiate all the objects.
  arma::uword des_node0;
  arma::uword des_node1;
  arma::vec ss = vec(k); // Number of traits.
  arma::cube Rs1 = cube(R);
  arma::cube Rs2 = cube(R);
  arma::mat Rinv = mat(k, k);
  arma::uvec des_node_vec = uvec(2);
  arma::uword nd;
  arma::uword nd_id;
  arma::uword key_id;
  arma::uword tip;
  arma::uword tip_id;
  arma::uword key_id0;
  arma::uword key_id1;
  int n_nodes = nodes.n_elem;
  arma::mat X0 = mat(k, n_nodes + 1);
  arma::cube V0 = cube(k, k, n_nodes + 1);
  arma::uvec key = uvec(n_nodes); // Not needed at the ROOT.
  arma::uword type;
  // Need to initialize ll as a 0 value.
  double ll = 0.0;
  // node_id are the nodes which ancestral is node 'i'.
  // This will **always** have length 2 assuming the tree is a bifurcating tree.
  arma::uvec node_id = uvec(2);

  // Note: The vector 'anc' is originally a named vector and the names are important.
  //       Here I am assuming a new vector 'names_anc' which is a integer vector (arma::uvec) with the names of the 'anc' cohersed into numeric.
  // 'names_anc' are important to find which type of contrast is done: node & node, tip & tip, or tip & node. Each has a different form of computing the quantities.

  // Loop to traverse the tree.
  // Will visit all the internal nodes including the ROOT.
  for(int i=0; i < n_nodes; i++) {
    
    // The index for the 'des', 'anc', and 'mapped_edge (lines)'.
    node_id = find( anc == nodes[i] );
    type = names_anc[node_id[0]];
    
    // Next is the contrast calculations.This will be different for each 'type' of ancestral node.
    // In this structure of if...else just a single clause will be evaluated.
    // This can reduce the number of if clauses that need to be tested.
    // Insert the order following the most abundant node types. This will minimize the number of if...else clauses to be evaluated.
    // 
    if(type == 1) {
      // Executes for node to tips contrast.
      // Former: NodeToTip function.
      des_node0 = des(node_id[0]) - 1; // des is a vector of indexes from R.
      des_node1 = des(node_id[1]) - 1;
      ss = X.col(des_node0) - X.col(des_node1);
    
      // Multiply each R matrix by the respective branch length (due to the regime) and sum the result. So this is a loop over the number of regimes 'p'.
      // Need to do this for each of the daughter lineages.
      for(arma::uword j = 0; j < p; j++) {
	// Note the stop condition for the loop. We can use the number of regimes here because the loop will stop when i=1 < p=2 !! And NOT when j=2 !!
	// Multiply each R matrix for the correspondent branch length.
	Rs1.slice(j) = R.slice(j) * mapped_edge(node_id[0],j);
	Rs2.slice(j) = R.slice(j) * mapped_edge(node_id[1],j);
      }
      // Join all slices together, sum everything into a single matrix.
      // Without copying it all again! Awesome!
      for(arma::uword z = 1; z < p; z++) {
	Rs1.slice(0) += Rs1.slice(z);
	Rs2.slice(0) += Rs2.slice(z);
      }
      Rinv = inv_sympd( Rs1.slice(0) + Rs2.slice(0) );

      ll = ll + logLikNode_C(ss, Rs1.slice(0) + Rs2.slice(0), Rinv, k);

      key[i] = anc(node_id[0]) - 1; // 'anc' is a vector from R, so need to change the indexing here.
      // Take care with the indexes starting from 0, need to reduce a unit from the length here.
      X0.col(i) = ((Rs1.slice(0) * Rinv) * X.col(des_node1)) + ((Rs2.slice(0) * Rinv)  * X.col(des_node0));

      V0.slice(i) = inv_sympd( inv_sympd(Rs1.slice(0)) + inv_sympd(Rs2.slice(0)) );
  
    } else if(type == 3) {
      // Executes for node to tip & node contrast.
      // Former NodeToTipNode function.
      // Here using a bifurcating tree always, so just need two objects.
      // des is a vector from R. Need to correct the indexes.
      des_node_vec[0] = des(node_id[0]) - 1;
      des_node_vec[1] = des(node_id[1]) - 1;
      // This will give the nodes associated with the descendants (from the key).
      // Here nd is the node of the tree.
      // Need to generate some indexes to be used by the function.
      nd = max( des_node_vec );
      tip = min( des_node_vec );
      nd_id = node_id[as_scalar( find( nd == des_node_vec ) )];
      tip_id = node_id[as_scalar( find( tip == des_node_vec ) )];
      // Here I am using key.head() to search only in the populated values of key.
      // The position 'i' will only be populated some lines down in the code.
      key_id = as_scalar( find( key.head(i) == nd ) );

      // Compute the contrast for the nodes.
      ss = X.col(tip) - X0.col(key_id);

      // 'Rs1' is relative to the branch leading to a tip 'tip_id' and 'Rs2' to the branch leading to a node 'nd_id'.
      for(arma::uword j = 0; j < p; j++) {
	// Note the stop condition for the loop. We can use the number of regimes here because the loop will stop when j=1 < p=2 !! And NOT when j=2 !!
	// Multiply each R matrix for the correspondent branch length.
	Rs1.slice(j) = R.slice(j) * mapped_edge(tip_id,j);
	Rs2.slice(j) = R.slice(j) * mapped_edge(nd_id,j);
      }
      // Rf_PrintValue( wrap( row2 ) );
      // Join all slices together, sum everything into a single matrix.
      // Without copying it all again! Awesome!
      for(arma::uword z = 1; z < p; z++) {
	Rs1.slice(0) += Rs1.slice(z);
	Rs2.slice(0) += Rs2.slice(z);
      }
      // Here need to add some addicional variance that carries over from the pruning.
      // Doing this just for the node. No additional variance associated with the tip.
      Rs2.slice(0) += V0.slice(key_id);

      Rinv = inv_sympd( Rs1.slice(0) + Rs2.slice(0) );

      ll = ll + logLikNode_C(ss, Rs1.slice(0) + Rs2.slice(0), Rinv, k);

      key[i] = anc(node_id[0]) - 1; // Need to decrease the value of 'anc' by 1. This comes from R.
  
      // Take care with the indexes starting from 0, need to reduce a unit from the length here.
      // Now we need to multiply the tip with the node and the node with the tip. That is why the relationship here is inverted. It is correct!
      X0.col(i) = ((Rs1.slice(0) * Rinv) * X0.col(key_id)) + ((Rs2.slice(0) * Rinv) * X.col(tip));
      V0.slice(i) = inv_sympd( inv_sympd(Rs1.slice(0)) + inv_sympd(Rs2.slice(0)) );
  
    } else {
      // Executes for node to nodes contrast.
      // Former NodeToNode function.
      // Here using a bifurcating tree always, so just need two objects.
      des_node0 = des(node_id[0]) - 1; // 'des' is a vector from R.
      des_node1 = des(node_id[1]) - 1;

      // This will give the nodes associated with the descendants (from the key).
      // This also depends on the correction of 'anc' when populating the key vector.
      // Using key.head(i) to assure that find is only looking to populated positions in the vector.
      key_id0 = as_scalar( find(key.head(i) == des_node0) );
      key_id1 = as_scalar( find(key.head(i) == des_node1) );
  
      // Compute the contrast for the nodes.
      ss = X0.col(key_id0) - X0.col(key_id1);
      // Multiply each R matrix by the respective branch length (due to the regime) and sum the result. So this is a loop over the number of regimes 'p'.
      // Need to do this for each of the daughter lineages.
      for(arma::uword j = 0; j < p; j++) {
	// Note the stop condition for the loop. We can use the number of regimes here because the loop will stop when j=1 < p=2 !! And NOT when j=2 !!
	// Multiply each R matrix for the correspondent branch length.
	Rs1.slice(j) = R.slice(j) * mapped_edge(node_id[0],j);
	Rs2.slice(j) = R.slice(j) * mapped_edge(node_id[1],j);
      }
      // Join all slices together, sum everything into a single matrix.
      // Without copying it all again! Awesome!
      for(arma::uword z = 1; z < p; z++) {
	Rs1.slice(0) += Rs1.slice(z);
	Rs2.slice(0) += Rs2.slice(z);
      }
      // Here need to add some addicional variance that carries over from the pruning.
      Rs1.slice(0) += V0.slice(key_id0);
      Rs2.slice(0) += V0.slice(key_id1);
  
      Rinv = inv_sympd( Rs1.slice(0) + Rs2.slice(0) );

      ll = ll + logLikNode_C(ss, Rs1.slice(0) + Rs2.slice(0), Rinv, k);
      
      // ## 'key is for both X0 and V0.
      // Here need to decrease in one unit the value of the 'anc' vector. This is a vector coming from R.
      key[i] = anc(node_id[0]) - 1;
  
      // Take care with the indexes starting from 0, need to reduce a unit from the length here.
      X0.col(i) = ((Rs1.slice(0) * Rinv) * X0.col(key_id1)) + ((Rs2.slice(0) * Rinv) * X0.col(key_id0));
      V0.slice(i) = inv_sympd( inv_sympd(Rs1.slice(0)) + inv_sympd(Rs2.slice(0)) );
    }

  }

  // At the end, need to compute the contrast for the ROOT.
  // Make the calculation for the root log-likelihood.
  // The index 'n_nodes' is correspondent to the position 'n_nodes + 1'. Remember the indexation starts from 0.
  ss = X0.col(n_nodes-1) - mu;
  ll = ll + logLikNode_C(ss, V0.slice(n_nodes-1), inv_sympd(V0.slice(n_nodes-1)), k);

  return(ll);
}

// #######################################################
// ###### The distribution functions
// #######################################################
// These are some distribution functions that are not available on Rcpp or on RcppArmadillo.

double logDensityIWish_C(arma::mat W, double v, arma::mat S){
  // Function for the density of a inverse Wishart distribution.
  // W the covariance matrix.
  // v the degrees of freedom parameter.
  // S the standard matrix for the distribution.
  double valS;
  double valW;
  double sign_sink;
  double lgammapart = 0;
    
  double k = S.n_cols;
  for(int i=0; i < S.n_cols; i++) {
    lgammapart = lgammapart + lgamma((v-i)/2);
  }
  log_det(valS, sign_sink, S);
  log_det(valW, sign_sink, W);
  
  double ldenom = lgammapart + ( ( (v*k)/2.0 ) * log( 2.0 ) ) + ( ( (k*(k-1.0))/4.0 ) * log( arma::datum::pi ) );
  // Need to make sure that we are doing 'double' operations here!
  double lnum = ( ( v/2.0 ) * valS ) + ( ( -(v + k + 1.0)/2.0 ) * valW ) + ( -0.5 * trace( S * inv_sympd(W) ) );
  return lnum - ldenom;
}

arma::mat riwish_C(int v, arma::mat S){
  // Generates a random draw from a inverse-Wishart distribution.
  // arma::mat CC = chol( inv_sympd(S) );
  arma::mat CC = chol( inv_sympd(S) );
  int p = S.n_cols;
  // Make a diagonal matrix with the elements:
  // R::rchisq( df ) // with df in a sequence v:(v - p + 1)
  arma::vec chi_sample(p);
  // regspace will make a sequence.
  arma::vec df = regspace(v, v-p+1); // The degrees of freedom when sampling.
  for( int i=0; i < p; i++ ) {
    chi_sample[i] = sqrt( R::rchisq( df[i] ) );
  }
  arma::mat Z = diagmat( chi_sample );
  // Need to fill the upper triangular elements with samples from a standard normal.
  for( int i=0; i < p-1; i++) {
    // randn uses a normal distribution to sample.
    Z(i,span((i+1), (p-1))) = trans(randn(p-(i+1)));
  }
  arma::mat out = Z * CC;
  return inv_sympd( trans(out) * out );
}

double hastingsDensity_C(arma::cube R, arma::cube R_prop, int k, arma::vec v, int Rp){
  // The hasting is only computed for the regime that is updated (Rp).
  arma::mat center_curr = (v[Rp]-k-1) * R.slice(Rp);
  arma::mat center_prop = (v[Rp]-k-1) * R_prop.slice(Rp);
  return logDensityIWish_C(R.slice(Rp), v[Rp], center_prop) - logDensityIWish_C(R_prop.slice(Rp), v[Rp], center_curr);
}  

// [[Rcpp::export]]
arma::mat cov2cor_C(arma::mat V){
  // This is a **brute force** function for the correlation matrix.
  arma::mat Vdiag = inv( sqrt( diagmat(V) ) );
  return Vdiag * V * Vdiag;
}

// #######################################################
// ###### The prior functions
// #######################################################

double priorRoot_C(arma::vec mu, arma::mat par_prior_mu, std::string den_mu){
  // Here the 'par_prior_mu' need to be a matrix constructed before.
  // The value for 'den_mu' will also be given before.
  // These objects will be constructed by 'makePrior' function.
  // If you have a problem with the vector type you can use:
  // NumericVector(a.begin(),a.end()) to transform into a Rcpp vector.
  double pp = 0.0;
  if( den_mu == "unif" ){
    for( int i=0; i < mu.n_elem; i++ ){
      pp = pp + R::dunif(mu[i], par_prior_mu(i,0), par_prior_mu(i,1), true);
    }
  } else{
    for( int i=0; i < mu.n_elem; i++ ){
      pp = pp + R::dnorm(mu[i], par_prior_mu(i,0), par_prior_mu(i,1), true);
    }
  }  
  return pp;
}

double priorSD_C(arma::mat sd, arma::mat par_prior_sd, std::string den_sd){
  // FUNCTION FOR 2 OR MORE RATE REGIMES.
  // 'sd' is a matrix with number of columns equal to 'p', number of regimes.
  // This function will work even if there is only a single regime, because Armadillo treats the vectors as column vector. So 'sd' will be a matrix with a single column.
  double pp = 0.0;
  if( den_sd == "unif" ){
    for( int i=0; i < sd.n_rows; i++ ){
      for( int j=0; j < sd.n_cols; j++){
	// Each line of the 'par_prior_sd' correspond to each of the sd vectors stored as the columns of the 'sd' matrix.
	pp = pp + R::dunif(sd(i,j), par_prior_sd(j,0), par_prior_sd(j,1), true);
      }
    }
  } else{
    for( int i=0; i < sd.n_rows; i++ ){
      for( int j=0; j < sd.n_cols; j++){
	// Each line of the 'par_prior_sd' correspond to each of the sd vectors stored as the columns of the 'sd' matrix.
	pp = pp + R::dlnorm(sd(i,j), par_prior_sd(j,0), par_prior_sd(j,1), true);
      }
    }
  }  
  return pp;
}

// [[Rcpp::export]]
double priorCorr_C(arma::cube corr, arma::vec nu, arma::cube sigma){
  // FUNCTION FOR 2 OR MORE RATE REGIMES.
  // nu is a vector with the correspondent degree of freedom for the rate regime.
  // This is the more complicated one. Need to deal with the Wishart distributions.
  // If the correlation prior was set to "uniform'. Then we just need to set sigma and v to the standard values when doing the 'makePrior' step. No need for a if test here.
  // Will treat the correlation and sigma as arrays (cube). This will work even if they are matrices. I think.
  int p = corr.n_slices;
  double pp = 0.0;
  for( int i=0; i < p; i++ ) {
    pp = pp + logDensityIWish_C(corr.slice(i), nu[i], sigma.slice(i)); // Need to define this one.
  }
  return pp;
}

// #######################################################
// ###### The proposal functions
// #######################################################

arma::vec multiplierProposal_C(int size, arma::vec w_sd){
  // A proposal that scales with the absolute value of the parameter. (Always positive.)
  // in R this is doing:
  // exp( 2 * log(w_sd) * runif(1, min=-0.5, max=0.5) )
  // Get this return factor and multiply by the current value for the proposal.
  // Return is a vector of length 'size', so it will work for multiple parameters.
  // The proposal ratio is the sum of the output vector.
  return exp( (randu(size) - 0.5) % (2.0 * log(w_sd)) );
}

arma::vec slideWindowLogSpace_C(arma::vec mu, arma::vec w_mu){
  // This will take the vector of root values mu and a vector of widths and make a draw.
  // Difference is that here we do the sliding window in log-space.
  // This proposal strategy is symmetrical in log-space but assymetrical in normal space.
  // Not a problem because the symmetry need to happen in the proposal step.
  // in R this is doing:
  // y <- sapply(1:length(mu), function(x) exp( log(mu) - (log(w_mu)/2) + (runif(1) * log(w_mu)) ) )
  return exp( log(mu) - (log(w_mu)/2) + (randu(mu.n_elem) % log(w_mu)) );
}

arma::vec slideWindow_C(arma::vec mu, arma::vec w_mu){
  // This will take the vector of root values mu and a vector of widths and make a draw.
  // in R this is doing:
  // y <- sapply(1:length(mu), function(x) runif(1, min = mu[i] - (w_mu[i]/2), max = mu[i] + (w_mu[i]/2) ) )
  return (mu - (w_mu/2)) + ( randu( mu.n_elem ) % w_mu );
}

// arma::vec slideWindowPositive_C(arma::vec sd, arma::vec w_sd){
//   // This will take the vector of root values mu and a vector of widths and make a draw.
//   arma::vec prop_sd = (sd - (w_sd/2)) + ( randu(sd.n_elem) % w_sd );
//   for(int i=0; i < sd.n_elem; i++){
//     prop_sd[i] = std::abs(prop_sd[i]);
//   }
//   return prop_sd;
// }

// [[Rcpp::export]]
arma::mat makePropIWish_C(arma::mat vcv, double k, double v){
  // Function to make a proposal for the correlation structure.
  // Proposal for a single R matrix, so vcv is not a cube.
  // Defining the parameters here as double because of the multiplication.
  // v here is the size if the step. This will usually be a large number.
  arma::mat center = (v-k-1) * vcv;
  // Need the riwish function here.
  return riwish_C(v, center);
}

// #######################################################
// ###### Function to write to file
// #######################################################

void writeToMultFile_C(std::ostream& mcmc_stream, int p, int k, arma::cube R, arma::vec mu){
  // Note the 'std::ostream&' argument here is the use of a reference.
  for( int i=0; i < p; i++ ){
    for( int j=0; j < k; j++ ){
      for( int z=0; z < k; z++ ){
	mcmc_stream << R.slice(i)(j,z);
	mcmc_stream << "; ";
      }
    }
  }
  for( int i=0; i < k-1; i++ ){
    mcmc_stream << mu[i];
    mcmc_stream << "; ";
  }
  mcmc_stream << mu.tail(1);
  // mcmc_stream << "\n";
}

// #######################################################
// ###### The MCMC function.
// #######################################################

// This is the MCMC function when there is only a single tree/regime configuration to be estimated. Allowing for multiple trees could be done in the same function. However, it is faster and simpler just to duplicate this function and assume that we are always working with a list of trees.
// The function 'runRatematrixMultiMCMC_C' is doing this.
// A similar need will happen for the MCMC with a single regime. In this case a bunch of computations are not needed. I can just simplify a lot this MCMC function. Using my own likelihood function will also means I can drop both 'mvMORPH' and 'phytools' as dependencies for the package. This will also make installation in a server much more user-friendly. [The 'rgl' dependecy of 'phytools' makes installation in servers difficult.]

// [[Rcpp::export]]
std::string runRatematrixMCMC_C(arma::mat X, int k, int p, arma::uvec nodes, arma::uvec des, arma::uvec anc, arma::uvec names_anc, arma::mat mapped_edge, arma::cube R, arma::vec mu, arma::mat sd, arma::cube Rcorr, arma::vec w_mu, arma::mat par_prior_mu, std::string den_mu, arma::mat w_sd, arma::mat par_prior_sd, std::string den_sd, arma::vec nu, arma::cube sigma, arma::vec v, std::string log_file, std::string mcmc_file, double prob_sample_root, double prob_sample_sd, int gen, int write_header){
  // The data parameters:
  // X, k, p, nodes, des, anc, names_anc, mapped_edge.
  // The starting point parameters. These are the objects to carry on the MCMC.
  // R, mu, sd, Rcorr
  // The starting point priors:
  // curr_root_prior, curr_sd_prior, Rcorr_curr_prior
  // The parameters for the priors:
  // par_prior_mu, den_mu, par_prior_sd, den_sd, nu, sigma
  // The starting point jacobian:
  // curr_jacobian
  // The parameters for the proposals:
  // w_mu, w_sd, v
  // The parameters to control the MCMC:
  // prob_sample_root, prob_sample_var, gen
  // write_header, wheather to write the header to file or just to append.

  // Open the files to write:
  // The log_file and mcmc_file arguments.
  std::ofstream log_stream (log_file, ios::out | ios::app);
  std::ofstream mcmc_stream (mcmc_file, ios::out | ios::app);

  // Write the header for the mcmc file.
  if(write_header == 1){
    for( int kk=1; kk < p+1; kk++ ){
      for( int ii=1; ii < k+1; ii++ ){
	for( int jj=1; jj < k+1; jj++ ){
	  mcmc_stream << "regime.p";
	  mcmc_stream << kk;
	  mcmc_stream << ".";
	  mcmc_stream << ii;
	  mcmc_stream << jj;
	  mcmc_stream << "; ";
	}
      }
    }
  
    for( int kk=1; kk < k; kk++ ){
      mcmc_stream << "trait.";
      mcmc_stream << kk;
      mcmc_stream << "; ";
    }
    mcmc_stream << "trait.";
    mcmc_stream << k;
    mcmc_stream << "\n";

    // Write the header for the log file.
    log_stream << "accepted; matrix.corr; matrix.sd; root; which.phylo; log.lik \n";
  } else{
    // Do nothing.
  }
  
  // Define the containers for the run:
  // The starting point log lik:
  double lik;
  // The starting point prior:
  double curr_root_prior;
  double curr_sd_prior;
  double Rcorr_curr_prior;
  // The jacobian for the MCMC. There are one for each regime.
  // Because need to track the jump separatelly.
  arma::vec curr_jacobian = vec(k, fill::zeros);

  int sample_root;
  int sample_sd;
  arma::vec prop_root;
  double prop_root_prior;
  double pp;
  double prop_root_lik;
  double ll;
  double r;
  double unif_draw;
  int Rp;
  arma::mat prop_sd;
  double prop_sd_prior;
  arma::mat prop_diag_sd;
  arma::cube R_prop;
  double prop_sd_lik;
  arma::cube Rcorr_prop;
  double Rcorr_prop_prior;
  arma::mat diag_sd;
  double prop_corr_lik;
  double hh;
  arma::vec decomp_var;
  double prop_jacobian = 0.0; // Need always to reset this one.
  double jj;
  arma::vec multi_factor; // Stores the multiplier proposal scaler.

  // Get starting priors, likelihood, and jacobian.
  lik = logLikPrunningMCMC_C(X, k, p, nodes, des, anc, names_anc, mapped_edge, R, mu);
  curr_root_prior = priorRoot_C(mu, par_prior_mu, den_mu);
  curr_sd_prior = priorSD_C(sd, par_prior_sd, den_sd);
  Rcorr_curr_prior = priorCorr_C(Rcorr, nu, sigma);

  Rcout << "Starting point Log-likelihood: " << lik << "\n";

  arma::mat var_vec = square(sd);
  // Jacobian for both the regimes:
  for( int j=0; j < p; j++ ){
    for( int i=0; i < k; i++ ){
      // The jacobian is computed on the variances!
      curr_jacobian[j] = curr_jacobian[j] + ( log( var_vec.col(j)[j] ) * log( (k-1.0)/2.0 ) );
    }
  }

  // Print starting point to files:
  log_stream << "1; 0; 0; 1; 1; ";
  log_stream << lik;
  log_stream << "\n"; 
  writeToMultFile_C(mcmc_stream, p, k, R, mu);

  Rcout << "Starting MCMC ... \n";
  // Starting the MCMC.
  for( int i=0; i < gen; i++ ){
    // Sample between root and matrix. Success here will be the update of the root.
    sample_root = R::rbinom(1, prob_sample_root);
    // If the matrix is updated, do we update the vector of variances?
    sample_sd = R::rbinom(1, prob_sample_sd);
    // If matrix or variance is updated, which regime?
    Rp = Rcpp::as<int >(Rcpp::sample(p, 1)) - 1; // Need to subtract 1 from the result here.
    // Rp = as_scalar( randi(1, distr_param(0, p-1)) ); // The index!
	  
    // The 'if...if...else' structure is lazy and will evaluate only the first entry.
    if(sample_root == 1){
      // Update the root state.
      prop_root = slideWindow_C(mu, w_mu);
      // Compute the prior.
      prop_root_prior = priorRoot_C(prop_root, par_prior_mu, den_mu);
      // The log prior ratio. curr_root_prior is the current prior density for the root.
      pp = prop_root_prior - curr_root_prior;
      prop_root_lik = logLikPrunningMCMC_C(X, k, p, nodes, des, anc, names_anc, mapped_edge, R, prop_root);
      // The likelihood ratio.
      ll = prop_root_lik - lik;
      // Get the ratio in log space.
      r = ll + pp;

      // Advance to the acceptance step.
      // Here we are only updating the root, so all other parameters are the same.
      unif_draw = as_scalar(randu(1)); // The draw from a uniform distribution.
      if( exp(r) > unif_draw ){ // Accept.
	log_stream << "1; 0; 0; 1; 1; ";
	log_stream << prop_root_lik;
	log_stream << "\n";
	mu = prop_root; // Update the value for the root. Need to carry over.
	curr_root_prior = prop_root_prior; // Update the root prior. Need to carry over.
	lik = prop_root_lik; // Update likelihood. Need to carry over.
      } else{ // Reject. Keep the values the same.
	log_stream << "0; 0; 0; 1; 1; ";
	log_stream << lik;
	log_stream << "\n";
      }
    } else if(sample_sd == 1){
      // Update the variance vector.
      // Draw which regime to update.
      // Rp = as_scalar( randi(1, distr_param(0, p-1)) ); // The index!
      prop_sd = sd; // The matrix of standard deviations.
      // prop_sd.col(Rp) = slideWindowPositive_C(prop_sd.col(Rp), w_sd.col(Rp));
      multi_factor = multiplierProposal_C(k, w_sd.col(Rp) ); // The factor for proposal. Also proposal ratio.
      prop_sd.col(Rp) = prop_sd.col(Rp) % multi_factor;
      prop_sd_prior = priorSD_C(prop_sd, par_prior_sd, den_sd);
      pp = prop_sd_prior - curr_sd_prior;
  
      // Need to rebuild the R matrix to compute the likelihood:
      prop_diag_sd = diagmat( prop_sd.col(Rp) );
      // The line below reconstructs the covariance matrix from the variance vector and the correlation matrix.
      R_prop = R; // R is the current R matrix.
      R_prop.slice(Rp) = prop_diag_sd * cov2cor_C( R.slice(Rp) ) * prop_diag_sd;
      prop_sd_lik = logLikPrunningMCMC_C(X, k, p, nodes, des, anc, names_anc, mapped_edge, R_prop, mu);
      ll = prop_sd_lik - lik;

      // Get the ratio i log space. Loglik, log prior and the proposal ratio (for the multiplier!).
      r = ll + pp + accu(multi_factor);

      // Advance to the acceptance step.
      // Here we are only updating the root, so all other parameters are the same.
      unif_draw = as_scalar(randu(1)); // The draw from a uniform distribution.
      if( exp(r) > unif_draw ){ // Accept.
	log_stream << "1; 0; ";
	log_stream << Rp+1; // Here is the regime.
	log_stream << "; 0; 1; ";
	log_stream << prop_sd_lik;
	log_stream << "\n";
	R = R_prop; // Update the evolutionary rate matrices. Need to carry over.
	sd = prop_sd; // Update the standard deviation.
	curr_sd_prior = prop_sd_prior;  // Update the prior. Need to carry over.
	lik = prop_sd_lik; // Update likelihood. Need to carry over.
      } else{ // Reject. Keep the values the same.
	log_stream << "0; 0; ";
	log_stream << Rp+1; // Here is the regime.
	log_stream << "; 0; 1; ";
	log_stream << lik;
	log_stream << "\n";
      }
      
    } else{
      // Update the correlation matrix.
      // Draw which regime to update.
      // Rp = as_scalar( randi(1, distr_param(0, p-1)) ); // The index!
      // IMPORTANT: R here needs to be the vcv used to draw the correlation only.
      // This is not the full VCV of the model.
      Rcorr_prop = Rcorr;
      // Here v is a vector. So the width can be controlled for each regime.
      Rcorr_prop.slice(Rp) = makePropIWish_C(Rcorr.slice(Rp), k, v[Rp]);
      // Computes the hastings density for this move.
      // Hastings is computed for the covariance matrix of the correlation and NOT the evolutionary rate matrix (the reconstructed).
      hh = hastingsDensity_C(Rcorr, Rcorr_prop, k, v, Rp);
      // Need the prior parameters 'nu' and 'sigma' here.
      Rcorr_prop_prior = priorCorr_C(Rcorr_prop, nu, sigma);
      pp = Rcorr_prop_prior - Rcorr_curr_prior;

      // Need to rebuild the R matrix to compute the likelihood:
      diag_sd = diagmat( sd.col(Rp) );
      // The line below reconstructs the covariance matrix from the variance vector and the correlation matrix.
      R_prop = R;
      // 'cor' is not the correct function to use here!
      R_prop.slice(Rp) = diag_sd * cov2cor_C( Rcorr_prop.slice(Rp) ) * diag_sd;
      prop_corr_lik = logLikPrunningMCMC_C(X, k, p, nodes, des, anc, names_anc, mapped_edge, R_prop, mu);
      ll = prop_corr_lik - lik;
      // This is the sdiance of the proposed vcv matrix.
      decomp_var = diagvec( Rcorr_prop.slice(Rp) ); // These are variances
      prop_jacobian = 0.0; // Always need to reset.
      // The Jacobian of the transformation.
      for( int i=0; i < k; i++ ){
	prop_jacobian = prop_jacobian + log(decomp_var[i]) * log( (k-1.0)/2.0 );
      }
      // The curr_jacobian is a vector. There are a jacobian for each regime.
      jj = prop_jacobian - curr_jacobian[Rp];

      // Get the ratio i log space. Loglik, log prior, log hastings and log jacobian.
      r = ll + pp + hh + jj;
      // Advance to the acceptance step.
      // Here we are only updating the root, so all other parameters are the same.
      unif_draw = as_scalar(randu(1)); // The draw from a uniform distribution.
      if( exp(r) > unif_draw ){ // Accept.
	// This line will write to the mcmc_file.
	// Instead of 'paste' I am using a line for each piece. Should have the same effect.
	log_stream << "1; ";
	log_stream << Rp+1; // Here is the regime.
	log_stream << "; 0; 0; 1; ";
	log_stream << prop_corr_lik;
	log_stream << "\n";
	R = R_prop; // Update the evolutionary rate matrices. Need to carry over.
	Rcorr = Rcorr_prop; // Update the correlation matrix.
	Rcorr_curr_prior = Rcorr_prop_prior; // Update the prior. Need to carry over.
	lik = prop_corr_lik; // Update likelihood. Need to carry over.
	curr_jacobian[Rp] = prop_jacobian; // Updates jacobian.
      } else{ // Reject. Keep the values the same.
	log_stream << "0; ";
	log_stream << Rp+1; // Here is the regime.
	log_stream << "; 0; 0; 1; ";
	log_stream << lik;
	log_stream << "\n";
      }
    }

    // Write the current state to the MCMC file.
    writeToMultFile_C(mcmc_stream, p, k, R, mu);
    
  }

  Rcout << "Closing files... \n";

  mcmc_stream.close();
  log_stream.close();

  return "Done.";
}

// The MCMC function for multiple regimes:

// [[Rcpp::export]]
std::string runRatematrixMultiMCMC_C(arma::mat X, int k, int p, arma::umat nodes, arma::umat des, arma::umat anc, arma::umat names_anc, arma::cube mapped_edge, arma::cube R, arma::vec mu, arma::mat sd, arma::cube Rcorr, arma::vec w_mu, arma::mat par_prior_mu, std::string den_mu, arma::mat w_sd, arma::mat par_prior_sd, std::string den_sd, arma::vec nu, arma::cube sigma, arma::vec v, std::string log_file, std::string mcmc_file, double prob_sample_root, double prob_sample_sd, int gen, int write_header){
  // The data parameters:
  // X, k, p, nodes, des, anc, names_anc, mapped_edge.
  // These parameters changed from the 'runRatematrixMCMC_C' function:
  // node, des, anc, and names_anc are now matrices. Each column of these matrices are the vector information for a given tree/regime configuration.
  // Also mapped_edge is now a cube. Each element of the cube is a different regime configuration for each of the trees.
  // So the function will, at each MCMC step, sample one of the trees and use that tree to run the MCMC.
  // The number of trees in the data will be computed from the number of elements in the mapped_edge cube.
  // The starting point parameters. These are the objects to carry on the MCMC.
  // R, mu, sd, Rcorr
  // The starting point priors:
  // curr_root_prior, curr_sd_prior, Rcorr_curr_prior
  // The parameters for the priors:
  // par_prior_mu, den_mu, par_prior_sd, den_sd, nu, sigma
  // The starting point jacobian:
  // curr_jacobian
  // The parameters for the proposals:
  // w_mu, w_sd, v
  // The parameters to control the MCMC:
  // prob_sample_root, prob_sample_var, gen
  // write_header, wheather to write the header to file or just to append.

  // Open the files to write:
  // The log_file and mcmc_file arguments.
  std::ofstream log_stream (log_file, ios::out | ios::app);
  std::ofstream mcmc_stream (mcmc_file, ios::out | ios::app);

  // Write the header for the mcmc file.
  if(write_header == 1){
    for( int kk=1; kk < p+1; kk++ ){
      for( int ii=1; ii < k+1; ii++ ){
	for( int jj=1; jj < k+1; jj++ ){
	  mcmc_stream << "regime.p";
	  mcmc_stream << kk;
	  mcmc_stream << ".";
	  mcmc_stream << ii;
	  mcmc_stream << jj;
	  mcmc_stream << "; ";
	}
      }
    }
  
    for( int kk=1; kk < k; kk++ ){
      mcmc_stream << "trait.";
      mcmc_stream << kk;
      mcmc_stream << "; ";
    }
    mcmc_stream << "trait.";
    mcmc_stream << k;
    mcmc_stream << "\n";

    // Write the header for the log file.
    log_stream << "accepted; matrix.corr; matrix.sd; root; which.phylo; log.lik \n";
  } else{
    // Do nothing.
  }
    
  // Define the containers for the run:
  // The starting point log lik:
  double lik;
  // The starting point prior:
  double curr_root_prior;
  double curr_sd_prior;
  double Rcorr_curr_prior;
  // The jacobian for the MCMC. There are one for each regime.
  // Because need to track the jump separatelly.
  arma::vec curr_jacobian = vec(k, fill::zeros);

  // Compute the number of trees in the data:
  int n_trees = mapped_edge.n_slices;
  int curr_phy = 0; // This is the current phy. It starts with the first tree.
  int prop_phy; // Container to store the sampled tree.
  
  int sample_root;
  int sample_sd;
  arma::vec prop_root;
  double prop_root_prior;
  double pp;
  double prop_root_lik;
  double ll;
  double r;
  double unif_draw;
  int Rp;
  arma::mat prop_sd;
  double prop_sd_prior;
  arma::mat prop_diag_sd;
  arma::cube R_prop;
  double prop_sd_lik;
  arma::cube Rcorr_prop;
  double Rcorr_prop_prior;
  arma::mat diag_sd;
  double prop_corr_lik;
  double hh;
  arma::vec decomp_var;
  double prop_jacobian = 0.0; // Need always to reset this one.
  double jj;
  arma::vec multi_factor; // The scale factor for the multiplier proposal.

  // Get starting priors, likelihood, and jacobian.
  // Here using the first tree (the curr_phy).
  lik = logLikPrunningMCMC_C(X, k, p, nodes.col(curr_phy), des.col(curr_phy), anc.col(curr_phy), names_anc.col(curr_phy), mapped_edge.slice(curr_phy), R, mu);
  curr_root_prior = priorRoot_C(mu, par_prior_mu, den_mu);
  curr_sd_prior = priorSD_C(sd, par_prior_sd, den_sd);
  Rcorr_curr_prior = priorCorr_C(Rcorr, nu, sigma);

  Rcout << "Starting point Log-likelihood (using first tree/regime in the list): " << lik << "\n";

  arma::mat var_vec = square(sd);
  // Jacobian for both the regimes:
  for( int j=0; j < p; j++ ){
    for( int i=0; i < k; i++ ){
      // The jacobian is computed on the variances!
      curr_jacobian[j] = curr_jacobian[j] + ( log( var_vec.col(j)[j] ) * log( (k-1.0)/2.0 ) );
    }
  }

  // Print starting point to files:
  log_stream << "1; 0; 0; 1; 1; ";
  log_stream << lik;
  log_stream << "\n"; 
  writeToMultFile_C(mcmc_stream, p, k, R, mu);

  Rcout << "Starting MCMC ... \n";
  // Starting the MCMC.
  for( int i=0; i < gen; i++ ){
    // Sample between root and matrix. Success here will be the update of the root.
    sample_root = R::rbinom(1, prob_sample_root);
    // If the matrix is updated, do we update the vector of variances?
    sample_sd = R::rbinom(1, prob_sample_sd);
    // If matrix or variance is updated, which regime?
    Rp = Rcpp::as<int >(Rcpp::sample(p, 1)) - 1; // Need to subtract 1 from the result here.
    // Rp = as_scalar( randi(1, distr_param(0, p-1)) ); // The index!
    // Which tree/regime configuration to use?
    prop_phy = Rcpp::as<int >(Rcpp::sample(n_trees, 1)) - 1; // This is an index too.
	  
    // The 'if...if...else' structure is lazy and will evaluate only the first entry.
    if(sample_root == 1){
      // Update the root state.
      // prop_root = slideWindow_C(mu, w_mu);
      prop_root = slideWindow_C(mu, w_mu);
      // Compute the prior.
      prop_root_prior = priorRoot_C(prop_root, par_prior_mu, den_mu);
      // The log prior ratio. curr_root_prior is the current prior density for the root.
      pp = prop_root_prior - curr_root_prior;
      prop_root_lik = logLikPrunningMCMC_C(X, k, p, nodes.col(prop_phy), des.col(prop_phy), anc.col(prop_phy), names_anc.col(prop_phy), mapped_edge.slice(prop_phy), R, prop_root);
      // The likelihood ratio.
      ll = prop_root_lik - lik;
      // Get the ratio in log space.
      r = ll + pp;

      // Advance to the acceptance step.
      // Here we are only updating the root, so all other parameters are the same.
      unif_draw = as_scalar(randu(1)); // The draw from a uniform distribution.
      if( exp(r) > unif_draw ){ // Accept.
	log_stream << "1; 0; 0; 1; ";
	log_stream << prop_phy+1; // Sum back to the number of the tree.
	log_stream << "; ";
	log_stream << prop_root_lik;
	log_stream << "\n";
	mu = prop_root; // Update the value for the root. Need to carry over.
	curr_root_prior = prop_root_prior; // Update the root prior. Need to carry over.
	lik = prop_root_lik; // Update likelihood. Need to carry over.
	curr_phy = prop_phy; // Update the current tree in the pool.
      } else{ // Reject. Keep the values the same.
	log_stream << "1; 0; 0; 1; ";
	log_stream << curr_phy+1; // Sum back to the number of the tree.
	log_stream << "; ";
	log_stream << lik;
	log_stream << "\n";
      }
    } else if(sample_sd == 1){
      // Update the variance vector.
      // Draw which regime to update.
      // Rp = as_scalar( randi(1, distr_param(0, p-1)) ); // The index!
      prop_sd = sd; // The matrix of standard deviations.
      // prop_sd.col(Rp) = slideWindowPositive_C(prop_sd.col(Rp), w_sd.col(Rp));
      multi_factor = multiplierProposal_C(k, w_sd.col(Rp) ); // The factor for proposal. Also proposal ratio.
      prop_sd.col(Rp) = prop_sd.col(Rp) % multi_factor;
      prop_sd_prior = priorSD_C(prop_sd, par_prior_sd, den_sd);
      pp = prop_sd_prior - curr_sd_prior;
  
      // Need to rebuild the R matrix to compute the likelihood:
      prop_diag_sd = diagmat( prop_sd.col(Rp) );
      // The line below reconstructs the covariance matrix from the variance vector and the correlation matrix.
      R_prop = R; // R is the current R matrix.
      R_prop.slice(Rp) = prop_diag_sd * cov2cor_C( R.slice(Rp) ) * prop_diag_sd;
      prop_sd_lik = logLikPrunningMCMC_C(X, k, p, nodes.col(prop_phy), des.col(prop_phy), anc.col(prop_phy), names_anc.col(prop_phy), mapped_edge.slice(prop_phy), R_prop, mu);
      ll = prop_sd_lik - lik;

      // Get the ratio i log space. Loglik, log prior and the proposal ratio (for the multiplier!).
      r = ll + pp + accu(multi_factor);

      // Advance to the acceptance step.
      // Here we are only updating the root, so all other parameters are the same.
      unif_draw = as_scalar(randu(1)); // The draw from a uniform distribution.
      if( exp(r) > unif_draw ){ // Accept.
	log_stream << "1; 0; ";
	log_stream << Rp+1; // Here is the regime.
	log_stream << "; 0; ";
	log_stream << prop_phy+1; // Add 1 to get number of phy.
	log_stream << "; ";
	log_stream << prop_sd_lik;
	log_stream << "\n";
	R = R_prop; // Update the evolutionary rate matrices. Need to carry over.
	sd = prop_sd; // Update the standard deviation.
	curr_sd_prior = prop_sd_prior;  // Update the prior. Need to carry over.
	lik = prop_sd_lik; // Update likelihood. Need to carry over.
	curr_phy = prop_phy; // Update the phylogeny from the pool.
      } else{ // Reject. Keep the values the same.
	log_stream << "0; 0; ";
	log_stream << Rp+1; // Here is the regime.
	log_stream << "; 0; ";
	log_stream << curr_phy+1; // Add 1 to get number of phy.
	log_stream << "; ";
	log_stream << lik;
	log_stream << "\n";
      }
      
    } else{
      // Update the correlation matrix.
      // Draw which regime to update.
      // Rp = as_scalar( randi(1, distr_param(0, p-1)) ); // The index!
      // IMPORTANT: R here needs to be the vcv used to draw the correlation only.
      // This is not the full VCV of the model.
      Rcorr_prop = Rcorr;
      // Here v is a vector. So the width can be controlled for each regime.
      Rcorr_prop.slice(Rp) = makePropIWish_C(Rcorr.slice(Rp), k, v[Rp]);
      // Computes the hastings density for this move.
      // Hastings is computed for the covariance matrix of the correlation and NOT the evolutionary rate matrix (the reconstructed).
      hh = hastingsDensity_C(Rcorr, Rcorr_prop, k, v, Rp);
      // Need the prior parameters 'nu' and 'sigma' here.
      Rcorr_prop_prior = priorCorr_C(Rcorr_prop, nu, sigma);
      pp = Rcorr_prop_prior - Rcorr_curr_prior;

      // Need to rebuild the R matrix to compute the likelihood:
      diag_sd = diagmat( sd.col(Rp) );
      // The line below reconstructs the covariance matrix from the variance vector and the correlation matrix.
      R_prop = R;
      // 'cor' is not the correct function to use here!
      R_prop.slice(Rp) = diag_sd * cov2cor_C( Rcorr_prop.slice(Rp) ) * diag_sd;
      prop_corr_lik = logLikPrunningMCMC_C(X, k, p, nodes.col(prop_phy), des.col(prop_phy), anc.col(prop_phy), names_anc.col(prop_phy), mapped_edge.slice(prop_phy), R_prop, mu);
      ll = prop_corr_lik - lik;
      // This is the sdiance of the proposed vcv matrix.
      decomp_var = diagvec( Rcorr_prop.slice(Rp) ); // These are variances
      prop_jacobian = 0.0; // Always need to reset.
      // The Jacobian of the transformation.
      for( int i=0; i < k; i++ ){
	prop_jacobian = prop_jacobian + log(decomp_var[i]) * log( (k-1.0)/2.0 );
      }
      // The curr_jacobian is a vector. There are a jacobian for each regime.
      jj = prop_jacobian - curr_jacobian[Rp];

      // Get the ratio i log space. Loglik, log prior, log hastings and log jacobian.
      r = ll + pp + hh + jj;
      // Advance to the acceptance step.
      // Here we are only updating the root, so all other parameters are the same.
      unif_draw = as_scalar(randu(1)); // The draw from a uniform distribution.
      if( exp(r) > unif_draw ){ // Accept.
	// This line will write to the mcmc_file.
	// Instead of 'paste' I am using a line for each piece. Should have the same effect.
	log_stream << "1; ";
	log_stream << Rp+1; // Here is the regime.
	log_stream << "; 0; 0; ";
	log_stream << prop_phy+1; // Add 1 to get the number of the phylo.
	log_stream << "; ";
	log_stream << prop_corr_lik;
	log_stream << "\n";
	R = R_prop; // Update the evolutionary rate matrices. Need to carry over.
	Rcorr = Rcorr_prop; // Update the correlation matrix.
	Rcorr_curr_prior = Rcorr_prop_prior; // Update the prior. Need to carry over.
	lik = prop_corr_lik; // Update likelihood. Need to carry over.
	curr_jacobian[Rp] = prop_jacobian; // Updates jacobian.
	curr_phy = prop_phy; // Updates the current phy from the pool.
      } else{ // Reject. Keep the values the same.
	log_stream << "0; ";
	log_stream << Rp+1; // Here is the regime.
	log_stream << "; 0; 0; ";
	log_stream << curr_phy+1; // Add 1 to get the number of the phylo.
	log_stream << "; ";
	log_stream << lik;
	log_stream << "\n";
      }
    }

    // Write the current state to the MCMC file.
    writeToMultFile_C(mcmc_stream, p, k, R, mu);
    
  }

  Rcout << "Closing files... \n";

  mcmc_stream.close();
  log_stream.close();

  return "Done.";
}
