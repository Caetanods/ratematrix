// Function to compute the likelihood of the mvBM model.

#ifndef likelihood_pruning_h
#define likelihood_pruning_h

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double logLikNode_C(arma::vec ss, arma::mat sigma_len, arma::mat sigma_len_inv, int k) {
  double val;
  double signal;
  arma::log_det(val, signal, sigma_len); // val is logdeterminant and signal should be 0.
  return -0.5 * ( k * log(2 * arma::datum::pi) + val + as_scalar(trans(ss) * sigma_len_inv * ss));
}

// Function to compute the log-likelihood for the model using the restricted likelihood.
// The main difference here is that the likelihood function is computed by averaging across all possible root states.

// Below is a function very similar to the one above. The difference is that we do not sample values for the inernal nodes. We just added the step of sampling tip values from a distribution. The idea here is that we are evaluating the likelihood of the model by marginalizing (integrating over) the values for the ancestral values of each of the nodes.

double logLikPrunningREML(arma::mat X, arma::uword k, arma::uword p, arma::vec nodes, arma::uvec des, arma::uvec anc, arma::uvec names_anc, arma::mat mapped_edge, arma::cube R) {

  // X: need to have the traits at the rows and the species at the columns.
  
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
  arma::uword n_nodes = nodes.n_elem;
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
  for(arma::uword i=0; i < n_nodes; i++) {
    
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

      Rinv = inv( Rs1.slice(0) + Rs2.slice(0) );

      ll = ll + logLikNode_C(ss, Rs1.slice(0) + Rs2.slice(0), Rinv, k);

      key[i] = anc(node_id[0]) - 1; // 'anc' is a vector from R, so need to change the indexing here.
      // Take care with the indexes starting from 0, need to reduce a unit from the length here.
      X0.col(i) = ((Rs1.slice(0) * Rinv) * X.col(des_node1)) + ((Rs2.slice(0) * Rinv)  * X.col(des_node0));

      V0.slice(i) = inv( inv(Rs1.slice(0)) + inv(Rs2.slice(0)) );
  
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

      Rinv = inv( Rs1.slice(0) + Rs2.slice(0) );

      ll = ll + logLikNode_C(ss, Rs1.slice(0) + Rs2.slice(0), Rinv, k);

      key[i] = anc(node_id[0]) - 1; // Need to decrease the value of 'anc' by 1. This comes from R.
  
      // Take care with the indexes starting from 0, need to reduce a unit from the length here.
      // Now we need to multiply the tip with the node and the node with the tip. That is why the relationship here is inverted. It is correct!
      X0.col(i) = ((Rs1.slice(0) * Rinv) * X0.col(key_id)) + ((Rs2.slice(0) * Rinv) * X.col(tip));
      V0.slice(i) = inv( inv(Rs1.slice(0)) + inv(Rs2.slice(0)) );
  
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
  
      Rinv = inv( Rs1.slice(0) + Rs2.slice(0) );

      ll = ll + logLikNode_C(ss, Rs1.slice(0) + Rs2.slice(0), Rinv, k);
      
      // ## 'key is for both X0 and V0.
      // Here need to decrease in one unit the value of the 'anc' vector. This is a vector coming from R.
      key[i] = anc(node_id[0]) - 1;
  
      // Take care with the indexes starting from 0, need to reduce a unit from the length here.
      X0.col(i) = ((Rs1.slice(0) * Rinv) * X0.col(key_id1)) + ((Rs2.slice(0) * Rinv) * X0.col(key_id0));
      V0.slice(i) = inv( inv(Rs1.slice(0)) + inv(Rs2.slice(0)) );
    }

  }

  // No contrast at the root because we are only computing the REML.
  return(ll);
}


#endif
