// Functions for making stochastic maps.

#ifndef stochastic_mapping_h
#define stochastic_mapping_h

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// #######################################################
// ###### The stochastic mapping functions
// #######################################################

arma::mat getReconStates(arma::uword n_nodes, int n_tips, int n_states, arma::vec edge_len, arma::mat edge_mat, arma::vec parents, arma::mat X, arma::mat Q, int root_node, int root_type) {
  // This is the same as the logLikFunction but it returns a matrix with the probabilities for each of the states at each node of the phylogeny.

  // n_nodes = number of nodes in phy
  // n_tips = number of tips in phy
  // n_states = number of states in the regime.
  // edge_len = vector with the edge length
  // edge_mat = the edge matrix for the tree.
  // parents = vector with the unique( edge.matrix[,1] )
  // X = data matrix. number of columns equal to the number of states. rows in the same order as the tip.labels. a 1 marks state is present and a 0 mark state absent.
  // Q = the transition matrix.
  // root_node = the node number for the root node
  // root_type = the type to compute the root probabilities. 0 = equal and 1 = madfitz

  arma::mat append_mat = mat(n_nodes, n_states, fill::zeros);
  arma::mat liks = join_vert(X, append_mat);
  arma::vec comp = vec(n_nodes + n_tips);
  arma::uword anc;
  arma::uvec ii; // This is a vector of indexes.
  arma::mat v = mat(n_states, 2); // Two descendant nodes.

  // Loop to traverse the tree.
  for(uword i=0; i < n_nodes; i++) {

    // Need to check the usage of 'anc'. Is it an index or a vector test?
    
    anc = parents[i] - 1; // This is an index. C++ starts from 0.
    ii = find( parents[i] == edge_mat.col(0) ); // More than one entry.
    
    uword des;
    arma::vec v_root = vec(n_states, fill::ones);
    for(uword j=0; j < 2; j++) {
      des = as_scalar( edge_mat(ii[j], 1) ) - 1; // This is an index
      v.col(j) = expmat(Q * edge_len[ ii[j] ]) * trans( liks.row( des ) );
      v_root = v_root % v.col(j);
    }
	  
    if( parents[i] == root_node ){ // The computations at the root
      if( root_type == 0 ){
	// This is the equal root probability model:
	arma::vec equal_pi = vec(n_states, fill::ones);
	equal_pi = equal_pi / n_states;
	comp[ anc ] = sum( v_root % equal_pi );
      } else{
	// if( root_type == 1) // This needs to be TRUE
	// This is the Maddison and Fitzjohn method.
	// arma::vec comp_unscaled = sum( v_root.col(0) );
	arma::vec liks_root = v_root / sum( v_root );
	arma::vec root_p = liks_root / sum( liks_root );
	comp[ anc ] = sum( v_root % root_p );
      }
    } else{
      comp[ anc ] = sum( v_root );
    }

    liks.row(anc) = trans( v_root / comp[ anc ] ); // Need row vector.
  }

  // Get the probabilities for the states at each node:
  // Make sure that they all sum to 1
  arma::mat recon_states = mat(liks);
  for(uword mr=0; mr < recon_states.n_rows; mr++) {
    recon_states.row(mr) = liks.row(mr) / sum( liks.row(mr) );
  }
  
  return recon_states;
}


#endif
