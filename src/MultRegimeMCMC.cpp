#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Source file calls multiple C++ headers.
// You can finde the headers with the C++ code on: inst/include/*.h
// All exported functions need to be defined here.

// #######################################################
// ###### Supporting functions
// #######################################################

// Include a series of supporting functions.
// These have both statistical distributions and some helping functions.
#include <distributions_and_helpers.h>

// #######################################################
// ###### Stochastic mappings
// #######################################################

// Include functions to make stochastic maps:
#include <stochastic_mapping.h>

// [[Rcpp::export]]
arma::mat makeSimmapMappedEdge(arma::uword n_nodes, arma::uword n_tips, arma::uword n_states, arma::vec edge_len, arma::mat edge_mat, arma::vec parents, arma::mat X, arma::mat Q, int root_node, bool root_type, int sims_limit) {
  // This function will not return the 'maps' element of the stochastic maps.
  // For this see the function below. This one is all that is needed for computing the likelihood of the BM model.
  // Function works by assuming that the order of the rows in sim_node_states is the same as in recon_states below. This follows because of the code in 'getReconStates' function.

  // NOTE: sims_limit controls the maximum number of times the stochastic map simulation will try to simulate for any particular branch (not cumulative!). If the value is 0, then no limit is stablished.

  // Define containers:
  int nrow_mapped_edge = edge_mat.n_rows;
  double time_chunk;
  arma::mat mapped_edge = mat(nrow_mapped_edge, n_states);
  arma::mat exp_Q = mat(Q);
  arma::uword anc_state;
  arma::uword des_state;
  arma::rowvec p_des_state;
  arma::uvec ii; // Index for many of the loops.
  arma::mat sim_node_states = mat(nrow_mapped_edge, 2);
  arma::vec prob_root = vec(n_states);
  arma::uword des_node;
  arma::uword curr_state;
  arma::uword sample_index;
  arma::umat trans_table_index = umat(n_states, n_states-1); // Take find output.
  arma::mat trans_table_prob = mat(n_states, n_states-1);
  arma::rowvec Q_row = rowvec(n_states);
  
  // Compute the probability for each state at each internal node given the model parameters. These are the conditional probabilities. Not the marginals.
  arma::mat recon_states = getReconStates(n_nodes, n_tips, n_states, edge_len, edge_mat, parents, X, Q, root_node, root_type);

  // Creates a table of transition probabilities given each starting state:
  // These only depend on Q and don't need to be computed more than once.
  for( arma::uword j=0; j < n_states; j++ ){
    arma::urowvec index_vec = urowvec(n_states);
    for( arma::uword jj=0; jj < n_states; jj++ ){
      index_vec[jj] = jj;
    }
    index_vec.shed_col(j);
    trans_table_index.row(j) = index_vec;
    // trans_table_index.row(j) = trans( find( Q_row > 0 ) );
    Q_row = Q.row(j); // This is a column vector type.
    Q_row.shed_col(j);
    trans_table_prob.row(j) = Q_row / sum( Q_row );
  }
  
  // Sample the state at the root.
  if( root_type ){ // 'madfitz'
    // Here the vector of probabilities is the same as the estimated under the model.
    arma::vec pi = trans( recon_states.row( root_node - 1 ) );
    prob_root = pi % pi;
  } else{ // root is 'equal'
    // Here the vector of probabilities is the same for all states.
    arma::vec pi = vec(n_states, fill::ones) / n_states;
    prob_root = trans( recon_states.row( root_node - 1 ) ) % pi;
  }
  anc_state = rMultinom( prob_root );
  ii = find( root_node == edge_mat.col(0) ); // Return length 2
  sim_node_states(ii[0],0) = anc_state;
  sim_node_states(ii[1],0) = anc_state;

  // Now we can sample all other states:
  // I will follow the reverse order of the post-order traversal edge matrix. This will make sure we are going from the root to the tips of the tree.
  for(int i=(nrow_mapped_edge-1); i >= 0; i-- ){
    
    exp_Q = expmat(Q * edge_len[i]);
    anc_state = sim_node_states(i,0); // This was already populated.
    des_node = edge_mat(i,1);
    p_des_state = exp_Q.row( anc_state ) % recon_states.row( des_node-1 ); // rowvec
    des_state = rMultinom( trans( p_des_state ) ); // needs a colvec
    sim_node_states(i,1) = des_state;
    if( des_node > n_tips ){ // If the des_node is not a tip node.
      ii = find( des_node == edge_mat.col(0) );
      sim_node_states(ii[0],0) = des_state;
      sim_node_states(ii[1],0) = des_state;
    }
	
    // With the start and end of the state we can perform the simulation. Forward simulate from the starting state, compute the changes along the branch, stop when the simulation passed the length of the branch. Then need to check if the arrival state is the same as the selected state for that node. Otherwise need to draw again.

    // Declare the dt vector:
    arma::vec dt = vec(n_states, fill::zeros);
    
    // Before simulating, check if any change is expected to happen at this node:
    double change_rate = -1.0 * Q(anc_state,anc_state);
      
    // bool equal_zero = approx_equal(change_rate, 0.0, "absdiff", 1.0e-12);
    if( change_rate < 0.0 || change_rate <= 1.0e-12 ){
      // Nothing will happen at this branch.
      // Time spent in the current state (anc_state) is the total branch length.
      dt[anc_state] = edge_len[i];
      // Double check if the ancestral state is the same of the descendant state:
      if( anc_state != sim_node_states(i,1) ){
	// Something is not good. Bad stochastic map.
	// Return empty matrix to signalize.
	arma::mat bad_maps = mat(size(mapped_edge), fill::zeros);
	return bad_maps;
      }
    } else{
      // Need to perform the simulation:
    
      // Reset the accept flag:
      bool loop_sims = true; // Set to keep simulating.
      // Set the tolerance for the edge length check:
      double edge_tol = 1.0e-12 * edge_len[i];
      // Track the number of trials. This will be used to break the simulation if necessary.
      int sims_trials = 0;
      
      while( loop_sims ){
	// The vector 'dt' tracks the 'dwelling time' on each state.
	// I am not creating the 'maps' history here.
	// Reset the dt vector:
	dt = vec(n_states, fill::zeros);
	// The starting state for the simulation.
	curr_state = anc_state;

	while( true ){
	  
	  // Here we record the number of times to make a valid simulation on the branch.
	  // We use the argument 'sims_limit' to break the stochastic map if the number of simulations
	  //    pass this limit. This particular proposal of stochastic map will be rejected by the MCMC.
	  sims_trials++;
	  if( sims_limit > 0 ){
	    // Trigger to check the limit.
	    if(sims_trials > sims_limit){
	      arma::mat bad_maps = mat(size(mapped_edge), fill::zeros);
	      return bad_maps;
	    }
	  }
	    
	  // Time until the next event:
	  // Note that this is 1/rate when compared with rexp in R.
	  time_chunk = R::rexp( 1/change_rate );
	  // Add to the time in the current state.
	  if( (sum(dt) + time_chunk) > (edge_len[i] - edge_tol) ){
	    // If the waiting time for the next event passes the edge length.
	    // Then add the remaining branch length and break.
	    time_chunk = edge_len[i] - sum(dt); // The rest of the time.
	    dt[curr_state] = dt[curr_state] + time_chunk;
	    break;
	  }
	  // Add to the time in the current state.
	  dt[curr_state] = dt[curr_state] + time_chunk;
	  // Updates the current state. Take a conditional sample based on Q:
	  sample_index = rMultinom( trans( trans_table_prob.row(curr_state) ) );
	  curr_state = trans_table_index(curr_state, sample_index);
	}
      
	// Check if the simulation is valid, conditioned on the state of the des node.
	if( curr_state == sim_node_states(i,1) ){
	  loop_sims = false;
	}
      }

    }

    // When done, store the result for this branch:
    mapped_edge.row(i) = trans(dt);
  }

  return mapped_edge;
}

// [[Rcpp::export]]
arma::mat makeSimmapMaps(arma::uword n_nodes, arma::uword n_tips, arma::uword n_states, arma::vec edge_len, arma::mat edge_mat, arma::vec parents, arma::mat X, arma::mat Q, int root_node, bool root_type, int max_nshifts) {
  // Same as the previous function. But this returns the 'maps' information that allow for reconstruction of the 'phytools' maps element.
  // Function works by assuming that the order of the rows in sim_node_states is the same as in recon_states below. This follows because of the code in 'getReconStates' function.
 
  // Define containers:
  int nrow_mapped_edge = edge_mat.n_rows;
  double time_chunk;
  arma::mat exp_Q = mat(Q);
  arma::uword anc_state;
  arma::uword des_state;
  arma::rowvec p_des_state;
  arma::uvec ii; // Index for many of the loops.
  arma::mat sim_node_states = mat(nrow_mapped_edge, 2);
  arma::vec prob_root = vec(n_states);
  arma::uword des_node;
  arma::uword curr_state;
  arma::uword sample_index;
  arma::umat trans_table_index = umat(n_states, n_states-1); // Take find output.
  arma::mat trans_table_prob = mat(n_states, n_states-1);
  arma::rowvec Q_row = rowvec(n_states);
  
  // Compute the probability for each state at each internal node given the model parameters. These are the conditional probabilities. Not the marginals.
  arma::mat recon_states = getReconStates(n_nodes, n_tips, n_states, edge_len, edge_mat, parents, X, Q, root_node, root_type);

  // Creates a table of transition probabilities given each starting state:
  // These only depend on Q and don't need to be computed more than once.
  for( arma::uword j=0; j < n_states; j++ ){
    arma::urowvec index_vec = urowvec(n_states);
    for( arma::uword jj=0; jj < n_states; jj++ ){
      index_vec[jj] = jj;
    }
    index_vec.shed_col(j);
    trans_table_index.row(j) = index_vec;
    // trans_table_index.row(j) = trans( find( Q_row > 0 ) );
    Q_row = Q.row(j); // This is a column vector type.
    Q_row.shed_col(j);
    trans_table_prob.row(j) = Q_row / sum( Q_row );
  }
  
  // The C++ code is not very efficient with growing vectors. So here we need to set a maximum number of per branch shifts events to work on the stochastic maps.
  // Rcout << "Max capacity for shifts on branches: " << max_nshifts << "\n";
  arma::umat maps_state = umat(nrow_mapped_edge, max_nshifts);
  maps_state.fill(n_states); // Help distinguish unnused entries.
  arma::mat maps_edge = mat(nrow_mapped_edge, max_nshifts, fill::zeros); // Zero entries are unnused.
  
  // Sample the state at the root.
  if( root_type ){ // 'madfitz'
    // Here the vector of probabilities is the same as the estimated under the model.
    arma::vec pi = trans( recon_states.row( root_node - 1 ) );
    prob_root = pi % pi;
  } else{ // root is 'equal'
    // Here the vector of probabilities is the same for all states.
    arma::vec pi = vec(n_states, fill::ones) / n_states;
    prob_root = trans( recon_states.row( root_node - 1 ) ) % pi;
  }
  anc_state = rMultinom( prob_root );
  ii = find( root_node == edge_mat.col(0) ); // Return length 2
  sim_node_states(ii[0],0) = anc_state;
  sim_node_states(ii[1],0) = anc_state;

  // Now we can sample all other states:
  // I will follow the reverse order of the post-order traversal edge matrix. This will make sure we are going from the root to the tips of the tree.
  for(int i=(nrow_mapped_edge-1); i >= 0; i-- ){
    
    exp_Q = expmat(Q * edge_len[i]);
    anc_state = sim_node_states(i,0); // This was already populated.
    des_node = edge_mat(i,1);
    p_des_state = exp_Q.row( anc_state ) % recon_states.row( des_node-1 ); // rowvec
    des_state = rMultinom( trans( p_des_state ) ); // needs a colvec
    sim_node_states(i,1) = des_state;
    if( des_node > n_tips ){ // If the des_node is not a tip node.
      ii = find( des_node == edge_mat.col(0) );
      sim_node_states(ii[0],0) = des_state;
      sim_node_states(ii[1],0) = des_state;
    }
	
    // With the start and end of the state we can perform the simulation. Forward simulate from the starting state, compute the changes along the branch, stop when the simulation passed the length of the branch. Then need to check if the arrival state is the same as the selected state for that node. Otherwise need to draw again.
    
    // Reset the accept flag:
    int accept = 0;
    while( accept == 0 ){
      // Reset the id for the shifts.
      uword id_shift = 0;
      // Reset the values for states and shift times.
      maps_state.row(i) = urowvec(max_nshifts, fill::ones) * n_states;
      maps_edge.row(i) = rowvec(max_nshifts, fill::zeros);
      // The starting state for the simulation is 'curr_state'
      curr_state = anc_state; // Because 'anc_state' is fixed and curr_state will change.
      maps_state(i, id_shift) = anc_state;

      while( true ){
	// Time until the next event:
	if( -1.0 * Q(curr_state,curr_state) <= 0 ){
	  // If the rate if 0, then nothing happens for the rest of the branch.
	  // This time_chunk will be larger than the branch. But the test below will make sure it breaks the loop and bounce it back to the correct length.
	  time_chunk = edge_len[i] + 1.0;
	} else{
	  time_chunk = R::rexp( 1/(-1.0 * Q(curr_state,curr_state)) );
	}
	// Add to the time in the current state.
	if( (sum( maps_edge.row(i) ) + time_chunk) > edge_len[i] ){
	  // If the waiting time for the next event passes the edge length.
	  // Then add the remaining branch length and break.
	  time_chunk = edge_len[i] - sum( maps_edge.row(i) ); // The rest of the time.
	  maps_edge(i,id_shift) = time_chunk;
	  break;
	}	   
	// Add to the time in the current state.
	maps_edge(i,id_shift) = time_chunk;
	// Updates the current state. Take a conditional sample based on Q:
	sample_index = rMultinom( trans( trans_table_prob.row( curr_state ) ) );
	curr_state = trans_table_index(curr_state, sample_index);
	
	// Check if the maximum number of events has been reached. Then return 0.
	if( id_shift+1 >= maps_state.n_cols ){
	  Rcout << "\n";
	  Rcout << "Max number of events on a branch reached! \n";
	  Rcout << "Please increase the value for the max_nshifts parameter. See Details. \n";
	  // Return a vector of 0.
	  arma::mat return_null = mat(2, 2, fill::zeros);
	  return return_null;
	}
	
	maps_state(i, id_shift+1) = curr_state; // Move to the next shift.
	id_shift++; // Advance the counter.
      }
      
      // Check if the simulation is valid, conditioned on the state of the des node.
      if( curr_state == sim_node_states(i,1) ){
	accept = 1;
      }
    }
  }

  // Need to consolidate the matrices in order to return.
  arma::mat maps_state_double = conv_to<mat>::from( maps_state ); // Force to double?
  arma::mat stack_maps = join_vert(maps_state_double, maps_edge);  
  return stack_maps;
}

// #######################################################
// ###### Prior distributions
// #######################################################

// Include prior densities.
#include <priors.h>

// #######################################################
// ###### The logLikelihood functions
// #######################################################

// [[Rcpp::export]]
double logLikMk_C(arma::uword n_nodes, arma::uword n_tips, arma::uword n_states, arma::vec edge_len, arma::mat edge_mat, arma::vec parents, arma::mat X, arma::mat Q, int root_node, int root_type) {
  // This is the log-lik function for a simple Mk model fitted to the tree. This will return the likelihood of the transition matrix for the rate regimes.

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

  // Get the log-lik for the model:
  return sum( log( comp.subvec(0 + n_tips, n_nodes + n_tips - 1) ) );
}

// Include log-likelihood for the mvBM model using pruning.
#include <likelihood_pruning.h>

// [[Rcpp::export]]
double logLikPrunningMCMC_C(arma::mat X, arma::uword k, arma::uword p, arma::vec nodes, arma::uvec des, arma::uvec anc, arma::uvec names_anc, arma::mat mapped_edge, arma::cube R, arma::vec mu) {
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

  // At the end, need to compute the contrast for the ROOT.
  // Make the calculation for the root log-likelihood.
  // The index 'n_nodes' is correspondent to the position 'n_nodes + 1'. Remember the indexation starts from 0.
  ss = X0.col(n_nodes-1) - mu;
  ll = ll + logLikNode_C(ss, V0.slice(n_nodes-1), inv(V0.slice(n_nodes-1)), k);

  return(ll);
}

// Define some proposal functions:
#include <proposals.h>

// #######################################################
// ###### The MCMC functions
// #######################################################

// Include functions to write to the log file:
#include <log_functions.h>

// This is the MCMC function when there is only a single tree/regime configuration to be estimated. Allowing for multiple trees could be done in the same function. However, it is faster and simpler just to duplicate this function and assume that we are always working with a list of trees.
// The function 'runRatematrixMultiMCMC_C' is doing this.
// There is a separated function also to deal with the case of joint inference of the Markov model for the predictor traits and the stochastic mapping for the tree.
// A similar need will happen for the MCMC with a single regime. In this case a bunch of computations are not needed. I can just simplify a lot this MCMC function. Using my own likelihood function will also means I can drop both 'mvMORPH' and 'phytools' as dependencies for the package. This will also make installation in a server much more user-friendly. [The 'rgl' dependecy of 'phytools' makes installation in servers difficult.]

// [[Rcpp::export]]
std::string runRatematrixMCMC_C(arma::mat X, int k, int p, arma::vec nodes, arma::uvec des, arma::uvec anc, arma::uvec names_anc, arma::mat mapped_edge, arma::cube R, arma::vec mu, arma::mat sd, arma::cube Rcorr, arma::vec w_mu, arma::mat par_prior_mu, std::string den_mu, arma::mat w_sd, arma::mat par_prior_sd, std::string den_sd, arma::vec nu, arma::cube sigma, arma::vec v, std::string log_file, std::string mcmc_file, double prob_sample_root, double prob_sample_sd, int gen, arma::vec post_seq, int write_header){
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
  // post_seq, a vector with the generations that will be written to file.

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
      curr_jacobian[j] = curr_jacobian[j] + ( log( var_vec(i,j) ) * log( (k-1.0)/2.0 ) );
    }
  }

  // A counter to help control when to write the sample to file.
  arma::uword post_seq_id = 0; // Keep track of the id for the gen seq to write.

  if( post_seq[post_seq_id] == 1 ){
    // Print starting point to files:
    log_stream << "1; 0; 0; 1; 1; ";
    log_stream << lik;
    log_stream << "\n"; 
    writeToMultFile_C(mcmc_stream, p, k, R, mu);
    post_seq_id++; // Updates the counter for the gen to write.
  }
  
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
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "1; 0; 0; 1; 1; ";
	  log_stream << prop_root_lik;
	  log_stream << "\n";
	}
	mu = prop_root; // Update the value for the root. Need to carry over.
	curr_root_prior = prop_root_prior; // Update the root prior. Need to carry over.
	lik = prop_root_lik; // Update likelihood. Need to carry over.
      } else{ // Reject. Keep the values the same.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "0; 0; 0; 1; 1; ";
	  log_stream << lik;
	  log_stream << "\n";
	}
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
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "1; 0; ";
	  log_stream << Rp+1; // Here is the regime.
	  log_stream << "; 0; 1; ";
	  log_stream << prop_sd_lik;
	  log_stream << "\n";
	}
	R = R_prop; // Update the evolutionary rate matrices. Need to carry over.
	sd = prop_sd; // Update the standard deviation.
	curr_sd_prior = prop_sd_prior;  // Update the prior. Need to carry over.
	lik = prop_sd_lik; // Update likelihood. Need to carry over.
      } else{ // Reject. Keep the values the same.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "0; 0; ";
	  log_stream << Rp+1; // Here is the regime.
	  log_stream << "; 0; 1; ";
	  log_stream << lik;
	  log_stream << "\n";
	}
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
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "1; ";
	  log_stream << Rp+1; // Here is the regime.
	  log_stream << "; 0; 0; 1; ";
	  log_stream << prop_corr_lik;
	  log_stream << "\n";
	}
	R = R_prop; // Update the evolutionary rate matrices. Need to carry over.
	Rcorr = Rcorr_prop; // Update the correlation matrix.
	Rcorr_curr_prior = Rcorr_prop_prior; // Update the prior. Need to carry over.
	lik = prop_corr_lik; // Update likelihood. Need to carry over.
	curr_jacobian[Rp] = prop_jacobian; // Updates jacobian.
      } else{ // Reject. Keep the values the same.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "0; ";
	  log_stream << Rp+1; // Here is the regime.
	  log_stream << "; 0; 0; 1; ";
	  log_stream << lik;
	  log_stream << "\n";
	}
      }
    }

    if( post_seq[post_seq_id] == (i+1) ){
      // Write the current state to the MCMC file.
      writeToMultFile_C(mcmc_stream, p, k, R, mu);
      // Finally, we also update the counter.
      post_seq_id++;
    }
    
  }

  Rcout << "Closing files... \n";

  mcmc_stream.close();
  log_stream.close();

  return "Done.";
}

// The MCMC function for multiple regimes:

// [[Rcpp::export]]
std::string runRatematrixMultiMCMC_C(arma::mat X, int k, int p, arma::mat nodes, arma::umat des, arma::umat anc, arma::umat names_anc, arma::cube mapped_edge, arma::cube R, arma::vec mu, arma::mat sd, arma::cube Rcorr, arma::vec w_mu, arma::mat par_prior_mu, std::string den_mu, arma::mat w_sd, arma::mat par_prior_sd, std::string den_sd, arma::vec nu, arma::cube sigma, arma::vec v, std::string log_file, std::string mcmc_file, double prob_sample_root, double prob_sample_sd, int gen, arma::vec post_seq, int write_header){
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
  // post_seq, vector that controls which of the generations are going to be saved to file.

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
      curr_jacobian[j] = curr_jacobian[j] + ( log( var_vec(i,j) ) * log( (k-1.0)/2.0 ) );
    }
  }

  // A counter to help control when to write the sample to file.
  arma::uword post_seq_id = 0; // Keep track of the id for the gen seq to write.
  
  if( post_seq[post_seq_id] == 1 ){
    // Print starting point to files:
    log_stream << "1; 0; 0; 1; 1; ";
    log_stream << lik;
    log_stream << "\n"; 
    writeToMultFile_C(mcmc_stream, p, k, R, mu);
    post_seq_id++;
  }

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
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "1; 0; 0; 1; ";
	  log_stream << prop_phy+1; // Sum back to the number of the tree.
	  log_stream << "; ";
	  log_stream << prop_root_lik;
	  log_stream << "\n";
	}
	mu = prop_root; // Update the value for the root. Need to carry over.
	curr_root_prior = prop_root_prior; // Update the root prior. Need to carry over.
	lik = prop_root_lik; // Update likelihood. Need to carry over.
	curr_phy = prop_phy; // Update the current tree in the pool.
      } else{ // Reject. Keep the values the same.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "1; 0; 0; 1; ";
	  log_stream << curr_phy+1; // Sum back to the number of the tree.
	  log_stream << "; ";
	  log_stream << lik;
	  log_stream << "\n";
	}
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
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "1; 0; ";
	  log_stream << Rp+1; // Here is the regime.
	  log_stream << "; 0; ";
	  log_stream << prop_phy+1; // Add 1 to get number of phy.
	  log_stream << "; ";
	  log_stream << prop_sd_lik;
	  log_stream << "\n";
	}
	R = R_prop; // Update the evolutionary rate matrices. Need to carry over.
	sd = prop_sd; // Update the standard deviation.
	curr_sd_prior = prop_sd_prior;  // Update the prior. Need to carry over.
	lik = prop_sd_lik; // Update likelihood. Need to carry over.
	curr_phy = prop_phy; // Update the phylogeny from the pool.
      } else{ // Reject. Keep the values the same.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "0; 0; ";
	  log_stream << Rp+1; // Here is the regime.
	  log_stream << "; 0; ";
	  log_stream << curr_phy+1; // Add 1 to get number of phy.
	  log_stream << "; ";
	  log_stream << lik;
	  log_stream << "\n";
	}
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
	if( post_seq[post_seq_id] == (i+1) ){
	  // This line will write to the mcmc_file.
	  // Instead of 'paste' I am using a line for each piece. Should have the same effect.
	  log_stream << "1; ";
	  log_stream << Rp+1; // Here is the regime.
	  log_stream << "; 0; 0; ";
	  log_stream << prop_phy+1; // Add 1 to get the number of the phylo.
	  log_stream << "; ";
	  log_stream << prop_corr_lik;
	  log_stream << "\n";
	}
	R = R_prop; // Update the evolutionary rate matrices. Need to carry over.
	Rcorr = Rcorr_prop; // Update the correlation matrix.
	Rcorr_curr_prior = Rcorr_prop_prior; // Update the prior. Need to carry over.
	lik = prop_corr_lik; // Update likelihood. Need to carry over.
	curr_jacobian[Rp] = prop_jacobian; // Updates jacobian.
	curr_phy = prop_phy; // Updates the current phy from the pool.
      } else{ // Reject. Keep the values the same.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "0; ";
	  log_stream << Rp+1; // Here is the regime.
	  log_stream << "; 0; 0; ";
	  log_stream << curr_phy+1; // Add 1 to get the number of the phylo.
	  log_stream << "; ";
	  log_stream << lik;
	  log_stream << "\n";
	}
      }
    }

    if( post_seq[post_seq_id] == (i+1) ){
      // Write the current state to the MCMC file.
      writeToMultFile_C(mcmc_stream, p, k, R, mu);
      // Update the index:
      post_seq_id++;
    }
    
  }

  Rcout << "Closing files... \n";

  mcmc_stream.close();
  log_stream.close();

  return "Done.";
}

// [[Rcpp::export]]
arma::mat buildQ(arma::vec vec_Q, arma::uword size, std::string model_Q){
  // Function to re-build the Q matrix.
  // February 2021: Function was exported to be able to use it to read in the MCMC.
  // Need to follow the same pattern used to extract the vector.
  arma::mat Q = mat(size, size, fill::zeros);
  
  if( model_Q == "ER" ){
    Q.fill(vec_Q[0]);
    // Now fill the diagonal.
    for( arma::uword i=0; i < size; i++ ){
      Q(i,i) = -1.0 * ( sum( Q.row(i) ) - vec_Q[0] );
    }
  } else if( model_Q == "SYM" ){
    Q.fill(0); // Fill the matrix with 0.
    arma::uword count = 0;
    // Go over the matrix and fill the upper and lower-tri.
    for( arma::uword i=0; i < size; i++ ){
      for( arma::uword j=0; j < size; j++ ){
	if( i >= j ) continue;
	Q(i,j) = vec_Q[count];
	Q(j,i) = vec_Q[count]; // The trick to fill the lower-tri.
	count++;
      }
    }
    // Now fill the diagonal.
    for( arma::uword i=0; i < size; i++ ){
      Q(i,i) = -1.0 * sum( Q.row(i) );
    }
  } else{ // model_Q == "ARD"
    Q.fill(0); // Fill the matrix with 0.
    arma::uword count = 0;
    // Go over the matrix and fill the upper and lower-tri.
    for( arma::uword i=0; i < size; i++ ){
      for( arma::uword j=0; j < size; j++ ){
	if( i == j ) continue;
	Q(i,j) = vec_Q[count];
	count++;
      }
    }
    // Now fill the diagonal.
    for( arma::uword i=0; i < size; i++ ){
      Q(i,i) = -1.0 * sum( Q.row(i) );
    }
  }

  return Q;
}

// [[Rcpp::export]]
std::string runRatematrixMCMC_jointMk_C(arma::mat X, arma::mat datMk, arma::uword k, arma::uword p, arma::vec nodes, arma::uword n_tips, arma::uvec des, arma::uvec anc, arma::uvec names_anc, arma::mat mapped_edge, arma::mat edge_mat, arma::uword n_nodes, arma::mat Q, double w_Q, std::string model_Q, int root_type, std::string den_Q, arma::vec par_prior_Q, arma::cube R, arma::vec mu, arma::mat sd, arma::cube Rcorr, arma::vec w_mu, arma::mat par_prior_mu, std::string den_mu, arma::mat w_sd, arma::mat par_prior_sd, std::string den_sd, arma::vec nu, arma::cube sigma, arma::vec v, std::string log_file, std::string mcmc_file, std::string Q_mcmc_file, arma::vec par_prob, arma::uword gen, arma::vec post_seq, int write_header, arma::uword sims_limit){

  // NOTE: 'sims_limit' is a parameter to reject the stochastic maps if it pass this limit.
  
  // The data parameters:
  // X, k, p, nodes, des, anc, names_anc, mapped_edge, datMk
  // The starting point parameters. These are the objects to carry on the MCMC.
  // R, mu, sd, Rcorr, Q, mapped_edge
  // Q is for the predictor regime and mapped_edge is the current stochastic mapping.
  // The starting point priors:
  // curr_root_prior, curr_sd_prior, Rcorr_curr_prior, curr_Q_prior
  // The parameters for the priors:
  // par_prior_mu, den_mu, par_prior_sd, den_sd, nu, sigma, par_prior_Q, den_Q
  // The starting point jacobian:
  // curr_jacobian
  // The parameters for the proposals:
  // w_mu, w_sd, v, w_Q
  // The parameters to control the MCMC:
  // prob_sample_root, prob_sample_var, prob_sample_Q, gen
  // write_header, wheather to write the header to file or just to append.

  // Open the files to write:
  // The log_file and mcmc_file arguments (one for the mvBM model and another for the Q matrix).
  std::ofstream log_stream (log_file, ios::out | ios::app);
  std::ofstream mcmc_stream (mcmc_file, ios::out | ios::app);
  std::ofstream Q_mcmc_stream (Q_mcmc_file, ios::out | ios::app);

  // Find the number of parameters for the Q matrix:
  arma::uword Q_npar;
  if( model_Q == "ER" ){
    Q_npar = 1;
  } else if( model_Q == "SYM" ){
    Q_npar = ( (p * p) - p ) / 2;
  } else {
    Q_npar = (p * p) - p;
  }
  
  // Write the header for the mcmc file.
  if(write_header == 1){
    for( arma::uword kk=1; kk < p+1; kk++ ){
      for( arma::uword ii=1; ii < k+1; ii++ ){
	for( arma::uword jj=1; jj < k+1; jj++ ){
	  mcmc_stream << "regime.p";
	  mcmc_stream << kk;
	  mcmc_stream << ".";
	  mcmc_stream << ii;
	  mcmc_stream << jj;
	  mcmc_stream << "; ";
	}
      }
    }
  
    for( arma::uword kk=1; kk < k; kk++ ){
      mcmc_stream << "trait.";
      mcmc_stream << kk;
      mcmc_stream << "; ";
    }
    mcmc_stream << "trait.";
    mcmc_stream << k;
    mcmc_stream << "\n";

    // Write the header for the Q mcmc file.
    for( arma::uword kk=1; kk < Q_npar; kk++ ){
      Q_mcmc_stream << "Q.par.";
      Q_mcmc_stream << kk;
      Q_mcmc_stream << "; ";
    }
    Q_mcmc_stream << "Q.par.";
    Q_mcmc_stream << Q_npar;
    Q_mcmc_stream << "\n";
    
    // Write the header for the log file.
    // Separate the log lik from the mvBM model from the one for the MK model.
  log_stream << "accepted; Q.matrix; stoch.map; smaps.limit; matrix.corr; matrix.sd; root; log.lik.BM; log.lik.MK \n";
  } else{
    // Do nothing.
    // This is the case for a continuing MCMC.
  }
  
  // Define the containers for the run:
  // Try to set the size of the objects for better memory management.
  
  // The starting point log lik for each of the models:
  double lik_mvBM;
  double lik_Mk;
  // The starting point prior:
  double curr_root_prior;
  double curr_sd_prior;
  double Rcorr_curr_prior;
  double curr_Q_prior;
  
  // The jacobian for the MCMC. There are one for each regime.
  // Because need to track the jump separatelly.
  arma::vec curr_jacobian = vec(k, fill::zeros);

  int sample_par;
  arma::vec prop_root;
  double prop_root_prior;
  double pp;
  double prop_root_lik;
  double ll;
  double r;
  double unif_draw;
  arma::uword Rp;
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
  arma::mat prop_Q = mat(Q);
  double prop_Q_prior;
  double prop_Q_lik;
  arma::vec vec_Q;
  arma::vec prop_vec_Q;
  arma::vec multi_Q_factor;
  arma::mat prop_mapped_edge = mat(mapped_edge);
  double prop_mapped_edge_lik;
  
  
  // Generate some additional information:

  // Set the vector of probabilities for the regimes:
  arma::vec regime_prob = vec(p);
  regime_prob.fill(1.0/p); // Same probability to sample any regime.
  
  // Branch lengths is the sum of the mapped_edges
  arma::vec edge_len = vec(mapped_edge.n_rows);
  for( arma::uword i=0; i < mapped_edge.n_rows; i++ ){
    edge_len[i] = sum( mapped_edge.row(i) );
  }
  
  int root_node = n_tips + 1; // The root node for the phylogeny.

  // Get starting priors, likelihood, and jacobian.
  // The likelihood for the model is the sum of the logliks of these two layers.
  lik_mvBM = logLikPrunningMCMC_C(X, k, p, nodes, des, anc, names_anc, mapped_edge, R, mu);
  lik_Mk = logLikMk_C(n_nodes, n_tips, p, edge_len, edge_mat, nodes, datMk, Q, root_node, root_type);

  curr_root_prior = priorRoot_C(mu, par_prior_mu, den_mu);
  curr_sd_prior = priorSD_C(sd, par_prior_sd, den_sd);
  Rcorr_curr_prior = priorCorr_C(Rcorr, nu, sigma);

  // The vectorized Q matrix for the starting state of the search.
  // Size of vector will depend on the model used for Q.
  vec_Q = extractQ(Q, p, model_Q);
  // Create vector with the multiplier factor for the Q matrix.
  // This will control the multiplier proposals.
  arma::vec w_Q_vec = vec(vec_Q);
  w_Q_vec.fill(w_Q); // Fill with the multiplier factor.
  
  curr_Q_prior = priorQ(vec_Q, par_prior_Q, den_Q); // The prior for the Q matrix.

  Rcout << "Starting point Log-likelihood: " << lik_mvBM + lik_Mk << "\n";

  arma::mat var_vec = square(sd);
  // Jacobian for both the regimes:
  for( arma::uword j=0; j < p; j++ ){
    for( arma::uword i=0; i < k; i++ ){
      // The jacobian is computed on the variances!
      curr_jacobian[j] = curr_jacobian[j] + ( log( var_vec(i,j) ) * log( (k-1.0)/2.0 ) );
    }
  }

  // A counter to help control when to write the sample to file.
  arma::uword post_seq_id = 0; // Keep track of the id for the gen seq to write.

  if( post_seq[post_seq_id] == 1 ){
    // Print starting point to files:
    log_stream << "1; 0; 0; 0; 0; 0; 0; ";
    log_stream << lik_mvBM;
    log_stream << ";";
    log_stream << lik_Mk;
    log_stream << "\n"; 
    writeToMultFile_C(mcmc_stream, p, k, R, mu);
    // Need to make a function to write the Q matrix to file.
    // Note that the length of the vector will change depending on k and model_Q
    writeQToFile(Q_mcmc_stream, vec_Q, p, model_Q);
    post_seq_id++;
  }

  Rcout << "Starting MCMC ... \n";
  
  // Starting the MCMC.
  for( arma::uword i=0; i < gen; i++ ){
    // Sample between the root, rate matrix, and the Mk matrix.
    // The result will be in the form 0: root, 1: sd, 2: corr, 3: Q, 4: stochastic map.
    // So a single vector control the probability to sample all the parameters.
    sample_par = rMultinom(par_prob); // par_prob needs to have length of 5 and sum to 1.
    // If matrix or variance is updated, which regime?
    Rp = rMultinom(regime_prob);
    // Rp is an arma::uword index.
	  
    // The 'if...if...else' structure is lazy and will evaluate only the first entry.
    if(sample_par == 0){
      // Update the root state.
      prop_root = slideWindow_C(mu, w_mu);
      // Compute the prior.
      prop_root_prior = priorRoot_C(prop_root, par_prior_mu, den_mu);
      // The log prior ratio. curr_root_prior is the current prior density for the root.
      // The prior for the other parameters are also constant. So they don't need to be here either.
      pp = prop_root_prior - curr_root_prior;
      prop_root_lik = logLikPrunningMCMC_C(X, k, p, nodes, des, anc, names_anc, mapped_edge, R, prop_root);
      // The likelihood ratio.
      // Because the lik for the mk model is constant here, it does not need to be computed.
      ll = prop_root_lik - lik_mvBM;
      // Get the ratio in log space.
      r = ll + pp;

      // Advance to the acceptance step.
      // Here we are only updating the root, so all other parameters are the same.
      unif_draw = as_scalar(randu(1)); // The draw from a uniform distribution.
      if( exp(r) > unif_draw ){ // Accept.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "1; 0; 0; 0; 0; 0; 1; ";
	  log_stream << prop_root_lik;
	  log_stream << ";";
	  log_stream << lik_Mk;
	  log_stream << "\n";
	}
	mu = prop_root; // Update the value for the root. Need to carry over.
	curr_root_prior = prop_root_prior; // Update the root prior. Need to carry over.
	lik_mvBM = prop_root_lik; // Update likelihood. Need to carry over.
      } else{ // Reject. Keep the values the same.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "0; 0; 0; 0; 0; 0; 1; ";
	  log_stream << lik_mvBM;
	  log_stream << ";";
	  log_stream << lik_Mk;
	  log_stream << "\n";
	}
      }
    } else if(sample_par == 1){
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
      ll = prop_sd_lik - lik_mvBM;

      // Get the ratio i log space. Loglik, log prior and the proposal ratio (for the multiplier!).
      r = ll + pp + accu(multi_factor);

      // Advance to the acceptance step.
      // Here we are only updating the root, so all other parameters are the same.
      unif_draw = as_scalar(randu(1)); // The draw from a uniform distribution.
      if( exp(r) > unif_draw ){ // Accept.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "1; 0; 0; 0; 0; ";
	  log_stream << Rp+1; // Here is the regime.
	  log_stream << "; 0; ";
	  log_stream << prop_sd_lik;
	  log_stream << ";";
	  log_stream << lik_Mk;
	  log_stream << "\n";
	}
	R = R_prop; // Update the evolutionary rate matrices. Need to carry over.
	sd = prop_sd; // Update the standard deviation.
	curr_sd_prior = prop_sd_prior;  // Update the prior. Need to carry over.
	lik_mvBM = prop_sd_lik; // Update likelihood. Need to carry over.
      } else{ // Reject. Keep the values the same.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "0; 0; 0; 0; 0; ";
	  log_stream << Rp+1; // Here is the regime.
	  log_stream << "; 0; ";
	  log_stream << lik_mvBM;
	  log_stream << ";";
	  log_stream << lik_Mk;
	  log_stream << "\n";
	}
      }
      
    } else if(sample_par == 2){
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
      ll = prop_corr_lik - lik_mvBM;
      // This is the sdiance of the proposed vcv matrix.
      decomp_var = diagvec( Rcorr_prop.slice(Rp) ); // These are variances
      prop_jacobian = 0.0; // Always need to reset.
      // The Jacobian of the transformation.
      for( arma::uword i=0; i < k; i++ ){
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
	if( post_seq[post_seq_id] == (i+1) ){
	  // This line will write to the mcmc_file.
	  // Instead of 'paste' I am using a line for each piece. Should have the same effect.
	  log_stream << "1; 0; 0; 0; ";
	  log_stream << Rp+1; // Here is the regime.
	  log_stream << "; 0; 0; ";
	  log_stream << prop_corr_lik;
	  log_stream << ";";
	  log_stream << lik_Mk;
	  log_stream << "\n";
	}
	R = R_prop; // Update the evolutionary rate matrices. Need to carry over.
	Rcorr = Rcorr_prop; // Update the correlation matrix.
	Rcorr_curr_prior = Rcorr_prop_prior; // Update the prior. Need to carry over.
	lik_mvBM = prop_corr_lik; // Update likelihood. Need to carry over.
	curr_jacobian[Rp] = prop_jacobian; // Updates jacobian.
      } else{ // Reject. Keep the values the same.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "0; 0; 0; 0; ";
	  log_stream << Rp+1; // Here is the regime.
	  log_stream << "; 0; 0; ";
	  log_stream << lik_mvBM;
	  log_stream << ";";
	  log_stream << lik_Mk;
	  log_stream << "\n";
	}
      }
    } else if(sample_par == 3){
      // Update the Q matrix, then make a new stochastic map.
      // When updating the Q matrix we need also to make a new draw for the stochastic map. Otherwise, the change in the Q matrix will not match the likelihood of the mvBM model. If I don't make a new stochastic map draw we could get a ugly disconnection between the Q matrix and the stochastic maps that the mvBM models are evaluated under.

      // UPDATE Q MATRIX
      // The current 'vec_Q' was already defined in the beginning of the MCMC function.
      // vec_Q = extractQ(Q, model_Q); // Get a vector with the parameters that describe the Q matrix.
      // Next line assumes 'vec_Q' is a column vector.
      multi_Q_factor = multiplierProposal_C(vec_Q.n_rows, w_Q_vec);
      prop_vec_Q = vec_Q % multi_Q_factor;
      prop_Q_prior = priorQ(prop_vec_Q, par_prior_Q, den_Q);
      
      // The prior only reflects on the Q matrix. The prior for the stochastic map (conditioned on the Q matrix) is flat.
      pp = prop_Q_prior - curr_Q_prior;
      prop_Q = buildQ(prop_vec_Q, p, model_Q); // Rebuild a matrix from the Q vector. Just to compute the lik and draw the map.

      // UPDATE MAPPED_EDGE
      // Need to check if the stochastic map is valid.
      // If the returned mapped matrix has accu of 0, then reject this move.
      prop_mapped_edge = makeSimmapMappedEdge(n_nodes, n_tips, p, edge_len, edge_mat, nodes, datMk, prop_Q, root_node, root_type, sims_limit);
      if( accu( prop_mapped_edge ) < max(edge_len) ){
	if( post_seq[post_seq_id] == (i+1) ){
	  // The mapped_edge returned an invalid matrix.
	  // Reject, mark the 'smaps.limit' column.
	  log_stream << "0; 1; 1; 1; 0; 0; 0; ";
	  log_stream << lik_mvBM;
	  log_stream << ";";
	  log_stream << lik_Mk;
	  log_stream << "\n";
	  // Write the current generation to files (note the continue flag after).
	  writeToMultFile_C(mcmc_stream, p, k, R, mu);
	  writeQToFile(Q_mcmc_stream, vec_Q, p, model_Q);
	}
	// Break generation.
	continue;
      }

      // Compute the liklihood for the Q matrix
      prop_Q_lik = logLikMk_C(n_nodes, n_tips, p, edge_len, edge_mat, nodes, datMk, prop_Q, root_node, root_type);
      
      // This the mvBM likelihood, just changed the mapped_edge.
      prop_mapped_edge_lik = logLikPrunningMCMC_C(X, k, p, nodes, des, anc, names_anc, prop_mapped_edge, R, mu);

      // Here the likelihood ratio need to take into account the whole model.
      ll = (prop_Q_lik + prop_mapped_edge_lik) - (lik_Mk + lik_mvBM);
      
      // Get the ratio i log space. Loglik, log prior and the proposal ratio (for the multiplier!).
      r = ll + pp + accu(multi_Q_factor);

      // Advance to the acceptance step.
      // Here we are only updating the root, so all other parameters are the same.
      unif_draw = as_scalar(randu(1)); // The draw from a uniform distribution.
      if( exp(r) > unif_draw ){ // Accept.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "1; 1; 1; 0; 0; 0; 0; ";
	  log_stream << prop_mapped_edge_lik;
	  log_stream << ";";
	  log_stream << prop_Q_lik;
	  log_stream << "\n";
	}
	curr_Q_prior = prop_Q_prior;  // Update the prior. Need to carry over.
	vec_Q = prop_vec_Q; // Passing the vectorized Q matrix.
	Q = prop_Q; // Passing the Q matrix.
	mapped_edge = prop_mapped_edge; // Update the stochastic map.
	lik_Mk = prop_Q_lik;
	lik_mvBM = prop_mapped_edge_lik; // Update likelihood. Need to carry over.
      } else{ // Reject. Keep the values the same.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "0; 1; 1; 0; 0; 0; 0; ";
	  log_stream << lik_mvBM;
	  log_stream << ";";
	  log_stream << lik_Mk;
	  log_stream << "\n";
	}
      }
      
    } else{
      // Update the ONLY the stochastic map.

      // UPDATE MAPPED_EDGE
      // Need to check if the stochastic map is valid.
      // If the returned mapped matrix has accu of 0, then reject this move.
      prop_mapped_edge = makeSimmapMappedEdge(n_nodes, n_tips, p, edge_len, edge_mat, nodes, datMk, prop_Q, root_node, root_type, sims_limit);
      if( accu( prop_mapped_edge ) < max(edge_len) ){
	if( post_seq[post_seq_id] == (i+1) ){
	  // The mapped_edge returned an invalid matrix.
	  // Reject, mark the 'smaps.limit' column.
	  log_stream << "0; 1; 1; 1; 0; 0; 0; ";
	  log_stream << lik_mvBM;
	  log_stream << ";";
	  log_stream << lik_Mk;
	  log_stream << "\n";
	  // Write the current generation to files (note the continue flag after).
	  writeToMultFile_C(mcmc_stream, p, k, R, mu);
	  writeQToFile(Q_mcmc_stream, vec_Q, p, model_Q);
	}
	// Break generation.
	continue;
      }
      
      // This the mvBM likelihood, just changed the mapped_edge.
      prop_mapped_edge_lik = logLikPrunningMCMC_C(X, k, p, nodes, des, anc, names_anc, prop_mapped_edge, R, mu);

      // Here the likelihood ratio is only of the mvBM model. Q does not change.
      // The prior ratio is log(1) . Because there is no prior for a draw of the stochastic map.
      // Thus 'r = ll;' and 'll = prop_mapped_edge_lik - lik_mvBM;'
      r = prop_mapped_edge_lik - lik_mvBM;

      // Advance to the acceptance step.
      // Here we are only updating the root, so all other parameters are the same.
      unif_draw = as_scalar(randu(1)); // The draw from a uniform distribution.
      if( exp(r) > unif_draw ){ // Accept.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "1; 0; 1; 0; 0; 0; 0; ";
	  log_stream << prop_mapped_edge_lik;
	  log_stream << ";";
	  log_stream << lik_Mk;
	  log_stream << "\n";
	}
	mapped_edge = prop_mapped_edge; // Update the stochastic map.
	lik_mvBM = prop_mapped_edge_lik; // Update likelihood. Need to carry over.
      } else{ // Reject. Keep the values the same.
	if( post_seq[post_seq_id] == (i+1) ){
	  log_stream << "0; 0; 1; 0; 0; 0; 0; ";
	  log_stream << lik_mvBM;
	  log_stream << ";";
	  log_stream << lik_Mk;
	  log_stream << "\n";
	}
      }
    }

    if( post_seq[post_seq_id] == (i+1) ){
      // Write the current state to the MCMC file.
      // Need to use a different output file to write the Q matrix.
      writeToMultFile_C(mcmc_stream, p, k, R, mu);
      // Q here is a vector, so we need size of matrix and model to write the correct quantity.
      writeQToFile(Q_mcmc_stream, vec_Q, p, model_Q);
      // Update the counter for writting to file.
      post_seq_id++;
    }
    
  }

  Rcout << "Closing files... \n";
  
  mcmc_stream.close();
  log_stream.close();
  Q_mcmc_stream.close();

  return "Done.";
}

