// Define multiple functions to log posterior samples to files.

#ifndef log_fn_h
#define log_fn_h

void writeToMultFile_C(std::ostream& mcmc_stream, arma::uword p, arma::uword k, arma::cube R, arma::vec mu){
  // Note the 'std::ostream&' argument here is the use of a reference.
  // This function will work fine when p = 1. No need to create exception for the case of a single regime.
  for( arma::uword i=0; i < p; i++ ){
    for( arma::uword j=0; j < k; j++ ){
      for( arma::uword z=0; z < k; z++ ){
	mcmc_stream << R.slice(i)(j,z);
	mcmc_stream << "; ";
      }
    }
  }
  for( arma::uword i=0; i < k-1; i++ ){
    mcmc_stream << mu[i];
    mcmc_stream << "; ";
  }
  mcmc_stream << mu.tail(1);
  // mcmc_stream << "\n";

}

void writeQToFile(std::ostream& Q_mcmc_stream, arma::vec vec_Q, arma::uword k, std::string model_Q){
  // Note the 'std::ostream&' argument here is the use of a reference.
  if( model_Q == "ER" ){
    Q_mcmc_stream << vec_Q;
  } else{
    arma::uword print_size = vec_Q.n_rows;
    for( arma::uword i=0; i < (print_size-1); i++ ){
      Q_mcmc_stream << vec_Q[i];
      Q_mcmc_stream << "; ";
    }
    Q_mcmc_stream << vec_Q[print_size-1];
    Q_mcmc_stream << "\n";
  }
}

// Some functions to work with polytopes.

void writeToMultFileNoRoot(std::ostream& mcmc_stream, arma::uword p, arma::uword k, arma::cube R){
  // Same as the 'writeToMultFile_C' function but without the root value.
  // Note the 'std::ostream&' argument here is the use of a reference.

  // In the case of a single regime we need to make a fix to the workaround below.
  if( p == 1 ){

    for( arma::uword j=0; j < k; j++ ){
      for( arma::uword z=0; z < k; z++ ){
	mcmc_stream << R.slice(0)(j,z);
	if( j == k-1 && z == k-1 ) {
	  // Finish the line on the last element.
	  mcmc_stream << "\n";
	} else{
	  // Not the last, add a separator.
	  mcmc_stream << "; ";
	}
      }
    }
    
  } else{ // Normal case with multiple regimes fitted to the tree.
    // Will need a work-around to be able to close the line.
    for( arma::uword i=0; i < p-1; i++ ){
      for( arma::uword j=0; j < k; j++ ){
	for( arma::uword z=0; z < k; z++ ){
	  mcmc_stream << R.slice(i)(j,z);
	  mcmc_stream << "; ";
	}
      }
    }

    for( arma::uword j=0; j < k; j++ ){
      for( arma::uword z=0; z < k; z++ ){
	mcmc_stream << R.slice(p-1)(j,z);
	if( j == k-1 && z == k-1 ) {
	  // Finish the line on the last element.
	  mcmc_stream << "\n";
	} else{
	  // Not the last, add a separator.
	  mcmc_stream << "; ";
	}
      }
    }
  }
  
}

void writePolySample(std::ostream& poly_stream, arma::mat poly_tips, arma::mat poly_nodes){
  // Write the sample for the tip states and for the internal nodes to file.

  // Need to write down by row.
  for( arma::uword i=0; i < poly_tips.n_rows; i++ ) {
    for( arma::uword j=0; j < poly_tips.n_cols; j++ ) {
      poly_stream << poly_tips(i,j);
      poly_stream << "; ";
    }
  }

  for( arma::uword i=0; i < poly_nodes.n_rows; i++ ) {
    for( arma::uword j=0; j < poly_nodes.n_cols; j++ ) {
      poly_stream << poly_nodes(i,j);
      if( i == (poly_nodes.n_rows-1) && j == (poly_nodes.n_cols-1) ) {
	// Last element, add a line end.
	poly_stream << "\n";
      } else{
	// Not the last element. Add a separator.
	poly_stream << "; ";
      }
    }
  }

}

void writePolySampleTipsOnly(std::ostream& poly_stream, arma::mat poly_tips){
  
  // Write the sample for the tips to the log file.
  for( arma::uword i=0; i < poly_tips.n_rows; i++ ) {
    for( arma::uword j=0; j < poly_tips.n_cols; j++ ) {
      poly_stream << poly_tips(i,j);
      if( i == (poly_tips.n_rows-1) && j == (poly_tips.n_cols-1) ){
	// Last element, add a line end.
	poly_stream << "\n";
      } else{
	poly_stream << "; ";
      }
    }
  }

}

// Write the samples for the Gamma rates to file.
void writeGammaSample(std::ostream& gamma_stream, double beta, arma::vec branches){
  
  // Write the sample for the tips to the log file.
  gamma_stream << beta;
  gamma_stream << "; ";
  for( arma::uword i=0; i < branches.n_rows-1; i++ ) {
      gamma_stream << branches(i);
      gamma_stream << "; ";
  }
  gamma_stream << branches( branches.n_rows-1 );
  gamma_stream << "\n";

}

void writeSDMat(std::ostream& sd_stream, arma::mat sd_mat){

  arma::uword nrows = sd_mat.n_rows;
  arma::uword ncols = sd_mat.n_cols;
  
  for( arma::uword i=0; i < nrows; i++ ){
    for( arma::uword j=0; j < ncols; j++ ){
      sd_stream << sd_mat(i,j);
      if(i == nrows-1 && j == ncols-1){
	// Last element. Write the end of the line.
	sd_stream << "\n";
      } else{
	sd_stream << "; ";
      }
    }
  }

}

#endif

// Local Variables:
// mode: c++
// End:
