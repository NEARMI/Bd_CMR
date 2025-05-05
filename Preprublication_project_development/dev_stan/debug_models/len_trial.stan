functions {
// ------------------------------ functions ------------------------------

// -----
// chi_t probability estimator
// -----

// Function to build the required chi_t probability estimator for the mark recapture model, which
 // is the probability of never recapturing an individual again after capturing them at time t
  // rewritten to run for an individual at a time given variable numbers of capture opportunities by individual (e.g. in different populations)
	
	real[] prob_uncaptured(int n_occ, vector p_sub, vector phi_sub) {

	real chi_sub[n_occ];          // chi for each capture date and individual 

        chi_sub[n_occ] = 1.0;         // on the last sampling date the probability is one

	for (t in 1:(n_occ - 1)) {    // loop over from the first to the second to last sampling date and multiply out the probabilities
         int t_curr = n_occ - t;
         int t_next = t_curr + 1;
	 chi_sub[t_curr] = (1 - phi_sub[t_curr]) + phi_sub[t_curr] * (1 - p_sub[t_next]) * chi_sub[t_next];		
      }

  return chi_sub;

 }

}


data {
// ------------------------------ data ------------------------------

// -----
// Notes
// -----
// Given long comments, this model is best read in full screen on an external monitor


  // dimensional and bookkeeping params (single vals)
	int<lower=1> n_pop;				    // Number of distinct populations (sampling areas)
	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years and all populations)
	int<lower=1> n_spec;				    // Total number of species
	int<lower=1> n_sex;				    // Total number of entries given for Sex (should be 3 if data cleaning worked [F, M, U]
	int<lower=1> n_col_mm;

	matrix[n_ind, n_sex] ind_sex;		  	    // Sex of each individual
	matrix[n_ind, n_spec] ind_spec;			    // The species of each individual
	matrix[n_ind, n_col_mm] ind_mm;

	int<lower=1> ind_in_pop[n_ind];		   	    // Which population each individual belongs to

  // covariates (length)
	int<lower=1> n_ind_len_have;			    // Number of individuals that we have length data	  
	int<lower=1> n_ind_len_mis;			    // Number of individuals with missing length data
	int<lower=1> ind_len_which_have[n_ind_len_have];    // Index of individuals that we have length data
	int<lower=1> ind_len_which_mis[n_ind_len_mis];      // Index of individuals with missing length data
	vector[n_ind_len_have] ind_len_have;		    // The actual length values that we have
 	int<lower=1> ind_len_spec_first_index[n_spec]; 	    // First size index associated with each unique species (for species-specific scaling of length values)
 	int<lower=1> ind_len_spec_size[n_spec];      	    // Number of individuals of each species with lengths (for species-specific scaling of length values)

}

parameters {
// ------------------------------ parameters ------------------------------

// -----
// imputed covariates: length
// -----

	vector[n_col_mm] beta_len;
	vector[n_ind_len_mis] ind_len_mis;		 // the imputed values of len
  	real<lower=0> sd_len[n_pop]; 	

}


transformed parameters {
// ------------------------------ transformed parameters ------------------------------

  // Individual lengths 
  	vector[n_ind_len_have] mu_len_have; 		  	 
  	vector[n_ind_len_mis] mu_len_mis; 
	vector[n_ind] ind_len;				 // all individual len (combining data and imputed values)
	vector[n_ind] ind_len_scaled;			 // all individual len scaled

// -----
// Imputed NA Length values
// -----

  // linear predictor for len regression on measured lengths
	mu_len_have   = ind_mm[ind_len_which_have, ] * beta_len;
 
  // linear predictor for len regression for imputing unknown lengths
  	mu_len_mis    = ind_mm[ind_len_which_mis, ] * beta_len;   

  // filling in the complete vector of ind_mehg with the data and missing values
	ind_len[ind_len_which_have] = ind_len_have;
	ind_len[ind_len_which_mis]  = ind_len_mis;

  // Scaling the predicted lengths within-species (loop over species)
	for (ns in 1:n_spec) {

  // Jump through a hoop to select out all of the length values for a given species 
	  vector[ind_len_spec_size[ns]] temp_ind_len = segment(ind_len, ind_len_spec_first_index[ns], ind_len_spec_size[ns]);

  // Scale the lengths of species ns and stick them in the complete long-form container
	  ind_len_scaled[ind_len_spec_first_index[ns]:(ind_len_spec_first_index[ns] + ind_len_spec_size[ns] - 1)] = (temp_ind_len - mean(temp_ind_len))/sd(temp_ind_len);

	}

}

model {
// ------------------------------ model ------------------------------

// -----
// Priors
// -----

  // Imputed Covariates Priors: length

	sd_len    ~ inv_gamma(8, 15);	
	beta_len  ~ normal(0, 3);


// -----
// Imputed NA Lengths
// -----

	ind_len_have ~ normal(mu_len_have, sd_len[ind_in_pop[ind_len_which_have]]);
	ind_len_mis  ~ normal(mu_len_mis, sd_len[ind_in_pop[ind_len_which_mis]]);

}

