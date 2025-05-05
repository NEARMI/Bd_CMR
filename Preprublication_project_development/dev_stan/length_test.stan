functions {
// ------------------------------ functions ------------------------------

}


data {
// ------------------------------ data ------------------------------

// -----
// Notes
// -----
// Given long comments, this model is best read in full screen on an external monitor


  // dimensional and bookkeeping params (single vals)
	int<lower=1> n_pop;				    // Number of distinct populations (sampling areas)
	int<lower=1> n_pop_year;			    // Index for pop*year
	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years and all populations)
	int<lower=1> ind_per_period_bd;	                    // n_ind * year (the unique periods of time in which each individual's bd is estimated)
	int<lower=0> ind_occ;			   	    // n_ind * all sampling periods (all events in which each individual could potentially have been captured)
	int<lower=0> ind_occ_min1;		 	    // n_ind * all sampling periods except the last 
	int<lower=0> n_days;				    // Number of sampling occasions
	int<lower=0> n_spec;				    // Total number of species
	int<lower=0> n_sex;				    // Total number of entries given for Sex (should be 3 if data cleaning worked [F, M, U]
	
  // short vector indexes (length of n_ind)	
	int<lower=1> ind_occ_size[n_ind];		    // Number of sampling periods for all individuals
	int<lower=1> ind_occ_min1_size[n_ind];		    // Number of sampling periods -1 for all individuals
	int<lower=1> phi_first_index[n_ind];		    // The indexes of phi corresponding to the first entry for each individual
	int<lower=1> p_first_index[n_ind];	            // The indexes of p corresponding to the first entry for each individual
	int<lower=0> ind_in_pop[n_ind];		   	    // Which population each individual belongs to
	int<lower=0> ind_sex[n_ind];			    // The sex of each individual
	int<lower=0> ind_spec[n_ind];			    // The species of each individual

  // covariates (length)
	int<lower=0> n_ind_len_have;			    // Number of individuals that we have length data	  
	int<lower=0> n_ind_len_mis;			    // Number of individuals with missing length data
	int<lower=0> ind_len_which_have[n_ind_len_have];    // Index of individuals that we have length data
	int<lower=0> ind_len_which_mis[n_ind_len_mis];      // Index of individuals with missing length data
	vector[n_ind_len_have] ind_len_have;		    // The actual length values that we have
 	int<lower=0> ind_len_spec_first_index[n_spec]; 	    // First size index associated with each unique species (for species-specific scaling of length values)
 	int<lower=0> ind_len_spec_size[n_spec];      	    // Number of individuals of each species with lengths (for species-specific scaling of length values)

	matrix[n_ind_len_have, n_spec] ind_len_spec_have;    
	matrix[n_ind_len_mis, n_spec] ind_len_spec_mis;	     
	matrix[n_ind_len_have, n_sex] ind_len_sex_have;	    
	matrix[n_ind_len_mis, n_sex] ind_len_sex_mis;	    

}


transformed data {

} 


parameters {
// ------------------------------ parameters ------------------------------

// -----
// imputed covariates: length
// -----

	real<lower=0> inverse_phi_len;		         // variance parameter for gamma regression
	vector[n_sex] beta_len_sex;			 // regression coefficient for len as a function of sex
	vector[n_spec] beta_len_spec;			 // regression coefficient for len as a function of species

	vector[n_ind_len_mis] ind_len_mis;		 // the imputed values of len

}


transformed parameters {
// ------------------------------ transformed parameters ------------------------------

  // Individual lengths 
	
  	vector[n_ind_len_have] mu_len_have; 		 // the expected values for the gamma regression
  	vector[n_ind_len_have] rate_len_have; 	 	 // rate parameter for the gamma distribution

  	vector[n_ind_len_mis] mu_len_mis; 		 // the expected values (linear predictor) for the missing len values
  	vector[n_ind_len_mis] rate_len_mis; 		 // rate parameter for the gamma distribution for the missing len values

	vector[n_ind] ind_len;				 // all individual len (combining data and imputed values)
	vector[n_ind] ind_len_scaled;			 // all individual len scaled

	vector[n_spec] ind_len_mean;			 // mean of ind_len
	vector[n_spec] ind_len_sd;			 // sd of ind_len
 
// -----
// Imputed NA Length values
// -----

  // linear predictor for len regression on measured lengths

  	mu_len_have   = exp(ind_len_sex_have * beta_len_sex + ind_len_spec_have * beta_len_spec);  
  	rate_len_have = rep_vector(inverse_phi_len, n_ind_len_have) ./ mu_len_have;

  // linear predictor for len regression for imputing unknown lengths
	mu_len_mis   = exp(ind_len_spec_mis * beta_len_spec + ind_len_sex_mis * beta_len_sex);  
  	rate_len_mis = rep_vector(inverse_phi_len, n_ind_len_mis) ./ mu_len_mis;

	ind_len[ind_len_which_have] = ind_len_have;	 // filling in the complete vector of ind_mehg with the data
	ind_len[ind_len_which_mis]  = ind_len_mis;       // filling in the complete vector of ind_mehg with the imputed values

  // Scaling the predicted lengths within-species
	for (ns in 1:n_spec) {

  // Jump through a hoop to select out all of the length values for a given species 
	  vector[ind_len_spec_size[ns]] temp_ind_len = segment(ind_len, ind_len_spec_first_index[ns], ind_len_spec_size[ns]);

	  ind_len_mean[ns] = mean(temp_ind_len);	 // take the mean of the lengths of all individuals of species ns
	  ind_len_sd[ns]   = sd(temp_ind_len);		 // take the sd of the lengths of all individuals of species ns

  // Scale the lengths of species ns and stick them in the complete long-form container
	  ind_len_scaled[ind_len_spec_first_index[ns]:(ind_len_spec_first_index[ns] + ind_len_spec_size[ns] - 1)] = (temp_ind_len - ind_len_mean[ns])/ind_len_sd[ns];

	}

}

model {
// ------------------------------ model ------------------------------

// -----
// Priors
// -----

  // Imputed Covariates Priors: length

	inverse_phi_len  ~ inv_gamma(8, 15);	
	beta_len_sex     ~ normal(0, 3);
	beta_len_spec    ~ normal(0, 3);

// -----
// Imputed NA Lengths
// -----

	ind_len_have ~ gamma(inverse_phi_len, rate_len_have);
	ind_len_mis  ~ gamma(inverse_phi_len, rate_len_mis);

}

generated quantities {
// ------------------------------ generated quantities ------------------------------
 
          
}


