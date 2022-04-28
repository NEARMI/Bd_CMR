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
	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years and all populations)	
	int<lower=1> ind_per_period_bd;			    // n_ind * year (the unique periods of time in which each individual's bd is estimated)
	int<lower=1> ind_occ;			   	    // n_ind * all sampling periods (all events in which each individual could potentially have been captured)
	int<lower=1> ind_occ_min1;		 	    // n_ind * all sampling periods except the last 
	int<lower=1> n_days;				    // number of sampling occasions
	int<lower=1> n_sex;					    // Number of sex entries (M, F, but possibly U)
	int<lower=0> n_pop_year;			    // Number of years in which sampling occurred
	
  // dimensional and bookkeeping params (vectors)	
	int<lower=1> ind_occ_size[n_ind];		    // Number of sampling periods for all individuals
	int<lower=1> ind_occ_min1_size[n_ind];		    // Number of sampling periods -1 for all individuals
	int<lower=1> phi_first_index[n_ind];		    // The indexes of phi corresponding to the first entry for each individual
	int<lower=1> p_first_index[n_ind];	            // The indexes of p corresponding to the first entry for each individual
	matrix[n_ind, n_sex] ind_sex;		  	    // Sex of each individual
	
  // long vector indices for observation model (p)
	int<lower=0> ind_occ_rep[ind_occ];		    // Index vector of all individuals (each individual repeated the number of sampling occasions)
	int<lower=0> p_day[ind_occ];			    // individual day identifier to try and estimate detection by day
  
  // long vector indices for survival model (phi)
	int<lower=0> ind_occ_min1_rep[ind_occ_min1];	    // Index vector of all individuals (each individual repeated the number of sampling occasions -1)
	int<lower=0> phi_bd_index[ind_occ_min1];	    // which entries of latent bd correspond to each entry of phi

  // long vector indices for bd (bd)
	int<lower=0> ind_bd_rep[ind_per_period_bd];	    // Index of which individual is associated with each estimated Bd value
	int<lower=0> bd_time[ind_per_period_bd];	    // Index of which year is associated with each estimated Bd value
	  
  // covariates (bd)
	int<lower=1> N_bd;				    // Number of defined values for bd
 	real X_bd[N_bd];			   	    // The bd values 
	int<lower=0> x_bd_index[N_bd];			    // entries of X (latent bd) that have a corresponding real measure to inform likelihood with

  // covariates (length)
	int<lower=0> n_ind_len_have;			    // Number of individuals that we have length data	  
	int<lower=0> n_ind_len_mis;			    // Number of individuals with missing length data
	int<lower=0> ind_len_which_have[n_ind_len_have];    // Index of individuals that we have length data
	int<lower=0> ind_len_which_mis[n_ind_len_mis];      // Index of individuals with missing length data
	vector[n_ind_len_have] ind_len_have;		    // The actual length values that we have

	matrix[n_ind_len_have, n_sex] ind_len_sex_have;	    // The sex of all individuals that we have lengths for, in model matrix form
	matrix[n_ind_len_mis, n_sex] ind_len_sex_mis;	    // The sex of all individuals that we don't have lengths for, in model matrix form

  // captures
	int<lower=1> N_y;				    // Number of defined values for captures
  	int<lower=0, upper=1> y[N_y];		            // The capture values 
	vector<lower=0>[n_days] n_capt_per_day;

  // indices of phi, p, and chi that are 0, 1, or estimated, and which entries inform the likelihood.
  // set up in R to avoid looping over the full length of phi and p here. See R code for details
	int<lower=1> n_phi_zero;  
	int<lower=1> n_phi_one;           
	int<lower=1> n_phi_off;    
	int<lower=1> phi_zero_index[n_phi_zero];  
	int<lower=1> phi_one_index[n_phi_one];
	int<lower=1> phi_off_index[n_phi_off]; 

	int<lower=1> n_p_zero;
	int<lower=1> n_p_est;
	int<lower=1> p_zero_index[n_p_zero];
	int<lower=1> p_est_index[n_p_est];

	int<lower=1> n_phi_ll;
	int<lower=1> n_p_ll;
	int<lower=1> which_phi_ll[n_phi_ll];
	int<lower=1> which_p_ll[n_p_ll];
	int<lower=1> which_chi_ll[n_ind];

}

parameters {
// ------------------------------ parameters ------------------------------

// -----
// bd submodel
// ----- 

	vector[n_pop_year] beta_bd_year;		 // Each year gets a unique Bd intercept
	real beta_bd_len;				 // individual-specific length effect on bd levels
	
	real<lower=0> bd_delta_sigma;			 // change in Bd by individual (normal random effect variance)		 
	real bd_delta_eps[n_ind];                        // the conditions modes of the random effect (each individual's intercept (for now))

	real<lower=0> bd_obs;    			 // observation noise for observed Bd compared to underlying state	

// -----
// survival
// -----

	vector[2] beta_offseason;  			 // survival as a function of bd stress
	vector[n_sex] beta_offseason_sex;		 // sex effect on survival

// -----
// detection
// -----
	
	vector[n_sex] beta_p_sex;
	real<lower=0> p_day_delta_sigma;
	real p_day_delta_eps[n_days];
	
// -----
// imputed covariates: length
// -----

	real<lower=0> inverse_phi_len;		         // variance parameter for gamma regression
	vector[n_sex] beta_len_sex;			 // regression coefficient len as a function of sex
	vector[n_ind_len_mis] ind_len_mis;		 // the imputed values of len

}

transformed parameters {
// ------------------------------ transformed parameters ------------------------------

	// Individual Lengths

  	vector[n_ind_len_have] mu_len_have; 		 // the expected values for the gamma regression
  	vector[n_ind_len_have] rate_len_have; 	 	 // rate parameter for the gamma distribution

  	vector[n_ind_len_mis] mu_len_mis; 		 // the expected values (linear predictor) for the missing len values
  	vector[n_ind_len_mis] rate_len_mis; 		 // rate parameter for the gamma distribution for the missing len values

	vector[n_ind] ind_len;				 // all individual len (combining data and imputed values)
	vector[n_ind] ind_len_scaled;			 // all individual len scaled

	// bd

	real bd_ind[n_ind];				 // individual random effect deviates
	vector[ind_per_period_bd] X;		         // each individual's estimated bd per year

	// Survival and detection processes

	vector<lower=0,upper=1>[ind_occ_min1] phi;       // survival from t to t+1, each individual repeated the number of times its population was measured
	vector<lower=0,upper=1>[ind_occ] p;              // detection at time t
	real<lower=0,upper=1> chi[ind_occ];              // probability an individual will never be seen again

	vector[n_days] p_day_dev;

	vector<lower=0,upper=1>[n_days] p_per_day;	 // average detection per day


// -----
// Imputed NA Data values
// -----

	// Individual Length

  	mu_len_have   = exp(ind_len_sex_have * beta_len_sex);  				// linear predictor for len regression 
  	rate_len_have = rep_vector(inverse_phi_len, n_ind_len_have) ./ mu_len_have;	// gamma parameter from mean

  	mu_len_mis    = exp(ind_len_sex_mis * beta_len_sex);   				// predict for missing using estimated coefficients 	
  	rate_len_mis = rep_vector(inverse_phi_len, n_ind_len_mis) ./ mu_len_mis;	// gamma parameter from mean
		
	ind_len[ind_len_which_have] = ind_len_have;				 	// filling in the complete vector of ind_mehg with the data
	ind_len[ind_len_which_mis]  = ind_len_mis;       				// filling in the complete vector of ind_mehg with the imputed values

	ind_len_scaled = (ind_len - mean(ind_len))/sd(ind_len);

// -----
// bd submodel, contained to estimating within-season bd
// -----

  // linear predictor for intercept for bd-response. Overall intercept + pop-specific intercept + individual random effect deviate
	for (i in 1:n_ind) {
  	  bd_ind[i]  = bd_delta_sigma * bd_delta_eps[i];  
	}

  // latent bd model before obs error
	for (t in 1:ind_per_period_bd) {
	  X[t] = beta_bd_year[bd_time[t]] + bd_ind[ind_bd_rep[t]] + beta_bd_len * ind_len_scaled[ind_bd_rep[t]];      
        }


// -----
// Survival probability over the whole period
// -----

	phi[phi_zero_index] = rep_vector(0, n_phi_zero);
	phi[phi_one_index]  = rep_vector(1, n_phi_one);

	phi[phi_off_index]  = inv_logit(
ind_sex[ind_occ_min1_rep[phi_off_index], ] * beta_offseason_sex + 
beta_offseason[1] * X[phi_bd_index[phi_off_index]] +
beta_offseason[2] * ind_len_scaled[ind_occ_min1_rep[phi_off_index]]
);


// -----
// Detection probability over the whole period
// -----

	for (t in 1:n_days) {
  	  p_day_dev[t]  = p_day_delta_sigma * p_day_delta_eps[t];
	  p_per_day[t] = inv_logit(p_day_dev[t]);  
	}

	p[p_zero_index] = rep_vector(0, n_p_zero);
	p[p_est_index]  = inv_logit(ind_sex[ind_occ_rep[p_est_index], ] * beta_p_sex + p_day_dev[p_day[p_est_index]]);
	
// -----
// Probability of never detecting an individual again after time t
// -----

  // For each individual calculate the probability it won't be captured again
	for (i in 1:n_ind) {
	 chi[p_first_index[i]:(p_first_index[i] + ind_occ_size[i] - 1)] = prob_uncaptured(ind_occ_size[i], 
              segment(p, p_first_index[i], ind_occ_size[i]), segment(phi, phi_first_index[i], ind_occ_min1_size[i]));
	}

}

model {
// ------------------------------ model ------------------------------

// -----
// Priors
// -----


// Bd Model Priors

	bd_delta_sigma  ~ inv_gamma(8, 15);
	bd_obs          ~ inv_gamma(10, 4);
	bd_delta_eps    ~ normal(0, 3);
	beta_bd_len     ~ normal(0, 3);

// Survival Priors

	beta_bd_year        ~ normal(0, 3);
	beta_offseason[1]   ~ normal(0, 1.45);
	beta_offseason[2]   ~ normal(0, 1.45);
	beta_offseason_sex  ~ normal(0, 1.45);

// Detection Priors

	beta_p_sex        ~ normal(0, 1.45);
	p_day_delta_sigma ~ inv_gamma(8, 15);
	p_day_delta_eps   ~ normal(0, 1.45);


// Imputed Covariates Priors: len

	inverse_phi_len  ~ inv_gamma(8, 15);	
	beta_len_sex     ~ normal(0, 3);


// -----
// Imputed NA Data values
// -----

	ind_len_have ~ gamma(inverse_phi_len, rate_len_have);
	ind_len_mis  ~ gamma(inverse_phi_len, rate_len_mis);

// -----
// Bd Process and Data Model
// -----

// observed bd is the linear predictor + some observation noise

	X_bd ~ normal(X[x_bd_index], bd_obs);
    
// -----
// Capture model
// -----

	1 ~ bernoulli(phi[which_phi_ll]);
	y[which_p_ll] ~ bernoulli(p[which_p_ll]);
	1 ~ bernoulli(chi[which_chi_ll]);

}


generated quantities {
// ------------------------------ generated quantities ------------------------------
 
  vector<lower=0>[n_days] pop_size;
  pop_size = n_capt_per_day ./ p_per_day;

}
