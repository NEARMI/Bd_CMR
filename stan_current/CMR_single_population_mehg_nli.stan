functions {
// ------------------------------ functions ------------------------------

// -----
// chi_t probability estimator
// -----

// Function to build the required chi_t probability estimator for the mark recapture model, which
 // is the probability of never recapturing an individual again after capturing them at time t
  // rewritten to run for an individual at a time given variable numbers of capture opportunities by individual (e.g. in different populations)
	
	real[] prob_uncaptured(int n_occ, real[] p_sub, real[] phi_sub) {

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
	int<lower=1> ind_per_period_bd;			    // n_ind * year (the unique periods of time in which each individual's bd is estimated)
	int<lower=0> ind_occ;			   	    // n_ind * all sampling periods (all events in which each individual could potentially have been captured)
	int<lower=0> ind_occ_min1;		 	    // n_ind * all sampling periods except the last 
	int<lower=0> n_days;				    // number of sampling occasions
	int n_sex;					    // Number of sex entries (M, F, but possibly U)
	int<lower=0> n_pop_year;			 // Number of years in which sampling occurred
	
  // dimensional and bookkeeping params (vectors)	
	int<lower=1> ind_occ_size[n_ind];		    // Number of sampling periods for all individuals
	int<lower=1> ind_occ_min1_size[n_ind];		    // Number of sampling periods -1 for all individuals
	int<lower=1> phi_first_index[n_ind];		    // The indexes of phi corresponding to the first entry for each individual
	int<lower=1> p_first_index[n_ind];	            // The indexes of p corresponding to the first entry for each individual
	
  // long vector indices for observation model (p)
	int<lower=0> ind_occ_rep[ind_occ];		    // Index vector of all individuals (each individual repeated the number of sampling occasions)
	int<lower=0> p_month[ind_occ];		            // Vector designating shorter periods (here month) for observational model (all occasions)
	int<lower=0> p_zeros[ind_occ];			    // Observation times for each individual in which we do not know if that individual is present
	int<lower=0> p_bd_index[ind_occ];		    // which entries of latent bd correspond to each entry of p

	int<lower=0> p_day[ind_occ];			    // individual day identifier to try and estimate detection by day
  
  // long vector indices for survival model (phi)
	int<lower=0> ind_occ_min1_rep[ind_occ_min1];	    // Index vector of all individuals (each individual repeated the number of sampling occasions -1)
	int<lower=0, upper=1> offseason[ind_occ_min1];	    // Vector indicating the last sampling periods of each season which gains offseason characteristics
	int<lower=0> phi_year[ind_occ_min1];		    // Same as p_year but without last sampling event
	int<lower=0> phi_zeros[ind_occ_min1];		    // Observation times for each individual in advance of first detecting that individual
	int<lower=0> phi_ones[ind_occ_min1];	            // Time periods where we force survival to be 1 (assuming a closed population)
	int<lower=0> phi_bd_index[ind_occ_min1];	    // which entries of latent bd correspond to each entry of phi

  // long vector indices for bd (bd)
	int<lower=0> ind_bd_rep[ind_per_period_bd];	    // Index of which individual is associated with each estimated Bd value
	int<lower=0> bd_time[ind_per_period_bd];	    // Index of which year is associated with each estimated Bd value
	  
  // covariates (bd)
	int<lower=1> N_bd;				    // Number of defined values for bd
 	real X_bd[N_bd];			   	    // The bd values 
	int<lower=1> X_ind[N_bd];			    // Individual associated with each bd measure
	int<lower=0> x_bd_index[N_bd];			    // entries of X (latent bd) that have a corresponding real measure to inform likelihood with

	int<lower=0> bd_first_index[ind_per_period_bd];	    // First entry of latent bd associated with each individual 'by' period
	int<lower=0> bd_last_index[ind_per_period_bd];	    // Last entry of latent bd associated with each individual 'by' period

  // covariates (length)
	vector[n_ind] ind_len_have;			    // Individual lengths already scaled (named with "have" for convenience for correspondence with other model)

  // covariates (MeHg)
	int<lower=0> n_ind_mehg_have;			    // Number of individuals that we have mehg data	  
	int<lower=0> n_ind_mehg_mis;			    // Number of individuals with missing mehg data
	int<lower=0> ind_mehg_which_have[n_ind_mehg_have];  // Index of individuals that we have mehg data
	int<lower=0> ind_mehg_which_mis[n_ind_mehg_mis];    // Index of individuals with missing mehg data
	vector[n_ind_mehg_have] ind_mehg_have;		    // The actual mehg values that we have
	
  // covariates (sex)
	matrix[n_ind, n_sex] ind_sex;		  	    // Sex of each individual

  // captures
	int<lower=1> N_y;				    // Number of defined values for captures
  	int<lower=0, upper=1> y[N_y];		            // The capture values 

  	int<lower=0> first[n_ind];         		    // Capture event in which each individual was first captured
  	int<lower=0> last[n_ind];         		    // Capture event in which each individual was last captured

	vector<lower=0>[n_days] n_capt_per_day;

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

	real beta_phi;                  		 // single background intercept for survival in the offseason
	vector[4] beta_offseason;  			 // survival as a function of bd stress
	vector[n_sex] beta_offseason_sex;		 // sex effect on survival

// -----
// detection
// -----
	
	vector[n_sex] beta_p_sex;
	
	real<lower=0> p_day_delta_sigma;
	real p_day_delta_eps[n_days];
	
// -----
// imputed covariates: MeHg
// -----

	real<lower=0> inverse_phi_mehg;		         // variance parameter for gamma regression
	vector[2] beta_mehg_have;			 // regression coefficient for mehg ~ length
	vector[n_ind_mehg_mis] ind_mehg_mis;		 // the imputed values of mehg
	

}

transformed parameters {
// ------------------------------ transformed parameters ------------------------------

	// Individual MeHg

  	vector[n_ind_mehg_have] mu_mehg_have; 		 // the expected values for the gamma regression
  	vector[n_ind_mehg_have] rate_mehg_have; 	 // rate parameter for the gamma distribution

  	vector[n_ind_mehg_mis] mu_mehg_mis; 		 // the expected values (linear predictor) for the missing mehg values
  	vector[n_ind_mehg_mis] rate_mehg_mis; 		 // rate parameter for the gamma distribution for the missing mehg values

	vector[n_ind] ind_mehg;				 // all individual mehg (combining data and imputed values)
	vector[n_ind] ind_mehg_scaled;			 // all individual mehg scaled
	real ind_mehg_mean;				 // mean of ind_mehg
	real ind_mehg_sd;				 // sd of ind_mehg	


	// bd

	real bd_ind[n_ind];				 // individual random effect deviates
	real X[ind_per_period_bd];		         // each individual's estimated bd per year


	// Survival and detection processes

	real<lower=0,upper=1> phi[ind_occ_min1];         // survival from t to t+1, each individual repeated the number of times its population was measured
	real<lower=0,upper=1> p[ind_occ];                // detection at time t
	real<lower=0,upper=1> chi[ind_occ];              // probability an individual will never be seen again

	real p_day_dev[n_days];

	vector<lower=0,upper=1>[n_days] p_per_day;	 // average detection per day

// -----
// Imputed NA Data values
// -----

	// Individual MeHg
  	
  	mu_mehg_have   = exp(beta_mehg_have[1] + beta_mehg_have[2] * ind_len_have[ind_mehg_which_have]);  // linear predictor for mehg regression	
  	rate_mehg_have = rep_vector(inverse_phi_mehg, n_ind_mehg_have) ./ mu_mehg_have;

  	mu_mehg_mis   = exp(beta_mehg_have[1] + beta_mehg_have[2] * ind_len_have[ind_mehg_which_mis]);    // linear predictor for mehg regression	
  	rate_mehg_mis = rep_vector(inverse_phi_mehg, n_ind_mehg_mis) ./ mu_mehg_mis;

	ind_mehg[ind_mehg_which_have] = ind_mehg_have;	    						    // filling in the complete vector of ind_mehg with the data and imputed values
	ind_mehg[ind_mehg_which_mis]  = ind_mehg_mis;

	ind_mehg_mean = mean(ind_mehg);			
	ind_mehg_sd   = sd(ind_mehg);

	ind_mehg_scaled = (ind_mehg - ind_mehg_mean)/ind_mehg_sd;


// -----
// bd submodel, contained to estimating within-season bd
// -----

	for (i in 1:n_ind) {
	    
		// linear predictor for intercept for bd-response. Overall intercept + pop-specific intercept + individual random effect deviate

  	  bd_ind[i]  = bd_delta_sigma * bd_delta_eps[i];  

	}

	for (t in 1:ind_per_period_bd) {

		// latent bd model before obs error

	  X[t] = beta_bd_year[bd_time[t]] + bd_ind[ind_bd_rep[t]] + beta_bd_len * ind_len_have[ind_bd_rep[t]];      

        }

// -----
// Survival probability over the whole period
// -----


	for (t in 1:ind_occ_min1) {

	 if (phi_zeros[t] == 1) {	 // phi_zeros is 1 before an individual is caught for the first time
           phi[t] = 0;			 // must be non-na values in stan, but the likelihood is only informed from first capture onward so the 0 here doesn't matter
	 } else {

	  if (offseason[t] == 0) {	 // in season survival process

	   if (phi_ones[t] == 1) { 	 // closed population assumption where survival is set to 1
	     phi[t] = 1;
           } else {
             phi[t] = inv_logit(beta_phi);
	   }

	  } else {			 // off season survival process
	     
	     phi[t] = inv_logit(
ind_sex[ind_occ_min1_rep[t], ] * beta_offseason_sex +
beta_offseason[1] * X[phi_bd_index[t]] + 
beta_offseason[2] * ind_len_have[ind_occ_min1_rep[t]] +
beta_offseason[3] * ind_mehg_scaled[ind_occ_min1_rep[t]] +
beta_offseason[4] * X[phi_bd_index[t]] * ind_mehg_scaled[ind_occ_min1_rep[t]]
);


	   }

	  }  

	 }

// -----
// Detection probability over the whole period
// -----

	for (i in 1:n_days) {
  	  p_day_dev[i]  = p_day_delta_sigma * p_day_delta_eps[i];  
	}

	for (t in 1:ind_occ) {   
	 if (p_zeros[t] == 0) {
	   p[t] = 0;
	 } else {       
           p[t] = inv_logit(ind_sex[ind_occ_rep[t], ] * beta_p_sex + p_day_dev[p_day[t]]);
	 }
	}

	for (t in 1:n_days) {
	  p_per_day[t] = inv_logit(beta_p + p_day_dev[t]);
	}
	 
	
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

	bd_delta_sigma ~ inv_gamma(8, 15);
	bd_obs         ~ inv_gamma(10, 4);

	beta_bd_len    ~ normal(0, 3);

	for (i in 1:n_ind) {
	  bd_delta_eps[i] ~ normal(0, 3);
	}

// Survival Priors

	beta_bd_year        ~ normal(0, 3);
	beta_phi            ~ normal(0, 1.95);
	beta_offseason[1]   ~ normal(0, 0.70);
	beta_offseason[2]   ~ normal(0, 0.70);
	beta_offseason[3]   ~ normal(0, 0.70);
	beta_offseason[4]   ~ normal(0, 0.70);
	beta_offseason_sex  ~ normal(0, 0.70);

// Detection Priors

	beta_p_sex        ~ normal(0, 1.45);
	p_day_delta_sigma ~ inv_gamma(8, 15);

	for (i in 1:n_days) {
	  p_day_delta_eps[i] ~ normal(0, 1.45);
	}

// Imputed Covariates Priors: MeHg

	inverse_phi_mehg  ~ inv_gamma(8, 15);	
	beta_mehg_have[1] ~ normal(0, 3);
	beta_mehg_have[2] ~ normal(0, 3);

// -----
// MeHg regression imputation
// -----

	ind_mehg_have ~ gamma(inverse_phi_mehg, rate_mehg_have);
	ind_mehg_mis  ~ gamma(inverse_phi_mehg, rate_mehg_mis);

// -----
// Bd Process and Data Model
// -----

// observed bd is the linear predictor + some observation noise

	for (t in 1:N_bd) {
          X_bd[t] ~ normal(X[x_bd_index[t]], bd_obs); 
	} 
    
// -----
// Capture model
// -----

	 for (i in 1:n_ind) {
	
	if (first[i] != last[i]) {
	  for (t in (first[i] + 1):last[i]) {			
	   1 ~ bernoulli(phi[phi_first_index[i] - 1 + t - 1]);    		   // Survival _to_ t (from phi[t - 1]) is 1 because we know the individual lived in that period 
	   y[p_first_index[i] - 1 + t] ~ bernoulli(p[p_first_index[i] - 1 + t]);   // Capture given detection
	  }
	}
	   1 ~ bernoulli(chi[p_first_index[i] - 1 + last[i]]);  		   // the probability of an animal never being seen again after the last time it was captured
	 }

}


generated quantities {
// ------------------------------ generated quantities ------------------------------
 
  vector<lower=0>[n_days] pop_size;
  pop_size = n_capt_per_day ./ p_per_day;

}
