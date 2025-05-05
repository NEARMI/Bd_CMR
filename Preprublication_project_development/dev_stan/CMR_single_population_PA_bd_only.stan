data {
// ------------------------------ data ------------------------------

	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years and all populations)
	int<lower=0> ind_occ;			   	    // n_ind * all sampling periods (all events in which each individual could potentially have been captured)
	int<lower=1> ind_per_period_bd;			    // n_ind * year (the unique periods of time in which each individual's bd is estimated)
	int<lower=0> ind_occ_rep[ind_occ];		    // Index vector of all individuals (each individual repeated the number of sampling occasions)
	vector[ind_occ] yday;			   	    // Scaled Julian day (for now a placeholder, will want to get to a meaningful measure of temp soon)
	vector[ind_occ] yday_sq;			    // Square of scaled Julian Day
	int<lower=0> X_first_index[ind_per_period_bd];	    // First index for each section of Bd values (occasions in a year)
	int<lower=0> X_gap[ind_per_period_bd];		    // Number of entries in X for each section of Bd values (occasions in a year)
	int<lower=1> N_bd;				    // Number of defined values for bd
	real X_bd[N_bd];			   	    // The bd values 	
	int<lower=0> x_bd_index_full[N_bd];		    // entries of X (latent bd) that have a corresponding real measure to inform likelihood with

}

parameters {
// ------------------------------ parameters ------------------------------
	
	real beta_bd;
	real beta_bd_day;				 // linear term for Bd over time
	real beta_bd_day_sq;				 // quadratic term for Bd over time
	real<lower=0> bd_delta_sigma;			 // change in Bd by individual (normal random effect variance)		 
	real bd_delta_eps[n_ind];                        // the conditions modes of the random effect (each individual's intercept (for now))
	real<lower=0> bd_obs;    			 // observation noise for observed Bd compared to underlying state	

}

transformed parameters {
// ------------------------------ transformed parameters ------------------------------

	// bd

	real bd_ind[n_ind];				 // individual random effect deviates
	real X[ind_occ];		    	         // each individual's estimated bd per year
	real X_max[ind_per_period_bd];			 // estimated max bd experienced by an individual in a given year

  // linear predictor for intercept for bd-response. Overall intercept + pop-specific intercept + individual random effect deviate
	for (i in 1:n_ind) {
  	  bd_ind[i]  = bd_delta_sigma * bd_delta_eps[i] + beta_bd;  
	}

  // latent bd model before obs error
	for (t in 1:ind_occ) {
	  X[t] = bd_ind[ind_occ_rep[t]] + beta_bd_day * yday[t] + beta_bd_day_sq * yday_sq[t];      
        }

  // maximum estimated bd load experienced by an individual between offseasons
	for (t in 1:ind_per_period_bd) {
	  X_max[t] = max(segment(X, X_first_index[t], X_gap[t]));
	} 
}

model {
// ------------------------------ model ------------------------------

// Bd Model Priors

	bd_delta_sigma    ~ inv_gamma(8, 15);
	bd_obs            ~ inv_gamma(10, 4);

	beta_bd		  ~ normal(0, 5);
	beta_bd_day	  ~ normal(0, 3);
	beta_bd_day_sq    ~ normal(0, 3);

	for (i in 1:n_ind) {
	  bd_delta_eps[i] ~ normal(0, 3);
	}

// observed bd is the linear predictor + some observation noise

	for (t in 1:N_bd) {
          X_bd[t] ~ normal(X[x_bd_index_full[t]], bd_obs); 
	} 

}
