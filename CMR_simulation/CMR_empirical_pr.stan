functions {

	// Function to build the required chi_t probability estimator for the mark recapture model, which
	 // is the probability of never recapturing an individual again after capturing them at time t
	// rewritten to run for an individual at a time given variable numbers of capture opportunities by individual (e.g. in different populations)
	
	real[] prob_uncaptured(int n_occ, real[] p_sub, real[] phi_sub) {

	real chi_sub[n_occ];          // chi for each capture date and individual 

        chi_sub[n_occ] = 1.0;         // on the last sampling date the probability is one

	for (t in 1:(n_occ - 1)) {    // loop over from the first to the second to last sampling date and multiply out the probabilities
         int t_curr = n_occ - t;
         int t_next = t_curr + 1;
	 chi_sub[t_curr] = (1 - phi_sub[t_curr]) + phi_sub[t_curr] * (1 - p_sub[t_next - 1]) * chi_sub[t_next];		
      }

  return chi_sub;

 }

}

data {

// -----
// Notes
// -----
// Given long comments, this model is best read in full screen on an external monitor


  // dimensional and bookkeeping params (single vals)
	int<lower=1> n_pop;				    // Number of distinct populations (sampling areas)
	int<lower=1> n_periods;				    // Total number of seasons/years (the "on" period where sampling occurs) over which individuals are captured
	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years)
	int<lower=n_periods> n_times;		    	    // Sum of total time points modeled within each season across all seasons/years 
 
	int<lower=1> times_within;			    // number of time periods in each season		

	int<lower=0> ind_occ;			   	    // n_ind * n_occasions, summed over the sampling of all populations
	int<lower=0> ind_occ_min1;		 	    // n_ind * n_occ_min1, summed over the sampling of all populations
	
  // dimensional and bookkeeping params (vectors)
	int<lower=0> time[n_times];		 	                          // Vector indicating time (e.g., weeks) *!within each season!*
	int<lower=0, upper=n_times> time_per_period[times_within, n_periods];     // Matrix of indices of time per period for subsetting X			 
	int<lower=0> periods[n_times];			                          // Vector designating periods for bd model (all times)
	
	int<lower=1> ind_occ_size[n_ind];					  // Number of sampling periods for all individuals
	int<lower=1> ind_occ_min1_size[n_ind];					  // Number of sampling periods -1 for all individuals
	
	int<lower=1> ind_in_pop[n_ind];						  // population in which each individual resides	

	int<lower=1> phi_first_index[n_ind];				          // The indexes of phi corresponding to the first entry for each individual
	int<lower=1> p_first_index[n_ind];				          // The indexes of p corresponding to the first entry for each individual
	
  // long vector indices for observation model (p)
	int<lower=0> ind_occ_rep[ind_occ];		    // Index vector of all individuals (each individual repeated the number of sampling occasions)
	int<lower=0> sampling_events_p[ind_occ];  	    // The date on which each sampling event occurred for each individual
	int<lower=0> periods_occ[ind_occ];		    // Vector designating periods for observational model (all occasions)
	int<lower=0> p_zeros[ind_occ];			    // Observation times for each individual in which we do not know if that individual is present
	int<lower=0> pop_p[ind_occ];			    // population index for detection predictors
  
  // long vector indices for survival model (phi)
	int<lower=0> ind_occ_min1_rep[ind_occ_min1];	    // Index vector of all individuals (each individual repeated the number of sampling occasions -1)
	int<lower=0> sampling_events_phi[ind_occ_min1];     // The date on which each sampling event occurred (minus the last one) for each individual
	int<lower=0, upper=1> offseason[ind_occ_min1];	    // Vector indicating the last sampling periods of each season which gains offseason characteristics
	int<lower=0> phi_zeros[ind_occ_min1];		    // Observation times for each individual in advance of first detecting that individual
	int<lower=0> pop_phi[ind_occ_min1];		    // population index for mortality predictors
	  
  // covariates
	int<lower=1> N_bd;				    // Number of defined values for bd
 	real X_bd[N_bd];			   	    // The bd values 
  	int<lower=1, upper=n_ind> ii_bd[N_bd];	            // individual index that defines each bd entry
  	int<lower=1, upper=n_times> tt_bd[N_bd];            // occasion index that defines each bd entry
	matrix[n_times, n_pop] temp;			    // Temperature covariate in each population
	int<lower=0> time_gaps[ind_occ_min1];  	 	    // Elapsed time between each sampling event 

  // captures
	int<lower=1> N_y;				    // Number of defined values for captures
  	int<lower=0, upper=1> y[N_y];		            // The capture values 
  	int<lower=0> first[n_ind];         		    // Capture event in which each individual was first captured
  	int<lower=0> last[n_ind];         		    // Capture event in which each individual was last captured

}

transformed data {

} 

parameters {

// -----
// bd submodel
// -----

	real beta_bd_int_pop;				 // population-specific intercepts in bd load
	vector[3] beta_bd;				 // two slope coefficients for grand mean change in bd over time + slope for temp

	real<lower=0> bd_delta_pop_sigma;		 // change in Bd by pop (normal random effect variance)
	real bd_delta_pop_eps[n_pop];			 // the conditions modes of the random effect (each populations intercept (for now))

	real<lower=0> bd_delta_sigma;			 // change in Bd by individual (normal random effect variance)		 
	real bd_delta_eps[n_ind];                        // the conditions modes of the random effect (each individual's intercept (for now))

	real<lower=0> bd_obs;    			 // observation noise for observed Bd compared to underlying state	


// -----
// survival
// -----

	real beta_phi;       	          	 	 // grand intercept and slope for survival
	real<upper=0> beta_phi_slope_pop;	    	 // population specific slopes for survival (as bd increases survival decreases)

	real<lower=0> phi_delta_pop_sigma;		 // change in Bd by individual (normal random effect variance)
	real phi_delta_pop_eps[n_pop];
        
	real<upper=0> beta_timegaps;			 // coefficient to control for the variable time between sampling events
	real<upper=0> beta_offseason;			 // season survival probability, maybe maybe not as a function of bd


// -----
// detection
// -----

	vector[2] beta_p;				 // intercept and slope coefficient for detection
	

// -----
// other
// -----	

	matrix<lower=0,upper=1>[n_ind, n_periods] gamma; // probability of each individual seen in subsequent periods actually having been in the population previously

}

transformed parameters {

	real<lower=0,upper=1> phi[ind_occ_min1];   // survival from t to t+1, each individual repeated the number of times its population was measured
	real<lower=0,upper=1> p[ind_occ];          // detection at time t
	real<lower=0,upper=1> chi[ind_occ];        // probability an individual will never be seen again
 
	matrix[n_ind, n_times] X;	   	   // Estimated "true" bd for all of the caught individuals with no bd measured
	matrix[n_ind, n_periods] X_max; 	   // summaries of X	
		
	real bd_ind[n_ind];                        // Individual random effect deviates
	real bd_pop[n_pop];
	real phi_pop[n_pop];

// -----
// bd submodel, contained to estimating within-season bd
// -----

	for (pp in 1:n_pop) {
	 bd_pop[pp]  = beta_bd_int_pop + bd_delta_pop_sigma * bd_delta_pop_eps[pp];
	 phi_pop[pp] = beta_phi_slope_pop + phi_delta_pop_sigma * phi_delta_pop_eps[pp];
	} 

	for (i in 1:n_ind) {
	    
		// linear predictor for intercept for bd-response. Overall intercept + pop-specific intercept + individual random effect deviate

  	  bd_ind[i] = bd_delta_sigma * bd_delta_eps[i];  

		// latent bd model before obs error

	 for (t in 1:n_times) {
	  X[i, t]   = (bd_pop[ind_in_pop[i]] + bd_ind[i])    +
		      beta_bd[1] * time[t]                   +
		      beta_bd[2] * square(time[t])           + 
		      beta_bd[3] * temp[t, ind_in_pop[i]];      
	 }

	 for (tp in 1:n_periods) {

		// calculation of the maximum bd load experienced in a year (could also be cumulative)

          X_max[i, tp] = max(X[i, time_per_period[1, tp]:time_per_period[times_within, tp]]);	
   
	 }

        }

// -----
// Survival probability over the whole period
// -----

	for (t in 1:ind_occ_min1) {

	 if (phi_zeros[t] == 1) {			// phi_zeros is 1 before an individual is caught for the first time
           phi[t] = 0;					// must be non-na values in stan, but the likelihood is only informed from first capture onward
	 } else {
           phi[t] = inv_logit(
                      beta_phi                       + 
                      beta_timegaps  * time_gaps[t]  +
                      beta_offseason * offseason[t]  +	
                      phi_pop[pop_phi[t]] * X[ind_occ_min1_rep[t], sampling_events_phi[t]]
                    );
	 }  

	}

// -----
// Detection probability over the whole period
// -----
	
	for (t in 1:ind_occ) {

		// p_zeros is = 1 in each season prior to an individual being caught for the first time
		// p gets scaled in these years in an attempt to scale the probability as a function of bd given that we don't know if the individual was there

	 if (p_zeros[t] == 1) {				
          p[t] = inv_logit(beta_p[1] + beta_p[2] * X[ind_occ_rep[t], sampling_events_p[t]]);
	 } else {
          p[t] = inv_logit(beta_p[1] + beta_p[2] * X[ind_occ_rep[t], sampling_events_p[t]]) * gamma[ind_occ_rep[t], periods_occ[t]];
	 }

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

// -----
// Priors
// -----

	beta_bd[1] ~ normal(0, 5);
	beta_bd[2] ~ normal(0, 5);
	beta_bd[3] ~ normal(0, 5);
	beta_phi   ~ normal(0, 5);
	beta_p[1]  ~ normal(0, 5);
	beta_p[2]  ~ normal(0, 5);

	beta_timegaps  ~ normal(0, 5);
	beta_offseason ~ normal(0, 5);

	beta_bd_int_pop    ~ normal(0, 5);
	beta_phi_slope_pop ~ normal(0, 5);

	bd_delta_sigma      ~ inv_gamma(1, 1);
	bd_delta_pop_sigma  ~ inv_gamma(1, 1);
	phi_delta_pop_sigma ~ inv_gamma(1, 1);
	
	bd_obs              ~ inv_gamma(1, 1);

	for (pp in 1:n_pop) {
	 bd_delta_pop_eps[pp]  ~ normal(0, 5);
	 phi_delta_pop_eps[pp] ~ normal(0, 5);
	}

	for (i in 1:n_ind) {
	  bd_delta_eps[i]   ~ normal(0, 3);
         for (j in 1:n_periods) {
          gamma[i, j]       ~ uniform(0, 1);
	 }
	}

// -----
// Bd Process and Data Model
// -----

		// observed bd is the linear predictor + some observation noise

	for (t in 1:N_bd) {
          X_bd[t] ~ normal(X[ii_bd[t], tt_bd[t]], bd_obs); 
	} 
    
// -----
// Capture model
// -----

	 for (i in 1:n_ind) {
	  for (t in (first[i] + 1):last[i]) {			
	   1 ~ bernoulli(phi[phi_first_index[i] - 1 + t - 1]);    		   // Survival _to_ t (from phi[t - 1]) is 1 because we know the individual lived in that period 
	  }
	  for (t in 1:last[i]) {
	   y[p_first_index[i] - 1 + t] ~ bernoulli(p[p_first_index[i] - 1 + t]);   // Capture given detection
	  }

	   1 ~ bernoulli(chi[p_first_index[i] - 1 + last[i]]);  		   // the probability of an animal never being seen again after the last time it was captured

	  }

}

generated quantities {
 
          
}


