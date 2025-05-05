functions {

	// Function to build the required chi_t probability estimator for the mark recapture model, which
	 // is the probability of never recapturing an individual again after capturing them at time t
	
	matrix prob_uncaptured(int n_ind, int n_occasions, matrix p, matrix phi) {

	matrix[n_ind, n_occasions] chi;    // chi for each capture date and individual 

	for (i in 1:n_ind) {
         chi[i, n_occasions] = 1.0;        // on the last sampling date the probability is one

	for (t in 1:(n_occasions - 1)) {   // loop over from the first to the second to last sampling date and multiply out the probabilities
         int t_curr = n_occasions - t;
         int t_next = t_curr + 1;
	 chi[i, t_curr] = (1 - phi[i, t_curr]) + phi[i, t_curr] * (1 - p[i, t_next - 1]) * chi[i, t_next];		
      }
    }

  return chi;
 }

}

data {

	// NOTE: to convert matrix to database-style structure anything with dimensions equal to the number of individuals or sampling occasions need to be melted
	//       as long as times can not vary by population it will be ok to just stack populations and have a vector that determines which individuals belong to which pop

  // dimensional and bookkeeping params (single vals)
	int<lower=1> n_periods;				    // Total number of seasons/years (the "on" period where sampling occurs) over which individuals are captured
	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years)
	int<lower=n_periods> n_times;		    	    // Sum of total time points modeled within each season across all seasons/years  
	int<lower=1> times_within;			    // number of time periods in each season		
	int<lower=n_periods> n_occasions;		    // Total number of sampling days (a subset of times) across all seasons/years 
	int<lower=1> n_occ_min1;			    // ^^ just throwing out the very last sampling day
	int<lower=0> ind_time;

  // dimensional and bookkeeping params (vectors)
	int<lower=0> time[n_times];		 	                          // Vector indicating time (e.g., weeks) *!within each season!*
	int<lower=0, upper=n_times> time_per_period[times_within, n_periods];     // Matrix of indices of time per period for subsetting X
	int<lower=0> sampling_events[n_occasions]; 	    // Indices from 1:n_times on which a sampling event occurred
	int<lower=0,upper=1> offseason[n_occ_min1];	    // Vector indicating the last sampling periods of each season which gains offseason characteristics
	int<lower=0> periods[n_times];			    // Vector designating periods for bd model (all times)
	int<lower=0> periods_occ[n_occasions];		    // Vector designating periods for observational model (all occasions)
	
	int<lower=0> ind_time_rep[ind_time];
	int<lower=0> time_rep[ind_time];
	int<lower=0> temp_rep[ind_time];
	  
  // covariates
	int<lower=1> N_bd;				    // Number of defined values for bd
  	array[N_bd] int<lower=1, upper=n_ind> ii_bd;	    // individual index that defines each bd entry
  	array[N_bd] int<lower=1, upper=n_occasions> tt_bd;  // occasion index that defines each bd entry
  	array[N_bd] int<lower=0> X_bd;			    // The bd values 

	real temp[n_times];				    // Temperature covariate
	int<lower=0> time_gaps[n_occ_min1];  	 	    // Elapsed time between each sampling event 

	int<lower=0> samp_indices[total_swabs];

  // captures
	int<lower=1> N_y;				    // Number of defined values for captures
  	array[N_y] int<lower=1, upper=n_ind> ii_y;	    // individual index
  	array[N_y] int<lower=1, upper=n_occasions> tt_y;    // occasion index 
  	array[N_y] int<lower=0, upper=1> y;		    // The capture values 

  	int<lower=0> first[n_ind];         		    // Capture event in which each individual was first captured
  	int<lower=0> last[n_ind];         		    // Capture event in which each individual was last captured
	int<lower=0> present[n_ind, n_periods];		    // Index to show in which period each individual was captured at least once

}

transformed data {

	int<lower=0> ind_occ_min1;
	int<lower=0> ind_occ;

	ind_occ_min1 = n_ind * n_occ_min1;
	ind_occ      = n_ind * n_occasions;

} 

parameters {

	vector[2] beta_phi;                  		 // intercept and slope coefficient for survival
        vector[2] beta_p;				 // intercept and slope coefficient for detection
	real<upper=0> beta_timegaps;			 // coefficient to control for the variable time between sampling events
	real<upper=0> beta_offseason;			 // season survival probability, maybe maybe not as a function of bd

	vector[3] beta_bd;				 // intercept and two slope coefficients for grand mean change in bd over time

	real<lower=0> bd_delta_sigma;			 // change in Bd by individual (normal random effect variance)
	real bd_delta_eps[n_ind];			 // the conditions modes of the random effect (each individual's intercept (for now))

	real<lower=0> bd_obs;    			 // observation noise for observed Bd compared to underlying state

	matrix<lower=0,upper=1>[n_ind, n_periods] gamma; // probability of each individual seen in subsequent periods actually having been in the population previously

}

transformed parameters {

	real<lower=0,upper=1> phi[ind_occ_min1];   // survival from t to t+1
	real<lower=0,upper=1> p[ind_occ];          // detection at time t
	real<lower=0,upper=1> chi[ind_occ];        // probability an individual will never be seen again
 
	matrix[n_ind, n_times] X;	   	   // Estimated "true" bd for all of the caught individuals with no bd measured
	matrix[n_ind, n_periods] X_max; 	   // summaries of X	
		
	real bd_ind[n_ind];                        // Individual random effect deviates


	//// bd submodel, contained to estimating within-season bd

	for (i in 1:n_ind) {
  	  bd_ind[i] = beta_bd[1] + bd_delta_sigma  * bd_delta_eps[i];
	}

	 for (t in 1:n_times) {
	  X[i, t]   = bd_ind[i] + beta_bd[2] * time[t] + beta_bd[3] * temp[t];
	 }

	 for (tp in 1:n_periods) {
          X_max[i, tp] = max(X[i, time_per_period[1, tp]:time_per_period[times_within, tp]]);
	 }

        }


	//// Survival and detection probability over the whole period

	// I think it would be best to form an indexing vector using first and last, stick it together, and refer to it with an indexing vector of individual
         // at which point the loop can be over the whole length of phi?

	for (i in 1:n_ind) {

	  phi[((i - 1) * n_occ_min1 + 1):((i - 1) * n_occ_min1 + first[i] - 1)] = 0; 
		
         for (t in first[i]:n_occ_min1) { // for all events after the first time an individual was caught, estimate its mortality probability
          
	  phi[((i - 1) * n_occ_min1 + first[i] + t)] = inv_logit(
beta_phi[1]                   + 
beta_timegaps  * time_gaps[t] +
beta_offseason * offseason[t] +	
beta_phi[2]    * X[i, sampling_events[t]]
);

	}

}

        // for all events estimate an individuals detection probability -- conditional on whether the individual was known to be present at that time or not
	
	for (i in 1:n_ind) {
	 for (t in 1:n_occasions) {

	   if (present[i, periods_occ[t]] == 1) {
	  p[(i - 1) * n_occasions + t] = inv_logit(beta_p[1] + beta_p[2] * X[i, sampling_events[t]]);
	   } else {
	  p[(i - 1) * n_occasions + t] = inv_logit(beta_p[1] + beta_p[2] * X[i, sampling_events[t]]) * gamma[i, periods_occ[t]];
	   }

	 }
	}

	//// Probability of never detecting an individual again after time t

	for (i in 1:ind) {
	  chi = prob_uncaptured(n_ind, n_occasions, p, phi);
	}
}

model {

	// Priors

	beta_phi[1] ~ normal(0, 5);
	beta_phi[2] ~ normal(0, 5);
	beta_p[1] ~ normal(0, 5);
	beta_p[2] ~ normal(0, 5);
	beta_bd[1] ~ normal(0, 5);
	beta_bd[2] ~ normal(0, 5);
	beta_bd[3] ~ normal(0, 5);
	beta_timegaps ~ normal(0, 5);
	beta_offseason ~ normal(0, 5);

	bd_delta_sigma ~ inv_gamma(1, 1);
	bd_delta_eps ~ normal(0, 2);
	bd_obs ~ inv_gamma(1, 1);

	for (i in 1:n_ind) {
         for (j in 1:n_periods) {
          gamma[i, j] ~ uniform(0, 1);
	 }
	}

	//// Bd Process and Data Model

	for (t in 1:N_bd) {
          X_bd[t] ~ normal(X[ii_bd[t], tt_bd[t]], bd_obs); 
	} 
    
  
	//// Capture model

	 for (i in 1:n_ind) {
	  for (t in (first[i] + 1):last[i]) {
	   // survived (the one) ~ bernoulli[phi] -- because we saw an individual again in the future we know it survived this period
	   1 ~ bernoulli(phi[((i - 1) * n_occ_min1 + t)]);     
	  }
	}

	
	  for (t in 1:N_y) {
   	   y[N_y] ~ bernoulli(p[N_y]);      // detection can inform y up until the last time the individual was caught, and then it is inseparable from chi
	  }

	   1 ~ bernoulli(chi[i, last[i]]);  // the probability of an animal never being seen again after the last time it was captured

	  }

}

generated quantities {
 

          
}


