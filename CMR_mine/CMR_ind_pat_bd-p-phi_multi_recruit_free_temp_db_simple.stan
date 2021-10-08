data {

  // dimensional and bookkeeping params (single vals)
	int<lower=1> n_periods;				    // Total number of seasons/years (the "on" period where sampling occurs) over which individuals are captured
	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years)
	int<lower=n_periods> n_times;		    	    // Sum of total time points modeled within each season across all seasons/years  
	int<lower=1> times_within;			    // number of time periods in each season		
	int<lower=n_periods> n_occasions;		    // Total number of sampling days (a subset of times) across all seasons/years 
	int<lower=1> n_occ_min1;			    // ^^ just throwing out the very last sampling day
	int<lower=0> ind_time;

	int<lower=0> ind_occ;
	int<lower=0> ind_occ_min1;
	
  // dimensional and bookkeeping params (vectors)
	int<lower=0> time[n_times];		 	                          // Vector indicating time (e.g., weeks) *!within each season!*
	int<lower=0, upper=n_times> time_per_period[times_within, n_periods];     // Matrix of indices of time per period for subsetting X
	int<lower=0> sampling_events_phi[ind_occ_min1]; 	
	int<lower=0> sampling_events_p[ind_occ];   			 
	int<lower=0,upper=1> offseason[ind_occ_min1];	    // Vector indicating the last sampling periods of each season which gains offseason characteristics
	int<lower=0> periods[n_times];			    // Vector designating periods for bd model (all times)
	int<lower=0> periods_occ[ind_occ];		    // Vector designating periods for observational model (all occasions)
	
  // long vector indices
	int<lower=0> ind_time_rep[ind_time];
	int<lower=0> time_rep[ind_time];
	real<lower=0> temp_rep[ind_time];
	  
  // covariates
	int<lower=1> N_bd;				    // Number of defined values for bd
 	real X_bd[N_bd];			   	    // The bd values 
  	int<lower=1, upper=n_ind> ii_bd[N_bd];	            // individual index that defines each bd entry
  	int<lower=1, upper=n_times> tt_bd[N_bd];            // occasion index that defines each bd entry

	real temp[n_times];				    // Temperature covariate
	int<lower=0> time_gaps[ind_occ_min1];  	 	    // Elapsed time between each sampling event 

	int<lower=0> samp_indices[N_bd];

  // captures
	int<lower=1> N_y;				    // Number of defined values for captures
  	int<lower=0, upper=1> y[N_y];		            // The capture values 

  	int<lower=0> first[n_ind];         		    // Capture event in which each individual was first captured
	int<lower=0> p_zeros[ind_occ];
	int<lower=0> phi_zeros[ind_occ_min1];


  	int<lower=0> last[n_ind];         		    // Capture event in which each individual was last captured
	int<lower=0> present[n_ind, n_periods];		    // Index to show in which period each individual was captured at least once

}

transformed data {

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

	real<lower=0,upper=1> phi[ind_occ_min1];   // survival from t to t+1, each individual repeated n_occ_min1 times 
	real<lower=0,upper=1> p[ind_occ];          // detection at time t
	real<lower=0,upper=1> chi[ind_occ];        // probability an individual will never be seen again
 
	matrix[n_ind, n_times] X;	   	   // Estimated "true" bd for all of the caught individuals with no bd measured
	matrix[n_ind, n_periods] X_max; 	   // summaries of X	
		
	real bd_ind[n_ind];                        // Individual random effect deviates


	//// bd submodel, contained to estimating within-season bd

	for (i in 1:n_ind) {
  	  bd_ind[i] = beta_bd[1] + bd_delta_sigma  * bd_delta_eps[i];

	 for (t in 1:n_times) {
	  X[i, t]   = bd_ind[i] + beta_bd[2] * time[t] + beta_bd[3] * temp[t];
	 }

	 for (tp in 1:n_periods) {
          X_max[i, tp] = max(X[i, time_per_period[1, tp]:time_per_period[times_within, tp]]);
	 }

        }


	//// Survival and detection probability over the whole period

	for (t in 1:ind_occ_min1) {

	 if (phi_zeros[t] == 1) {
           phi[t] = 0;
	 } else {
           phi[t] = inv_logit(
                      beta_phi[1]                   + 
                      beta_timegaps  * time_gaps[t] +
                      beta_offseason * offseason[t] +	
                      beta_phi[2]    * X[ind_time_rep[t], sampling_events_phi[t]]
                    );
	 }  

	}

        // for all events estimate an individuals detection probability -- conditional on whether the individual was known to be present at that time or not
	
	for (t in 1:ind_occ) {

	 if (p_zeros[t] == 0) {
          p[t] = inv_logit(beta_p[1] + beta_p[2] * X[ind_time_rep[t], sampling_events_p[t]]);
	 } else {
          p[t] = inv_logit(beta_p[1] + beta_p[2] * X[ind_time_rep[t], sampling_events_p[t]]) * gamma[ind_time_rep[t], periods_occ[t]];
	 }

	}
	
	//// Probability of never detecting an individual again after time t

	for (i in 1:n_ind) {
	
	int cind   = n_occasions * i;	// gives the final index for individual i
        int phiind = n_occ_min1 * i;	// gives the final index for individual i

	chi[cind] = 1.0;                // on the last sampling date the probability is one for all individuals

	for (t in 1:(n_occasions - 1)) {   // loop over from the first to the second to last sampling date and multiply out the probabilities
	
	 int t_curr     = cind - t;
         int t_next     = t_curr + 1;
	 int t_curr_phi = phiind - t + 1;

	 chi[t_curr] = (1 - phi[t_curr_phi]) + 
                                      phi[t_curr_phi] * 
                                      (1 - p[t_next - 1]) * 
                                      chi[t_next];		
	}

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

print(phi[(174+2):(174+28)]);

	 for (i in 1:n_ind) {
	  for (t in (first[i] + 1):last[i]) {
	   1 ~ bernoulli(phi[(((i - 1) * n_occ_min1) + t)]);     
	  }
	
	  for (t in 1:N_y) {
   	   y[t] ~ bernoulli(p[t]);          // detection can inform y up until the last time the individual was caught, and then it is inseparable from chi
	  }

	   1 ~ bernoulli(chi[((i - 1) * n_occasions) + last[i]]);  // the probability of an animal never being seen again after the last time it was captured

	  }

}

generated quantities {
 

          
}


