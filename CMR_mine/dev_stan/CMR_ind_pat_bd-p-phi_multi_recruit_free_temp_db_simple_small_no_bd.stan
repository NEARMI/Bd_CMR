functions {

	// Function to build the required chi_t probability estimator for the mark recapture model, which
	 // is the probability of never recapturing an individual again after capturing them at time t
	
	real[] prob_uncaptured(int n_occ, real[] p_sub, real[] phi_sub) {

	real chi_sub[n_occ];    // chi for each capture date and individual 

        chi_sub[n_occ] = 1.0;         // on the last sampling date the probability is one

	for (t in 1:(n_occ - 1)) {   // loop over from the first to the second to last sampling date and multiply out the probabilities
         int t_curr = n_occ - t;
         int t_next = t_curr + 1;
	 chi_sub[t_curr] = (1 - phi_sub[t_curr]) + phi_sub[t_curr] * (1 - p_sub[t_next - 1]) * chi_sub[t_next];		
      }

  return chi_sub;

 }

}


data {

  // dimensional and bookkeeping params (single vals)
	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years)
	int<lower=1> n_times;		    	            // Sum of total time points modeled within each season across all seasons/years  	
	int<lower=1> n_occasions;		            // Total number of sampling days (a subset of times) across all seasons/years 
	int<lower=1> n_occ_min1;			    // ^^ just throwing out the very last sampling day
	
	int<lower=0> ind_occ;
	int<lower=0> ind_occ_min1;
	
  // dimensional and bookkeeping params (vectors)
	int<lower=0> time[n_times];		 	             // Vector indicating time (e.g., weeks) *!within each season!*
	
  // long vector indices
	int<lower=0> sampling_events_phi[ind_occ_min1]; 	
	int<lower=0> sampling_events_p[ind_occ];   
	int<lower=0> ind_occ_min1_rep[ind_occ_min1];
	int<lower=0> ind_occ_rep[ind_occ];
	  
  // covariates
	row_vector[n_times] X_bd[n_ind];		    // Covariate
	int<lower=0> time_gaps[ind_occ_min1];  	 	    // Elapsed time between each sampling event 

  // captures
	int<lower=1> N_y;				    // Number of defined values for captures
  	int<lower=0, upper=1> y[N_y];		            // The capture values 

  	int<lower=0> first[n_ind];         		    // Capture event in which each individual was first captured
	int<lower=0> last[n_ind];         		    // Capture event in which each individual was last captured

	int<lower=0> phi_zeros[ind_occ_min1];

}

transformed data {

} 

parameters {

	vector[2] beta_phi;                  		 // intercept and slope coefficient for survival
        vector[2] beta_p;				 // intercept and slope coefficient for detection
	real<upper=0> beta_timegaps;			 // coefficient to control for the variable time between sampling events

}

transformed parameters {

	real<lower=0,upper=1> phi[ind_occ_min1];   // survival from t to t+1, each individual repeated n_occ_min1 times 
	real<lower=0,upper=1> p[ind_occ];          // detection at time t
	real<lower=0,upper=1> chi[ind_occ];        // probability an individual will never be seen again


	//// Survival and detection probability over the whole period

	for (t in 1:ind_occ_min1) {

	 if (phi_zeros[t] == 1) {
           phi[t] = 0;
	 } else {
           phi[t] = inv_logit(
                      beta_phi[1]                   + 
                      beta_timegaps  * time_gaps[t] +
                      beta_phi[2]    * X_bd[ind_occ_min1_rep[t], sampling_events_phi[t]]
                    );
	 }  

	}

        // for all events estimate an individuals detection probability -- conditional on whether the individual was known to be present at that time or not
	
	for (t in 1:ind_occ) {
          p[t] = inv_logit(beta_p[1] + beta_p[2] * X_bd[ind_occ_rep[t], sampling_events_p[t]]);
	}

	
	//// Probability of never detecting an individual again after time t

	for (i in 1:n_ind) {
	  chi[((i - 1) * n_occasions + 1):(i * n_occasions)] = prob_uncaptured(n_occasions, segment(p, (i-1)*n_occasions + 1, 10), segment(phi, (i-1)*n_occ_min1 + 1, 9));
	}

}

model {

	// Priors

	beta_phi[1] ~ normal(0, 5);
	beta_phi[2] ~ normal(0, 5);
	beta_p[1] ~ normal(0, 5);
	beta_p[2] ~ normal(0, 5);
	beta_timegaps ~ normal(0, 5);
  

	//// Capture model

	 for (i in 1:n_ind) {
	  for (t in (first[i] + 1):last[i]) {			
	   1 ~ bernoulli(phi[(((i - 1) * n_occ_min1) + t - 1)]);    // Survival _to_ t (from phi[t - 1]) is 1 because we know the individual lived in that period 
	  }
	  for (t in 1:last[i]) {
	   y[((i - 1) * n_occasions) + t] ~ bernoulli(p[((i - 1) * n_occasions) + t]);
	  }
	   1 ~ bernoulli(chi[((i - 1) * n_occasions) + last[i]]);  // the probability of an animal never being seen again after the last time it was captured	
	 }

}

generated quantities {
 

          
}


