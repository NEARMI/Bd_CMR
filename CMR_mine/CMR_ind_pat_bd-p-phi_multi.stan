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

	int<lower=0> n_periods;				    // Total number of primary n_periods over which individuals are captured
	int<lower=0> n_ind;				    // Total number of individuals caught (ever, over all years)
	int<lower=2> n_times[n_periods];		    // Number of discrete time points in each season
	int<lower=2> n_occasions[n_periods];		    // Number of capture occasions on a subset of times
	int<lower=1> n_oc_min1[n_periods];		
	
	int<lower=0> time[n_times[1]];		 	    		  // Vector indicating time
	int<lower=0> sampling[n_times[1], n_periods]; 	    		  // Length of n_times; 1 if sampling, 0 if no sampling occurred
	int<lower=0> sampling_events[n_occasions[1], n_periods]; 	  // Indices of times on which a sampling event occurred
	int<lower=0> time_gaps[n_oc_min1[1], n_periods];  	 	  // Number of time n_periods in-between sampling events
	  
	int y[n_ind, n_occasions[1], n_periods];		          // Capture-history observation matrix of bd-unmeasured individuals
  	int<lower=0> first[n_ind, n_periods];         			  // Capture time that each individual was first captured
  	int<lower=0> last[n_ind, n_periods];         			  // Capture time that each individual was last captured

	matrix[n_ind, n_occasions[1]] X_bd[n_periods];	   		  // Covariate
	matrix[n_ind, n_occasions[1]] X_measured[n_periods];    	  // Captures during which Bd was taken
	matrix[n_ind, n_occasions[1]] periods[n_periods];

}

transformed data {


} 

parameters {

	vector[2] beta_phi;                  		// intercept and slope coefficient for survival
        vector[2] beta_p;				// intercept and slope coefficient for detection
	real beta_timegaps;				// coefficient to control for the variable time between sampling events

	vector[3] beta_bd;				// intercept and two slope coefficients for grand mean change in bd over time
	real beta_period;				// difference in load between primary n_periods (years)

	real<lower=0> bd_delta_sigma;			// change in Bd by individual (normal random effect variance)
	real bd_delta_eps[n_ind];			// the conditions modes of the random effect (each individual's intercept (for now))

	real<lower=0> bd_obs;    			// observation noise for observed Bd compared to underlying state

}

transformed parameters {

	// per sample * per individual mortality and detection probability
	matrix<lower=0,upper=1>[n_ind, n_oc_min1[1]] phi[n_periods];
	matrix<lower=0,upper=1>[n_ind, n_oc_min1[1]] p[n_periods];
	matrix<lower=0,upper=1>[n_ind, n_occasions[1]] chi[n_periods];

	matrix[n_ind, n_times[1]] X[n_periods];	   	// Estimated "true" bd for all of the caught individuals with no bd measured		
	real bd_ind[n_ind];                             // Individual random effect deviates

	for (i in 1:n_ind) {
  	  bd_ind[i]  = beta_bd[1] + bd_delta_sigma  * bd_delta_eps[i];

	for (k in 1:n_periods) {
	 for (t in 1:n_times[1]) {
	  X[k][i, t] = bd_ind[i] + beta_period * periods[k][i, t] + beta_bd[2] * time[t] + beta_bd[3] * square(time[t]);
	  }
	 }
	}

	// Constraints on survival and detection based on knowns about the data
	for (i in 1:n_ind) {
	 for (k in 1:n_periods) {

	// for all events prior to the first catch set individuals mortality and detection to 0
	 for (t in 1:(first[i, k] - 1)) {
	  phi[k][i, t] = 0;
          p[k][i, t] = 0;      
	 }
		
         // linear predictor for survival for individual i and time t based on the covariate X. Prior to first capture (see above loop)
          // the probabilities are 0, after first capture the probabilities are determined by their covariates
           for (t in first[i, k]:n_oc_min1[1]) {
          
	  phi[k][i, t] = inv_logit(beta_phi[1] + beta_timegaps * time_gaps[t, k] * beta_phi[2] * X[k][i, sampling_events[t, k]]);
	  p[k][i, t] = inv_logit(beta_p[1] + beta_p[2] * X[k][i, sampling_events[t, k]]);

	 }

	 }
	}

       for (k in 1:n_periods) {
	chi[k] = prob_uncaptured(n_ind, n_occasions[k], p[k], phi[k]);
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
	beta_period ~ normal(0, 5);
	beta_timegaps ~ normal(0, 5);

	bd_delta_sigma ~ inv_gamma(1, 1);
	bd_delta_eps ~ normal(0, 1);
	bd_obs ~ inv_gamma(1, 1);


	// Bd Process and Data Model
    
	for (k in 1:n_periods) {
         for (i in 1:n_ind) {
          for (t in 1:n_occasions[1]) {

        // measured bd (X_bd) only informs latent bd (X) on those occasions where bd was sampled
           if (X_measured[k][i, t] == 1) {
             X_bd[k][i, t] ~ normal(X[k][i, sampling_events[t, k]], bd_obs); 
	   }

   	  }            
	 }
	}

	// Capture model
	for (k in 1:n_periods) {
	 for (i in 1:n_ind) {
	  if (first[i, k] > 0) {
	  
	  for (t in (first[i, k] + 1):last[i, k]) {
	   1 ~ bernoulli(phi[k][i, t - 1]);
           y[i, t, k] ~ bernoulli(p[k][i, t - 1]);
	  }
	   1 ~ bernoulli(chi[k][i, last[i]]);

	  }
	 }
	}

}

generated quantities {
           
}


