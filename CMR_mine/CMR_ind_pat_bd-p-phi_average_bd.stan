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

	int<lower=0> n_ind;				    // Total number of individuals caught 
	int<lower=2> n_times;				    // Number of discrete time points over the whole season
	int<lower=2> n_occasions;		            // Number of capture occasions on a subset of times
	int<lower=1> n_oc_min1;
	row_vector[n_times] time;		 	    // Vector indicating time
	int<lower=0> sampling[n_times];		    	    // Vector indicating the indices of time on which a sampling event occurred
	int<lower=0> sampling_events[n_occasions];

	int<lower=0,upper=1> y[n_ind, n_occasions];	    // Capture-history observation matrix of bd-unmeasured individuals
  	int<lower=0,upper=n_occasions> first[n_ind];        // Capture time that each individual was first captured
  	int<lower=0,upper=n_occasions> last[n_ind];         // Capture time that each individual was last captured

	row_vector[n_occasions] X_bd[n_ind];		    // Covariate
	matrix[n_ind, n_occasions] X_measured;		    // Captures during which Bd was taken
	int<lower=0,upper=n_occasions> time_gaps[n_oc_min1];
	int<lower=0,upper=n_times> bd_after_gap[n_occasions];

}

transformed data {


} 

parameters {

	vector[2] beta_phi;                  		// intercept and slope coefficient for survival
        vector[2] beta_p;				// intercept and slope coefficient for detection
	vector[3] beta_bd;				// intercept and two slope coefficients for grand mean change in bd over time
	real<lower=0> beta_timegaps;

	real<lower=0> bd_delta_sigma;			// change in Bd by individual (normal random effect variance)
	real bd_delta_eps[n_ind];			// the conditions modes of the random effect (each individual's intercept (for now))

	real<lower=0> bd_obs;    			// observation noise for observed Bd compared to underlying state

}

transformed parameters {

	// per sample * per individual mortality and detection probability
	matrix<lower=0,upper=1>[n_ind, n_oc_min1] phi;	// probability of surviving to the next period
	matrix<lower=0,upper=1>[n_ind, n_occasions] p;	// probability of being captured in a period
	matrix<lower=0,upper=1>[n_ind, n_occasions] chi;

        row_vector[n_times] X[n_ind];			// Estimated "true" bd for all of the caught individuals with no bd measured
	real bd_ind[n_ind];                             // Individual random effect deviates
	row_vector[n_occasions] X_avg[n_ind];		// Average bd load over the time period in-between measures

	for (i in 1:n_ind) {
  	  bd_ind[i]  = beta_bd[1] + bd_delta_sigma  * bd_delta_eps[i];

	for (t in 1:n_times) {
	  X[i, t] = bd_ind[i] + beta_bd[2] * time[t] + beta_bd[3] * square(time[t]);
	}

	// Average estimated bd load from each sampling day through to the next sampling day
	for (t in 1:n_occasions) {
	  X_avg[i, t] = mean(X[i, t:bd_after_gap[t]]);
	}

	}

	// Constraints on survival and detection based on knowns about the data
	for (i in 1:n_ind) {

	// for all events prior to the first catch set individuals mortality and detection to 0
	 for (t in 1:(first[i] - 1)) {
	  phi[i, t] = 0;
          p[i, t] = 0;      
	 }
		
         // linear predictor for survival for individual i and time t based on the covariate X. Prior to first capture (see above loop)
          // the probabilities are 0, after first capture the probabilities are determined by their covariates
           for (t in first[i]:n_occasions) {
         
	    p[i, t] = inv_logit(beta_p[1] + beta_p[2] * X[i, sampling_events[t]]);

	 }
	
	  for (t in first[i]:n_oc_min1){

	//  phi[i, t] = inv_logit(beta_phi[1] + beta_timegaps * time_gaps[t] + beta_phi[2] * X[i, sampling_events[t]]);
	    phi[i, t] = inv_logit(beta_phi[1] + beta_timegaps * time_gaps[t] + beta_phi[2] * X_avg[i, t]);

	  }

	}

	chi = prob_uncaptured(n_ind, n_occasions, p, phi);
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

	bd_delta_sigma ~ inv_gamma(1, 1);
	bd_delta_eps ~ normal(0, 1);
	bd_obs ~ inv_gamma(1, 1);

	// Bd Process and Data Model
    
        for (i in 1:n_ind) {
         for (t in 1:n_occasions) {

        // measured bd (X_bd) only informs latent bd (X) on those occasions where bd was sampled
           if (X_measured[i, t] == 1) {
             X_bd[i, t] ~ normal(X[i, sampling_events[t]], bd_obs); 
	   }

   	 }            
	}

	// Capture model
	for (i in 1:n_ind) {
//	 if (first[i] > 0) {
	  
	  for (t in (first[i] + 1):last[i]) {
	   1 ~ bernoulli(phi[i, t - 1]);
           y[i, t] ~ bernoulli(p[i, t - 1]);
	  }
	   1 ~ bernoulli(chi[i, last[i]]);

//	  }
	}

}

generated quantities {
  

}


