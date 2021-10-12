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

  // dimensional and bookkeeping params (single vals)
	int<lower=0> n_ind;				    // Total number of individuals caught 
	int<lower=2> n_times;				    // Number of discrete time points over the whole season
	int<lower=2> n_occasions;		            // Number of capture occasions on a subset of times
	int<lower=1> n_oc_min1;

  // dimensional and bookkeeping params (vectors)	
	row_vector[n_times] time;		 	    // Vector indicating time
	int<lower=0> sampling[n_times];		    	    // Vector indicating the indices of time on which a sampling event occurred
	int<lower=0> sampling_events[n_occasions];

	row_vector[n_times] X_bd[n_ind];		    // Covariate
	int<lower=0,upper=n_occasions> time_gaps[n_oc_min1];

	int<lower=0,upper=1> y[n_ind, n_occasions];	    // Capture-history observation matrix of bd-unmeasured individuals
  	int<lower=0,upper=n_occasions> first[n_ind];        // Capture time that each individual was first captured
  	int<lower=0,upper=n_occasions> last[n_ind];         // Capture time that each individual was last captured

}

transformed data {


} 

parameters {

	vector[2] beta_phi;                  		// intercept and slope coefficient for survival
        vector[2] beta_p;				// intercept and slope coefficient for detection
	real beta_timegaps;

}

transformed parameters {

	// per sample * per individual mortality and detection probability
	matrix<lower=0,upper=1>[n_ind, n_oc_min1] phi;	// probability of surviving to the next period
	matrix<lower=0,upper=1>[n_ind, n_occasions] p;	// probability of being captured in a period
	matrix<lower=0,upper=1>[n_ind, n_occasions] chi;

	// Constraints on survival and detection based on knowns about the data
	for (i in 1:n_ind) {

	// for all events prior to the first catch set individuals mortality and detection to 0
	 for (t in 1:(first[i] - 1)) {
	  phi[i, t] = 0;     
	 }
		
         // linear predictor for survival for individual i and time t based on the covariate X. Prior to first capture (see above loop)
          // the probabilities are 0, after first capture the probabilities are determined by their covariates
           for (t in 1:n_occasions) {
	    p[i, t] = inv_logit(beta_p[1] + beta_p[2] * X_bd[i, sampling_events[t]]);
	  }
	
	  for (t in first[i]:n_oc_min1){
	    phi[i, t] = inv_logit(beta_phi[1] + beta_timegaps * time_gaps[t] + beta_phi[2] * X_bd[i, sampling_events[t]]);
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
	beta_timegaps ~ normal(0, 5);


	// Capture model
	for (i in 1:n_ind) {
	  
	  for (t in (first[i] + 1):last[i]) {
	   1 ~ bernoulli(phi[i, t - 1]);
           y[i, t] ~ bernoulli(p[i, t]);
	  }
	   1 ~ bernoulli(chi[i, last[i]]);

	}

}

generated quantities {
  

}


