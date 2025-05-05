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

	int<lower=1> n_periods;				    // Total number of seasons/years (the "on" period where sampling occurs) over which individuals are captured
	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years)
	int<lower=n_periods> n_times;		    	    // Sum of total time points modeled within each season across all seasons/years  		
	int<lower=n_periods> n_occasions;		    // Total number of sampling days (a subset of times) across all seasons/years 
	int<lower=1> n_occ_min1;			    // ^^ just throwing out the very last sampling day
	
	int<lower=0> time[n_times];		 	    // Vector indicating time (e.g., weeks) *!within each season!*

	int<lower=0> sampling_events[n_occasions]; 	    // Indices from 1:n_times on which a sampling event occurred
	int<lower=0> time_gaps[n_occ_min1];  	 	    // Elapsed time between each sampling event 

	int<lower=0,upper=1> offseason[n_occ_min1];	    // Vector indicating the last sampling periods of each season which gains offseason characteristics
	  
	int y[n_ind, n_occasions];		    	    // Capture-history observation matrix of bd-unmeasured individuals
  	int<lower=0> first[n_ind];         		    // Capture event in which each individual was first captured
  	int<lower=0> last[n_ind];         		    // Capture event in which each individual was last captured

	matrix[n_ind, n_occasions] X_bd;	   	    // Covariate
	matrix[n_ind, n_occasions] X_measured;    	    // Captures during which Bd was taken
	int<lower=0> periods[n_times];			    // Vector designating periods

}

transformed data {


} 

parameters {

	vector[2] beta_phi;                  		// intercept and slope coefficient for survival
        vector[2] beta_p;				// intercept and slope coefficient for detection
	real beta_timegaps;				// coefficient to control for the variable time between sampling events
	real beta_offseason;				// season survival probability

	vector[3] beta_bd;				// intercept and two slope coefficients for grand mean change in bd over time
	real beta_period;				// difference in load between primary n_periods (years)

	real<lower=0> bd_delta_sigma;			// change in Bd by individual (normal random effect variance)
	real bd_delta_eps[n_ind];			// the conditions modes of the random effect (each individual's intercept (for now))

	real<lower=0> bd_obs;    			// observation noise for observed Bd compared to underlying state

}

transformed parameters {

	matrix<lower=0,upper=1>[n_ind, n_occ_min1] phi;    // survival from t to t+1
	matrix<lower=0,upper=1>[n_ind, n_occasions] p;     // detection at time t
	matrix<lower=0,upper=1>[n_ind, n_occasions] chi;   // probability an individual will never be seen again
 
	matrix[n_ind, n_times] X;	   		   // Estimated "true" bd for all of the caught individuals with no bd measured		
	real bd_ind[n_ind];                                // Individual random effect deviates


	// bd submodel, contained to estimating within-season bd

	for (i in 1:n_ind) {
  	  bd_ind[i]  = beta_bd[1] + bd_delta_sigma  * bd_delta_eps[i];

	 for (t in 1:n_times) {
	  X[i, t] = bd_ind[i] + beta_period * periods[t] + beta_bd[2] * time[t] + beta_bd[3] * square(time[t]);
	 }
	}

	// Survival and probability over the whole period

	for (i in 1:n_ind) {

	// for all events prior to the first catch set individuals mortality and detection to 0
	 for (t in 1:(first[i] - 1)) {
	  phi[i, t] = 0;
          p[i, t] = 0;      
	 }
		
         // linear predictor for survival for individual i and time t based on the covariate X. Prior to first capture (see above loop)
          // the probabilities are 0, after first capture the probabilities are determined by their covariates
           for (t in first[i]:n_occ_min1) {
          
	  phi[i, t] = inv_logit(
beta_phi[1]                   + 
beta_timegaps  * time_gaps[t] +
beta_offseason * offseason[t] +	
beta_phi[2]    * X[i, sampling_events[t]]
);

	}

          for (t in first[i]:n_occasions) {

	  p[i, t]   = inv_logit(
beta_p[1] + 
beta_p[2] * X[i, sampling_events[t]]
);

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
	beta_period ~ normal(0, 5);
	beta_timegaps ~ normal(0, 5);
	beta_offseason ~ normal(0, 5);

	bd_delta_sigma ~ inv_gamma(1, 1);
	bd_delta_eps ~ normal(0, 2);
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
         // Is there a way to add an additional process here that is the probability the individual was in the population?
         // i.e., pres[n_ind, n_times] 
         // 1 ~ bernoulli(pres[i, t])
         // where pres is informed by _something_

	 for (i in 1:n_ind) {
	  if (first[i] > 0) {
	  
	  for (t in (first[i] + 1):last[i]) {
	   1 ~ bernoulli(phi[i, t - 1]);     // survived (the one) ~ bernoulli[phi] -- because we saw an individual again in the future we know it survived this period
           y[i, t] ~ bernoulli(p[i, t]);
	  }
	   1 ~ bernoulli(chi[i, last[i]]);  // the probability of an animal never being seen again after the last time it was captured

	  }
	 }

}

generated quantities {
           
}


