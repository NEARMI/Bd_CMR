functions {

	// Function to build the required X_t (chi_t) probability estimator for the mark recapture model, which
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


// Important note about this model is that it seems reasonable set up to answer questions about the effects of bd on survival, but because
 // it is a regression-form and not dynamics, it may be less suited for a number of other things

// With this updated model (SEP 22 NOTE the vector "time" is used to model the latent bd process and the vector "sampling" is used to designate
 // a time on which sampling occurred at a specific location

	int<lower=0> n_ind;				    // Total number of individuals caught 
	int<lower=2> n_times;				    // Number of discrete time points over the whole season
	int<lower=2> n_occasions;		            // Number of capture occasions on a subset of times
	row_vector[n_times] time;		 	    // Vector indicating time
	int<lower=1> sampling[n_occasions];		    // Vector indicating the indices of time on which a sampling event occurred

	int<lower=0,upper=1> y[n_ind, n_occasions];	    // Capture-history observation matrix of bd-unmeasured individuals
  	int<lower=0,upper=n_occasions> first[n_ind];        // Capture time that each individual was first captured
  	int<lower=0,upper=n_occasions> last[n_ind];         // Capture time that each individual was last captured
	vector[n_occasions] n_captured;                     // Total number captured at each sampling event

	row_vector[n_occasions] X_bd[n_ind];		    // Covariate
	matrix[n_ind, n_occasions] X_measured;		    // Captures during which Bd was taken

}

transformed data {

} 

parameters {

	vector[2] beta_phi;                  		// intercept and slope coefficient for survival
        vector[2] beta_p;				// intercept and slope coefficient for detection
	vector[3] beta_bd;				// intercept and two slope coefficients for grand mean change in bd over time
	
	real<lower=0> bd_obs;    			// observation noise for observed Bd compared to underlying state
  	
}

transformed parameters {

	// per sample * per individual mortality and detection probability
	matrix<lower=0,upper=1>[n_ind, n_occasions] phi;
	matrix<lower=0,upper=1>[n_ind, n_occasions] p;
	matrix<lower=0,upper=1>[n_ind, n_occasions] chi;

	row_vector[n_times] mu;
        row_vector[n_times] X[n_ind];			// Estimated "true" bd for all of the caught individuals with no bd measured
   

	// Constraints given when each individual was first captured
	for (i in 1:n_ind) {

       	// Because we know individuals did not die or were not captured in advance of their first capture, set these to 0
         // Though TBH I find this a bit confusing because these p and phi are based on latent covariates that are modeled, so why do we set them to 0?
	for (t in 1:(first[i] - 1)) {
	  phi[i, t] = 0;
          p[i, t] = 0;      
	 }
		
         // linear predictor for survival for individual i and time t based on the covariate X. Prior to first capture (see above loop)
          // the probabilities are 0, after first capture the probabilities are determined by their covariates
          
	 for (t in first[i]:n_occasions) {
          
	// the individual survival and capture probabilities are based on an individuals nature (random effect of individual for intercept)
         // and a grand mean relationship between bd and their bd (which does vary by individual)

	  phi[i, t] = inv_logit(beta_phi[1] + beta_phi[2] * X[i, sampling[t]]);
	  p[i, t] = inv_logit(beta_p[1] + beta_p[2] * X[i, sampling[t]]);

	 }

	}


	chi = prob_uncaptured(n_ind, n_occasions, p, phi);

	mu = beta_bd[1] + beta_bd[2] * time + beta_bd[3] * square(time);

        for (i in 1:n_ind) {
	  for (j in 1:n_times) {
		
	    X[i, j] = mu[j];

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

	bd_obs ~ inv_gamma(1, 1);


	// Bd Data Model
    
        for (i in 1:n_ind) {
         for (t in 1:n_occasions) {

        // measured bd (X_bd) only informs latent bd (X) on those occasions where bd was sampled
           if (X_measured[i, t] == 1) {
             X_bd[i, t] ~ normal(X[i, sampling[t]], bd_obs); 
	   }

   	 }            
	}

	// Capture model
	for (i in 1:n_ind) {
	 if (first[i] > 0) {
	  
	  for (t in (first[i] + 1):last[i]) {
	   1 ~ bernoulli(phi[i, t - 1]);
           y[i, t] ~ bernoulli(p[i, t - 1]);
	  }
	   1 ~ bernoulli(chi[i, last[i]]);

	  }
	}

}

generated quantities {
        
 //     vector<lower=0>[n_occasions] pop;      // Estimate of the full population
 //     vector<lower=0>[n_occasions] captures; // Estimate of the number captured at each sampling event
 //     vector<lower=0>[n_occasions] avg_p;    // Average detection probability on each sampling period
 //     vector[n_occasions] total_bd;          // Total Bd load in the population

 //     for (k in 1:n_occasions){
 //      avg_p[k] = mean(p[, k]);
 //     }

 //     pop = n_captured[1:n_occasions] ./ avg_p;
 //     captures = pop[1:n_occasions] .* avg_p;

  //    for (k in 1:n_occasions){
  //     total_bd[k] = pop[k] * mean(X[, k]);
  //    }
      
   
}


