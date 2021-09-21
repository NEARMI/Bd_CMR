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

	int<lower=0> n_ind;				    // Total number of individuals caught 
	int<lower=2> n_occasions;		            // Number of capture occasions

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

	real bd_delta_mu;				// grand mean change in Bd
	real<lower=0> bd_delta_sigma;			// change in Bd by individual (normal random effect variance)
	real bd_delta_eps[n_ind];			// the conditions modes of the random effect (each individual's slope in bd over time)
	real<lower=0> bd_add;    			// additive noise in changes in Bd over time
	real<lower=0> bd_obs;    			// observation noise for observed Bd compared to underlying state
  
        real start_mean;
	real<lower=0> start_var; 

        row_vector[n_occasions] X[n_ind];		// Estimated "true" bd for all of the caught individuals with no bd measured

}

transformed parameters {

	// per sample * per individual mortality and detection probability
	matrix<lower=0,upper=1>[n_ind, n_occasions] phi;
	matrix<lower=0,upper=1>[n_ind, n_occasions] p;
	matrix<lower=0,upper=1>[n_ind, n_occasions] chi;

	real bd_ind[n_ind];                                    // Individual random effect deviates

	// Individual random effect deviates for individual bd change and just general random effects for survival and detection
         // for now no random slopes, just random intercepts (each individual has some characteristic about them that makes them more
          // likely to die or to be detected)
	for (i in 1:n_ind) {
  	  bd_ind[i]  = bd_delta_mu + bd_delta_sigma  * bd_delta_eps[i];

	}
	
	// Constraints given when each individual was first captured
	for (i in 1:n_ind) {

       	// Because we know individuals did not die or were not captured in advance of their first capture, set these to 0
	for (t in 1:(first[i] - 1)) {
	  phi[i, t] = 0;
          p[i, t] = 0;      
	 }
		
         // linear predictor for survival for individual i and time t based on the covariate X. Prior to first capture (see above loop)
          // the probabilities are 0, after first capture the probabilities are determined by their covariates
           for (t in first[i]:n_occasions) {
          
	  phi[i, t] = inv_logit(beta_phi[1] + beta_phi[2] * X[i, t]);
	  p[i, t] = inv_logit(beta_p[1] + beta_p[2] * X[i, t]);

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

	bd_delta_mu ~ normal(0, 4);
	bd_delta_sigma ~ inv_gamma(1, 1);
	bd_delta_eps ~ normal(0, 1);
	bd_add ~ inv_gamma(1, 1);
	bd_obs ~ inv_gamma(1, 1);

	start_mean ~ normal(5, 2);
	start_var ~ inv_gamma(1, 1);


	// Bd Process Model

        for (i in 1:n_ind) {
	 X[i, 1] ~ normal(start_mean, start_var);
	}

        for (t in 2:n_occasions) {
	 for (j in 1:n_ind) {

        // individuals vary in their bd slopes over time 
          X[j, t] ~ normal(X[j, t-1] + bd_ind[j], bd_add);  

	 }
	}  

	// Bd Data Model
    
        for (t in 1:n_occasions) {
          for (i in 1:n_ind) {

        // measured bd (X_bd) only informs latent bd (X) on those occasions where bd was sampled
           if (X_measured[i, t] == 1) {
          X_bd[i, t] ~ normal(X[i, t], bd_obs); 
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


