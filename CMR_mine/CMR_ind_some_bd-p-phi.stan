functions {

	// Function to build the required X_t (chi_t) probability estimator for the mark recapture model, which
	 // is the probability of never recapturing an individual again after capturing them at time t
	matrix prob_uncaptured(int n_ind, int n_occasions, matrix p, matrix phi) {

	matrix[n_ind, n_occasions] chi;    // chi for each capture date and individual 

	for (i in 1:n_ind) {
         chi[i, n_occasions] = 1.0;       // on the last sampling date the probability is one

	for (t in 1:(n_occasions - 1)) {  // loop over from the first to the second to last sampling date and multiply out the probabilities
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
	int<lower=0> n_ind_bd;				    // Total number of individuals caught WITH bd measured
	int<lower=2> n_occasions;		            // Number of capture occasions
	int n_occ_minus_1; 				    // Number of capture occasions minus 1 (create a data entry for it given that n_occasions - 1 is often used)

	int<lower=0,upper=1> y[n_ind, n_occasions];	    // Capture-history observation matrix of bd-unmeasured individuals
  	int<lower=0,upper=n_occasions> first[n_ind];        // Capture time that each individual was first captured
  	int<lower=0,upper=n_occasions> last[n_ind];         // Capture time that each individual was last captured
	vector[n_occasions] n_captured;                     // Total number captured at each sampling event

	row_vector[n_occ_minus_1] X_bd[n_ind_bd];	    // Covariate
	int<lower=0> X_which[n_ind_bd];		   	    // Individuals for which Bd was taken

//	int<lower=1,upper=n_ind> ind_id[n_ind];		    // Individual ID

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

        row_vector[n_occ_minus_1] X[n_ind];		// Estimated "true" bd for all of the caught individuals with no bd measured

//	real<lower=0> sigma_alpha_phi;			// random effect variance in the intercept for survival
//	real<lower=0> sigma_beta_phi;			// random effect variance in the slope for survival over time
//	vector[n_ind] eps_alpha_phi;			// conditional modes of the random effect
//	vector[n_ind] eps_beta_phi;

//	real<lower=0> sigma_alpha_p;			// random effect variance in the intercept for detection
//	real<lower=0> sigma_beta_p;			// random effect variance in the slope for detection over time
//	vector[n_ind] eps_alpha_p;			// conditional modes of the random effect
//	vector[n_ind] eps_beta_p;

}

transformed parameters {

	// per sample * per individual mortality and detection probability
	matrix<lower=0,upper=1>[n_ind, n_occ_minus_1] phi;
	matrix<lower=0,upper=1>[n_ind, n_occ_minus_1] p;
	matrix<lower=0,upper=1>[n_ind, n_occasions] chi;
	real bd_ind[n_ind];                                    // Individual random effect deviates for individual bd change

	// Individual random effect deviates for individual bd change
	for (i in 1:n_ind) {
  	  bd_ind[i] = bd_delta_mu + bd_delta_sigma * bd_delta_eps[i];
	}
	
	// Constraints given when each individual was first captured
	for (i in 1:n_ind) {

       	// For all of those individuals that were caught at least once, set these individuals' capture and mortality 
         // probabilities prior to being captured the first time to 0  
	for (t in 1:(first[i] - 1)) {
	  if (first[i] > 0) {
	   phi[i, t] = 0;
           p[i, t] = 0;
	  }

	 }
		
         // linear predictor for survival for individual i and time t based on the covariate X
          // I am a bit unclear here what should be looped over -- all times or just the times post first capture
          // That is, I am a bit unsure how the current model is handling individuals that were never caught (which I guess shouldn't be in the data at all)
      // for (t in first[i]:n_occ_minus_1) { 
         for (t in 1:n_occ_minus_1) {

//	  mu = inv_logit((beta[1] + eps_alpha_phi[ind_id[i]] * sigma_alpha_phi) + 
//		(beta[2] + eps_beta_phi[ind_id[i]] * sigma_beta_phi) * X[i, t]);
          
	// for the random effect add the eps_alpha_phi and eps_alpha_p below (see above)
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

//	sigma_alpha_phi ~ inv_gamma(1, 1);
//	sigma_beta_phi ~ inv_gamma(1, 1);
//	eps_alpha_phi ~ normal(0, 1);
//      eps_beta_phi ~ normal(0, 1);

//	sigma_alpha_p ~ inv_gamma(1, 1);
//	sigma_beta_p ~ inv_gamma(1, 1);
//	eps_alpha_p ~ normal(0, 1);
//      eps_beta_p ~ normal(0, 1);
		
	// Bd Process Model

        for (i in 1:n_ind) {
	 X[i, 1] ~ normal(start_mean, start_var);
	}

        for (t in 2:n_occ_minus_1) {
	 for (j in 1:n_ind) {
          X[j, t] ~ normal(X[j, t-1] + bd_ind[j], bd_add);  
	 }
	}  

	// Bd Data Model
         // Slightly confusing way to write this but works because order is conserved. That is,
         // the length(X_which) individuals of the full X line up with the subset X_bd
        for (t in 2:n_occ_minus_1) {
          X_bd[, t] ~ normal(X[X_which, t], bd_obs); 
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
        
      vector<lower=0>[n_occ_minus_1] pop;      // Estimate of the full population
      vector<lower=0>[n_occ_minus_1] captures; // Estimate of the number captured at each sampling event
      vector<lower=0>[n_occ_minus_1] avg_p;    // Average detection probability on each sampling period

      for (k in 1:n_occ_minus_1){
       avg_p[k] = mean(p[, k]);
      }
      pop = n_captured[1:n_occ_minus_1] ./ avg_p;
      captures = pop[1:n_occ_minus_1] .* avg_p;
   
}


