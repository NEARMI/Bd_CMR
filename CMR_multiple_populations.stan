functions {

	// Function to build the required chi_t probability estimator for the mark recapture model, which
	 // is the probability of never recapturing an individual again after capturing them at time t
	// rewritten to run for an individual at a time given variable numbers of capture opportunities by individual (e.g. in different populations)
	
	real[] prob_uncaptured(int n_occ, real[] p_sub, real[] phi_sub) {

	real chi_sub[n_occ];          // chi for each capture date and individual 

        chi_sub[n_occ] = 1.0;         // on the last sampling date the probability is one

	for (t in 1:(n_occ - 1)) {    // loop over from the first to the second to last sampling date and multiply out the probabilities
         int t_curr = n_occ - t;
         int t_next = t_curr + 1;
	 chi_sub[t_curr] = (1 - phi_sub[t_curr]) + phi_sub[t_curr] * (1 - p_sub[t_next]) * chi_sub[t_next];		
      }

  return chi_sub;

 }

}


data {

// -----
// Notes
// -----
// Given long comments, this model is best read in full screen on an external monitor


  // dimensional and bookkeeping params (single vals)
	int<lower=1> n_pop;				    // Number of distinct populations (sampling areas)
	int<lower=1> n_pop_year;
	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years)
	int<lower=1> n_spec;	
	int<lower=1> ind_per_period_p;			   
	int<lower=1> ind_per_period_bd;
	int<lower=0> ind_occ;			   	    // n_ind * n_occasions, summed over the sampling of all populations
	int<lower=0> ind_occ_min1;		 	    // n_ind * n_occ_min1, summed over the sampling of all populations
	int<lower=0> ind_time;	
	
  // dimensional and bookkeeping params (vectors)	
	int<lower=1> ind_occ_size[n_ind];		    // Number of sampling periods for all individuals
	int<lower=1> ind_occ_min1_size[n_ind];		    // Number of sampling periods -1 for all individuals
	int<lower=1> phi_first_index[n_ind];		    // The indexes of phi corresponding to the first entry for each individual
	int<lower=1> p_first_index[n_ind];	            // The indexes of p corresponding to the first entry for each individual
	int<lower=0> ind_in_pop[n_ind];
	
  // long vector indices for observation model (p)
	int<lower=0> ind_occ_rep[ind_occ];		    // Index vector of all individuals (each individual repeated the number of sampling occasions)
	int<lower=0> periods_occ[ind_occ];		    // Vector designating periods for observational model (all occasions)
	int<lower=0> p_month[ind_occ];	
	int<lower=0> p_zeros[ind_occ];			    // Observation times for each individual in which we do not know if that individual is present
	int<lower=0> pop_p[ind_occ];			    // population index for detection predictors
	int<lower=0> spec_p[ind_occ];
	int<lower=0> p_bd_index[ind_occ];		    // which entries of latent bd correspond to each entry of p
	int<lower=1> gamma_index[ind_occ];		    // gamma value associated with each entry of p
  
  // long vector indices for survival model (phi)
	int<lower=0> ind_occ_min1_rep[ind_occ_min1];	    // Index vector of all individuals (each individual repeated the number of sampling occasions -1)
	int<lower=0, upper=1> offseason[ind_occ_min1];	    // Vector indicating the last sampling periods of each season which gains offseason characteristics
	int<lower=0> phi_year[ind_occ_min1];			
	int<lower=0> phi_zeros[ind_occ_min1];		    // Observation times for each individual in advance of first detecting that individual
	int<lower=0> phi_ones[ind_occ_min1];
	int<lower=0> spec_phi[ind_occ_min1];
	int<lower=0> pop_phi[ind_occ_min1];		    // population index for mortality predictors
	int<lower=0> phi_bd_index[ind_occ_min1];	    // which entries of the summarized stat correspond to each period

  // long vector indices for bd (bd)
	int<lower=0> ind_bd_rep[ind_per_period_bd];
	int<lower=0> bd_time[ind_time];
	int<lower=0> ind_in_pop_year[ind_per_period_bd];
	  
  // covariates
	int<lower=1> N_bd;				    // Number of defined values for bd
 	real X_bd[N_bd];			   	    // The bd values 
	int<lower=1> X_ind[N_bd];			    // Individual associated with each bd measure
	int<lower=0> x_bd_index[N_bd];			    // entries of X (latent bd) that have a corresponding real measure to inform likelihood with

	int<lower=0> bd_first_index[ind_per_period_bd];	    // First entry of latent bd associated with each individual 'by' period
	int<lower=0> bd_last_index[ind_per_period_bd];	    // Last entry of latent bd associated with each individual 'by' period

	real ind_size[n_ind];

  // captures
	int<lower=1> N_y;				    // Number of defined values for captures
  	int<lower=0, upper=1> y[N_y];		            // The capture values 

  	int<lower=0> first[n_ind];         		    // Capture event in which each individual was first captured
  	int<lower=0> last[n_ind];         		    // Capture event in which each individual was last captured

}


transformed data {

} 


parameters {

// -----
// bd submodel
// -----

	real beta_bd;
	
	vector<lower=0>[n_pop] bd_ind_sigma;		 // change in Bd by individual (normal random effect variance)		 
	real bd_ind_eps[n_ind];                          // the conditions modes of the random effect (each individual's intercept (for now))

	real<lower=0> bd_pop_sigma;			 // change in Bd by pop (normal random effect variance)
	real bd_pop_eps[n_pop_year];			 // the conditions modes of the random effect (each populations intercept (for now))

	real<lower=0> bd_obs;    			 // observation noise for observed Bd compared to underlying state	

// -----
// survival
// -----

	vector[3] beta_phi;                  		 // survival between seasons as a function of bd
	vector[2] beta_offseason;  			 // survival as a function of bd stress
	vector[3] beta_offseason_year;			 // each year has its own survival prob
	vector[6] beta_spec;

	real<lower=0> inseason_pop_sigma;		 // change in Bd by individual (normal random effect variance)
	real inseason_pop_eps[n_pop];

	real<lower=0> offseason_pop_sigma;		 // change in Bd by individual (normal random effect variance)
	real offseason_pop_eps[n_pop];

// -----
// detection
// -----

	vector[3] beta_p;				 // intercept and slope coefficient for detection

	real<lower=0> p_pop_sigma;			 // change in Bd by individual (normal random effect variance)
	real p_pop_eps[n_pop];
	
// -----
// other
// -----	

	vector<lower=0,upper=1>[ind_per_period_p] gamma;

}


transformed parameters {

	real<lower=0,upper=1> phi[ind_occ_min1];   // survival from t to t+1, each individual repeated the number of times its population was measured
	real<lower=0,upper=1> p[ind_occ];          // detection at time t
	real<lower=0,upper=1> chi[ind_time];       // probability an individual will never be seen again

	real X[ind_per_period_bd];		   // each individual's estimated bd per year 

	// individual random effect deviates
	real bd_ind[n_ind];

	// population random effect deviates     
	real bd_pop_year[n_pop_year];
	real inseason_pop[n_pop];
	real offseason_pop[n_pop];
	real p_pop[n_pop];          
	

// -----
// bd submodel, contained to estimating within-season bd
// -----

	for (pp in 1:n_pop_year) {
	  bd_pop_year[pp]  = bd_pop_sigma * bd_pop_eps[pp];
	} 

	for (i in 1:n_ind) {
	    
		// linear predictor for intercept for bd-response. Overall intercept + pop-specific intercept + individual random effect deviate
  	  bd_ind[i]  = bd_ind_sigma[ind_in_pop[i]] * bd_ind_eps[i];  

	}

	for (t in 1:ind_per_period_bd) {

		// latent bd model before obs error
	  X[t] = (beta_bd + bd_ind[ind_bd_rep[t]] + bd_pop_year[ind_in_pop_year[t]]);      

        }


// -----
// Survival probability over the whole period
// -----

	for (pp in 1:n_pop) {
	 inseason_pop[pp]  = inseason_pop_sigma  * inseason_pop_eps[pp];
	 offseason_pop[pp] = offseason_pop_sigma * offseason_pop_eps[pp];
	} 

	for (t in 1:ind_occ_min1) {

	 if (phi_zeros[t] == 1) {			// phi_zeros is 1 before an individual is caught for the first time
           phi[t] = 0;					// must be non-na values in stan, but the likelihood is only informed from first capture onward
	 } else {

	  if (phi_ones[t] == 1) {
           phi[t] = 1;
	  } else {

	   if (offseason[t] == 0) {

             phi[t] = inv_logit(beta_phi[1] + inseason_pop[pop_phi[t]] +
beta_phi[2] * X[phi_bd_index[t]] + 
beta_phi[3] * ind_size[ind_occ_min1_rep[t]]);

	   } else {

	     phi[t] = inv_logit(beta_offseason_year[phi_year[t]] + 
(beta_offseason[1] + offseason_pop[pop_phi[t]] + beta_spec[spec_phi[t]]) * X[phi_bd_index[t]] + 
beta_offseason[2]  * ind_size[ind_occ_min1_rep[t]]);
	
	   }

	   }

	  }  

	 }

// -----
// Detection probability over the whole period
// -----
	
	for (pp in 1:n_pop) {
	  p_pop[pp]   = p_pop_sigma * p_pop_eps[pp];
	}

	for (t in 1:ind_occ) {

		// p_zeros is = 1 in each season prior to an individual being caught for the first time
		// p gets scaled in these years in an attempt to scale the probability as a function of bd given that we don't know if the individual was there

	 if (p_zeros[t] == 1) {				
          p[t] = inv_logit(beta_p[1] + p_pop[pop_p[t]] +
beta_p[2] * X[p_bd_index[t]] + 
beta_p[3] * ind_size[ind_occ_rep[t]]);
	 } else {
          p[t] = inv_logit(beta_p[1] + p_pop[pop_p[t]] +
beta_p[2] * X[p_bd_index[t]] + 
beta_p[3] * ind_size[ind_occ_rep[t]]) * gamma[gamma_index[t]];
	 }

	}
	
// -----
// Probability of never detecting an individual again after time t
// -----

		// For each individual calculate the probability it won't be captured again

	for (i in 1:n_ind) {
	 chi[p_first_index[i]:(p_first_index[i] + ind_occ_size[i] - 1)] = prob_uncaptured(ind_occ_size[i], 
              segment(p, p_first_index[i], ind_occ_size[i]), segment(phi, phi_first_index[i], ind_occ_min1_size[i]));
	}

}

model {

// -----
// Priors
// -----

	beta_bd ~ normal(0, 3);

	beta_p[1]   ~ normal(0, 1.2);
	beta_p[2]   ~ normal(0, 1.2);
	beta_p[3]   ~ normal(0, 1.2);

	beta_phi[1] ~ normal(0, 1.0);
	beta_phi[2] ~ normal(0, 1.0);
	beta_phi[3] ~ normal(0, 1.0);

	beta_offseason[1] ~ normal(0, 1.0);
	beta_offseason[2] ~ normal(0, 1.0);

	beta_offseason_year ~ normal(0, 1.0);

	bd_obs         ~ inv_gamma(10, 4);

	bd_pop_sigma        ~ inv_gamma(8, 15);
	offseason_pop_sigma ~ inv_gamma(8, 15);
	inseason_pop_sigma  ~ inv_gamma(8, 15);
	p_pop_sigma         ~ inv_gamma(8, 15);

	for (i in 1:n_ind) {
	  bd_ind_eps[i] ~ normal(0, 3);
	}

	for (i in 1:n_pop) {
	  bd_ind_sigma[i]      ~ inv_gamma(8, 15);
	  inseason_pop_eps[i]  ~ normal(0, 3);
	  offseason_pop_eps[i] ~ normal(0, 3);
	  p_pop_eps[i]         ~ normal(0, 3);
	}

	for (i in 1:n_pop_year) {
	  bd_pop_eps[i]        ~ normal(0, 3);
	}
         
        gamma ~ uniform(0, 1);


// -----
// Bd Process and Data Model
// -----

		// observed bd is the linear predictor + some observation noise

	for (t in 1:N_bd) {
          X_bd[t] ~ normal(X[x_bd_index[t]], bd_obs); 
	} 
    
// -----
// Capture model
// -----

	 for (i in 1:n_ind) {
	  for (t in (first[i] + 1):last[i]) {			
	   1 ~ bernoulli(phi[phi_first_index[i] - 1 + t - 1]);    		   // Survival _to_ t (from phi[t - 1]) is 1 because we know the individual lived in that period 
	  }
	  for (t in 1:last[i]) {
	   y[p_first_index[i] - 1 + t] ~ bernoulli(p[p_first_index[i] - 1 + t]);   // Capture given detection
	  }

	   1 ~ bernoulli(chi[p_first_index[i] - 1 + last[i]]);  		   // the probability of an animal never being seen again after the last time it was captured

	 }

}

generated quantities {
 
          
}


