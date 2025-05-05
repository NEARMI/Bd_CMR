functions {
// ------------------------------ functions ------------------------------

// -----
// chi_t probability estimator
// -----

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
// ------------------------------ data ------------------------------

// -----
// Notes
// -----
// Given long comments, this model is best read in full screen on an external monitor


  // dimensional and bookkeeping params (single vals)
	int<lower=1> n_pop;				    // Number of distinct populations (sampling areas)
	int<lower=1> n_pop_year;			    // Index for pop*year
	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years and all populations)
	int<lower=1> ind_per_period_bd;	                    // n_ind * year (the unique periods of time in which each individual's bd is estimated)
	int<lower=0> ind_occ;			   	    // n_ind * all sampling periods (all events in which each individual could potentially have been captured)
	int<lower=0> ind_occ_min1;		 	    // n_ind * all sampling periods except the last 
	int<lower=0> n_days;				    // Number of sampling occasions

	int<lower=0> n_spec;				    // Total number of species
	
  // short vector indexes (length of n_ind)	
	int<lower=1> ind_occ_size[n_ind];		    // Number of sampling periods for all individuals
	int<lower=1> ind_occ_min1_size[n_ind];		    // Number of sampling periods -1 for all individuals
	int<lower=1> phi_first_index[n_ind];		    // The indexes of phi corresponding to the first entry for each individual
	int<lower=1> p_first_index[n_ind];	            // The indexes of p corresponding to the first entry for each individual

	int<lower=0> ind_in_pop[n_ind];		   	    // Which population each individual belongs to

  // short vector indexes (length of n_days)
	int<lower=1> day_which_pop[n_days];
	int<lower=1> spec_which_pop[n_days];
		
  // long vector indices for observation model (p)
	int<lower=0> ind_occ_rep[ind_occ];		    // Index vector of all individuals (each individual repeated the number of sampling occasions)
	int<lower=0> p_month[ind_occ];		            // Vector designating shorter periods (here month) for observational model (all occasions)
	int<lower=0> p_zeros[ind_occ];			    // Observation times for each individual in which we do not know if that individual is present
	int<lower=0> p_bd_index[ind_occ];		    // Which entries of latent bd correspond to each entry of p
	int<lower=0> p_day[ind_occ];			    // Individual day identifier to try and estimate detection by day

	int<lower=0> pop_p[ind_occ];			    // Population index for detection predictors
	int<lower=0> spec_p[ind_occ];			    // Species identity of each individual
  
  // long vector indices for survival model (phi)
	int<lower=0> ind_occ_min1_rep[ind_occ_min1];	    // Index vector of all individuals (each individual repeated the number of sampling occasions -1)
	int<lower=0, upper=1> offseason[ind_occ_min1];	    // Vector indicating the last sampling periods of each season which gains offseason characteristics
	int<lower=0> phi_year[ind_occ_min1];		    // Same as p_year but without last sampling event
	int<lower=0> phi_zeros[ind_occ_min1];		    // Observation times for each individual in advance of first detecting that individual
	int<lower=0> phi_ones[ind_occ_min1];		    // Time periods where we force survival to be 1 (assuming a closed population)
	int<lower=0> phi_bd_index[ind_occ_min1];	    // Which entries of latent bd correspond to each entry of phi
	
	int<lower=0> pop_phi[ind_occ_min1];		    // Population index for mortality predictors
	int<lower=0> spec_phi[ind_occ_min1];		    // Species identity of each individual
	  
  // individual-level covariates
	int<lower=1> N_bd;				    // Number of defined values for bd
 	real X_bd[N_bd];			   	    // The bd values 
	int<lower=1> X_ind[N_bd];			    // Individual associated with each bd measure
	int<lower=0> bd_first_index[ind_per_period_bd];	    // First entry of latent bd associated with each individual 'by' period
	int<lower=0> bd_last_index[ind_per_period_bd];	    // Last entry of latent bd associated with each individual 'by' period
	int<lower=0> x_bd_index[N_bd];			    // entries of X (latent bd) that have a corresponding real measure to inform likelihood with

  // long vector indexes: bd stuff (bd)
	int<lower=0> ind_bd_rep[ind_per_period_bd];	    // Index of individual for individual bd estimates (as each individual gets one estimate per year)    
	int<lower=0> ind_in_pop_year[ind_per_period_bd];    // Index of pop*year for individual bd estimates

	int<lower=0> spec_for_bd[ind_per_period_bd];	    // Index of species identity for individual bd estimates
	int<lower=0> pop_for_bd[ind_per_period_bd];         // Index of population for individual bd estimates

  // covariates (length)
	int<lower=0> n_ind_len_have;			    // Number of individuals that we have length data	  
	int<lower=0> n_ind_len_mis;			    // Number of individuals with missing length data
	int<lower=0> ind_len_which_have[n_ind_len_have];    // Index of individuals that we have length data
	int<lower=0> ind_len_which_mis[n_ind_len_mis];      // Index of individuals with missing length data
	vector[n_ind_len_have] ind_len_have;		    // The actual length values that we have

  // site-level covariates, categorical
	int<lower=0> pop_drawdown[n_pop];		    // population specific covariate for proportion drawdown

  // site-level covariates, continuous
	real pop_temp[n_pop_year];		 	    // population*year specific covariate for temperature
	
  // captures
	int<lower=1> N_y;				    // Number of defined values for captures
  	int<lower=0, upper=1> y[N_y];		            // The capture values 

  	int<lower=0> first[n_ind];         		    // Capture event in which each individual was first captured
  	int<lower=0> last[n_ind];         		    // Capture event in which each individual was last captured

	vector<lower=0>[n_days] n_capt_per_day;	   	    // Number of captures on each day 

}


transformed data {

} 


parameters {
// ------------------------------ parameters ------------------------------

// -----
// bd submodel
// -----

  // fixed
	real beta_bd[n_spec];				 // species-level average bd level
	real beta_bd_temp;				 // population-level temperature effect on bd levels
	
  // random
	vector<lower=0>[n_pop] bd_ind_sigma;		 // variation in Bd among individuals (normal random effect variance) (nested in pops)		 
	real bd_ind_eps[n_ind];                          // the conditions modes of the random effect (each individual's intercept)

	real<lower=0> bd_pop_sigma;			 // variation in Bd among pop*year (normal random effect variance)
	real bd_pop_eps[n_pop_year];			 // the conditions modes of the random effect (each populations intercept)

	real<lower=0> bd_obs;    			 // observation noise for observed Bd compared to underlying state	


// -----
// survival
// -----

  // fixed
	vector[2] beta_offseason;  			 // survival as a function of bd stress
	real beta_inseason[n_spec];			 // in season survival intercept

  // random
	real<lower=0> offseason_pop_sigma;		 // variation in offseason survival by population (intercept)
	real offseason_pop_eps[n_pop];

	real<lower=0> offseason_pop_bd_sigma;		 // variation in offseason survival by population (slope over bd)
	real offseason_pop_bd_eps[n_pop];

	real<lower=0> inseason_pop_sigma;		 // variation in inseason survival by population (intercept)
	real inseason_pop_eps[n_pop];

// -----
// detection
// -----

  // fixed
	real beta_p[n_spec];				 // species-level average detection

  // random
	real<lower=0> p_pop_sigma;			 // variation in detection by population
	real p_pop_eps[n_pop];
	
	real<lower=0> p_day_delta_sigma[n_pop];	         // variation in detection by day (nested within each population)
	real p_day_delta_eps[n_days];

// -----
// imputed covariates ** need to make this variable by species or population
// -----

  // length
	vector[n_ind_len_mis] ind_len_mis;		 // the imputed values of length
	real<lower=0> ind_len_alpha; 			 // estimated gamma parameter of length distribution
	real<lower=0> ind_len_beta;			 // estimated gamma parameter of length distribution

}


transformed parameters {
// ------------------------------ transformed parameters ------------------------------

  // Individual lengths ** need to make this variable by species or population
	vector[n_ind] ind_len;				 // all individual lengths (combining data and imputed values)
	vector[n_ind] ind_len_scaled;			 // all individual lengths scaled
	real ind_len_mean;				 // mean of ind_len
	real ind_len_sd;				 // sd of ind_len

  // bd
	real bd_ind[n_ind];				 // individual random effect deviates
	real X[ind_per_period_bd];		         // each individual's estimated bd per year 
	
  // population:year random effect deviates
	real bd_pop_year[n_pop_year];

  // Survival
	real<lower=0,upper=1> phi[ind_occ_min1];         // survival from t to t+1, each individual repeated the number of times its population was measured
	real inseason_pop[n_pop];
	real offseason_pop[n_pop];
	real offseason_pop_bd[n_pop];

  // Detection 
	real<lower=0,upper=1> p[ind_occ];                // detection at time t
	real p_pop[n_pop];   				 // population-level detection deviates

	real p_day_dev[n_days];				 // day (nested in population) level detection deviates

	vector<lower=0,upper=1>[n_days] p_per_day;	 // detection per day at the level of the population

  // Chi estimator 
	real<lower=0,upper=1> chi[ind_occ];              // probability an individual will never be seen again

// -----
// Imputed NA Data values
// -----

  // Individual Length ** need to make this variable by species or population
	ind_len[ind_len_which_have] = ind_len_have;	 // filling in the complete vector of ind_lengths with the data
	ind_len[ind_len_which_mis]  = ind_len_mis;       // filling in the complete vector of ind_lengths with the imputed values

	ind_len_mean = mean(ind_len);			
	ind_len_sd   = sd(ind_len);

	ind_len_scaled = (ind_len - ind_len_mean)/ind_len_sd;


// -----
// bd submodel, contained to estimating within-season bd
// -----

	for (pp in 1:n_pop_year) {
	  bd_pop_year[pp] = bd_pop_sigma * bd_pop_eps[pp];
	} 

	for (i in 1:n_ind) {
	    
  // linear predictor for intercept for bd-response. Overall intercept + pop-specific intercept + individual random effect deviate
  	  bd_ind[i] = bd_ind_sigma[ind_in_pop[i]] * bd_ind_eps[i];  

	}

	for (t in 1:ind_per_period_bd) {

  // latent bd model before obs error (species effect + pop temp effect + ind deviate + pop*year deviate)
	  X[t] = (beta_bd[spec_for_bd[t]] + beta_bd_temp * pop_temp[pop_for_bd[t]] + bd_ind[ind_bd_rep[t]] + bd_pop_year[ind_in_pop_year[t]]);      

        }


// -----
// Survival probability over the whole period
// -----

	for (pp in 1:n_pop) {
	 inseason_pop[pp]     = inseason_pop_sigma  * inseason_pop_eps[pp];
	 offseason_pop[pp]    = offseason_pop_sigma * offseason_pop_eps[pp];
	 offseason_pop_bd[pp] = offseason_pop_bd_sigma * offseason_pop_bd_eps[pp];
	} 

	for (t in 1:ind_occ_min1) {

	 if (phi_zeros[t] == 1) {			// phi_zeros is 1 before an individual is caught for the first time
           phi[t] = 0;					// must be non-na values in stan, but the likelihood is only informed from first capture onward
	 } else {

	  if (offseason[t] == 0) {			// in season survival process

	   if (phi_ones[t] == 1) { 
	     phi[t] = 1;				// closed population assumption where survival is set to 1
           } else {

  // within-season survival given by a simple population-level intercept ** could increase this complexity eventually 
             phi[t] = inv_logit(beta_inseason[spec_phi[t]] + inseason_pop[pop_phi[t]]);

	   }

	  } else {

  // offseason survival given by a population-level intercept, a bd-effect, and a size effect ** MeHg and other site-level covariates to be added
	     phi[t] = inv_logit(
offseason_pop[pop_phi[t]] + 
(beta_offseason[1] + offseason_pop_bd[pop_phi[t]]) * X[phi_bd_index[t]] + 
beta_offseason[2] * ind_len_scaled[ind_occ_min1_rep[t]]
);

	  }

	 }  

	}


// -----
// Detection probability over the whole period
// -----
	
  // population level random effect
	for (pp in 1:n_pop) {
	  p_pop[pp]    = p_pop_sigma * p_pop_eps[pp];
	}

	for (i in 1:n_days) {
  	  p_day_dev[i]  = p_day_delta_sigma[day_which_pop[i]] * p_day_delta_eps[i];  
	}

	for (t in 1:ind_occ) {
	 if (p_zeros[t] == 0) {       // just a placeholder because stan can't have NA -- these are periods before an animal was captured that doesn't impact the likelihood
	   p[t] = 0;
	 } else {
           p[t] = inv_logit(p_pop[pop_p[t]] + beta_p[spec_p[t]] * ind_len_scaled[ind_occ_rep[t]] + p_day_dev[p_day[t]]);
	 }
	}

	for (t in 1:n_days) {
	  p_per_day[t] = inv_logit(p_pop[day_which_pop[t]] + beta_p[spec_which_pop[t]] + p_day_dev[t]);
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
// ------------------------------ model ------------------------------

// -----
// Priors
// -----

  // Bd Model Priors

	bd_pop_sigma   ~ inv_gamma(8, 15);
	bd_obs         ~ inv_gamma(10, 4);

	beta_bd        ~ normal(0, 3);
	beta_bd_temp   ~ normal(0, 3);

	for (i in 1:n_pop) {
	  bd_ind_sigma[i] ~ inv_gamma(8, 15);
	}

	for (i in 1:n_pop_year) {
	  bd_pop_eps[i]   ~ normal(0, 3);
	}

  // Survival Priors

	beta_offseason[1]   ~ normal(0, 1.15);
	beta_offseason[2]   ~ normal(0, 1.15);

	for (i in 1:n_spec) {
	  beta_inseason[i]  ~ normal(0, 1.15);
	}

	offseason_pop_sigma    ~ inv_gamma(8, 15);
	offseason_pop_bd_sigma ~ inv_gamma(8, 15);
	inseason_pop_sigma     ~ inv_gamma(8, 15);

	for (i in 1:n_ind) {
	  bd_ind_eps[i]   ~ normal(0, 3);
	}

	for (i in 1:n_pop) {
	  p_pop_eps[i]            ~ normal(0, 1.15);
	  inseason_pop_eps[i]     ~ normal(0, 1.45);
	  offseason_pop_eps[i]    ~ normal(0, 1.15);
	  offseason_pop_bd_eps[i] ~ normal(0, 1.15);
	}

  // Detection Priors

	p_pop_sigma    ~ inv_gamma(8, 15);

	for (i in 1:n_spec) {
	  beta_p[i]    ~ normal(0, 1.15);
	}

	for (i in 1:n_pop) {
	  p_day_delta_sigma[i] ~ inv_gamma(8, 15);
	}

	for (i in 1:n_days) {
	  p_day_delta_eps[i]   ~ normal(0, 1.15);
	}

         
  // Imputed Covariates Priors: length

	ind_len_alpha ~ inv_gamma(10, 4);	
	ind_len_beta  ~ inv_gamma(10, 4);


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

	  if (first[i] != last[i]) {

	   for (t in (first[i] + 1):last[i]) {			
	    1 ~ bernoulli(phi[phi_first_index[i] - 1 + t - 1]);    		   // Survival _to_ t (from phi[t - 1]) is 1 because we know the individual lived in that period 
	    y[p_first_index[i] - 1 + t] ~ bernoulli(p[p_first_index[i] - 1 + t]);  // Capture given detection
	   }

          }

	    1 ~ bernoulli(chi[p_first_index[i] - 1 + last[i]]);  		   // the probability of an animal never being seen again after the last time it was captured

	 }

}

generated quantities {
// ------------------------------ generated quantities ------------------------------
 
  vector<lower=0>[n_days] pop_size;
  pop_size = n_capt_per_day ./ p_per_day;
          
}


