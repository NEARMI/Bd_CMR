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
	int<lower=0> n_sex;				    // Total number of entries given for Sex (should be 3 if data cleaning worked [F, M, U]
	
  // short vector indexes (length of n_ind)	
	int<lower=1> ind_occ_size[n_ind];		    // Number of sampling periods for all individuals
	int<lower=1> ind_occ_min1_size[n_ind];		    // Number of sampling periods -1 for all individuals
	int<lower=1> phi_first_index[n_ind];		    // The indexes of phi corresponding to the first entry for each individual
	int<lower=1> p_first_index[n_ind];	            // The indexes of p corresponding to the first entry for each individual
	int<lower=0> ind_in_pop[n_ind];		   	    // Which population each individual belongs to
	int<lower=0> ind_sex[n_ind];			    // The sex of each individual
	int<lower=0> ind_spec[n_ind];			    // The species of each individual

  // short vector indexes (length of n_pop)
	int<lower=1> spec_pop[n_pop];			    // Which species is found in each pop_spec

  // short vector indexes (length of n_days)
	int<lower=1> day_which_pop[n_days];		    // Which pop is associated with each unique sampling day (for day-level detection deviates)
	int<lower=1> spec_which_pop[n_days];		    // Which species is associated with each sampling day
		
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
 	int<lower=0> ind_len_spec_first_index[n_spec]; 	    // First size index associated with each unique species (for species-specific scaling of length values)
 	int<lower=0> ind_len_spec_size[n_spec];      	    // Number of individuals of each species with lengths (for species-specific scaling of length values)

  // covariates (MeHg)
 	int<lower=0> n_ind_mehg;			    // Number of individuals with measured MeHg
 	vector<lower=0>[n_ind_mehg] ind_mehg;		    // Measured values of MeHg
 	int<lower=0> ind_mehg_pop[n_ind_mehg];	   	    // Populations associated with each measure of MeHg
 	int<lower=0> ind_mehg_spec[n_ind_mehg];	   	    // Species associated with each measure of MeHg

  // site-level covariates, forced categorical
  	int<lower=0> pop_sub[n_pop];      		    // population specific covariate for substrate
  	int<lower=0> pop_region[n_pop];   		    // population specific covariate for region of the country
  	int<lower=0> pop_hydro[n_pop];    		    // population specific covariate for hydro period

  // site-level covariates, categorical but potentially continuous
	real<lower=0> pop_drawdown[n_pop];		    // population specific covariate for proportion drawdown    
  	
  // site-by-day level covariates, categorical but potentially continuous
	real<lower=0> p_drawdown[n_days];   		    // average value of drawdown among all of the sampled SubSites on a given day (for daily detection) 
  	real<lower=0> p_veg[n_days];        		    // average value of vegetation among all of the sampled SubSites on a given day (for daily detection)

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
	real beta_bd_spec[n_spec];			 // species-level average bd level
	real beta_bd_temp;				 // population-level temperature effect on bd levels
	real beta_bd_len;				 // individual-specific length effect on bd levels

  // error
	real<lower=0> bd_obs;    			 // observation noise for observed Bd compared to underlying state	
	
  // random: variance
	vector<lower=0>[n_pop] bd_ind_sigma;		 // variation in Bd among individuals (normal random effect variance) (nested in pops)		 
	real<lower=0> bd_pop_sigma;			 // variation in Bd among pop*year (normal random effect variance)

  // random: deviates (conditional modes of the random effects)
	real bd_ind_eps[n_ind];                         
	real bd_pop_eps[n_pop_year];			


// -----
// survival: offseason
// -----

  // fixed
	real beta_offseason_int[n_spec];		 // Intercept for between season survival
	real beta_offseason_bd[n_spec];			 // Bd effect on between season survival
	real beta_offseason_len[n_spec];		 // Length effect on between season survival 
	real beta_offseason_mehg[n_spec];		 // MeHg effect on between season survival
	real beta_offseason_sex[3];			 // Sex effect on survival 

  // random: variance
	real<lower=0> offseason_pop_sigma;		 // variation in offseason survival by population (intercept)
	real<lower=0> offseason_pop_bd_sigma;		 // variation in offseason survival by population (slope over bd)
	real<lower=0> offseason_pop_len_sigma;		 // variation in offseason survival by population (slope over animal length)

  // random: deviates
	real offseason_pop_eps[n_pop];
	real offseason_pop_bd_eps[n_pop];
	real offseason_pop_len_eps[n_pop];


// -----
// survival: inseason
// -----

  // fixed
	real beta_inseason[n_spec];			 // in season survival intercept

  // random: variance
	real<lower=0> inseason_pop_sigma;		 // variation in inseason survival by population (intercept)

  // random: deviates
	real inseason_pop_eps[n_pop];


// -----
// detection
// -----

  // fixed
	real beta_p_spec[n_spec];			 // species-level average detection
	real beta_p_drawdown;				 // daily detection probably as a function of drawdown amount
	real beta_p_veg;			         // daily detection probably as a function of vegetation amount

  // random: variance
	real<lower=0> p_pop_sigma;			 // variation in detection by population
	real<lower=0> p_day_delta_sigma[n_pop];	         // variation in detection by day (nested within each population)

  // random: deviates
	real p_day_delta_eps[n_days];
	real p_pop_eps[n_pop];


// -----
// imputed covariates: length
// -----

	real<lower=0> inverse_phi_len;		         // variance parameter for gamma regression
	vector[n_sex] beta_len_sex;			 // regression coefficient for len as a function of sex
	vector[n_spec] beta_len_spec;			 // regression coefficient for len as a function of species
	

	real<lower=0> len_pop_sigma;			 // variation in len by population
	real len_pop_eps[n_pop];

	vector[n_ind_len_mis] ind_len_mis;		 // the imputed values of len


// -----
// imputed covariates: MeHg
// -----

  // fixed
	real<lower=0> inverse_phi_mehg;		         // variance parameter for gamma regression
	real beta_mehg_spec[n_spec];			 // species-specific MeHg means 
	real beta_mehg_drawdown;			 // effect of drawdown on MeHg

  // random: variance
	real<lower=0> mehg_pop_sigma;			 // variation in mean MeHg by population

  // random: deviates
	real mehg_pop_eps[n_pop];

}


transformed parameters {
// ------------------------------ transformed parameters ------------------------------

  // Individual lengths 
  	vector[n_ind_len_have] mu_len_have; 		 // the expected values for the gamma regression
  	vector[n_ind_len_have] rate_len_have; 	 	 // rate parameter for the gamma distribution

  	vector[n_ind_len_mis] mu_len_mis; 		 // the expected values (linear predictor) for the missing len values
  	vector[n_ind_len_mis] rate_len_mis; 		 // rate parameter for the gamma distribution for the missing len values

	vector[n_pop] pop_len_dev;			 // population specific length deviate

	vector[n_ind] ind_len;				 // all individual len (combining data and imputed values)
	vector[n_ind] ind_len_scaled;			 // all individual len scaled
	vector[n_spec] ind_len_mean;			 // mean of ind_len
	vector[n_spec] ind_len_sd;			 // sd of ind_len

  // mehg
  	vector[n_ind_mehg] mu_mehg;	 		 // the expected values for the gamma regression
  	vector[n_ind_mehg] rate_mehg;	 	         // rate parameter for the gamma distribution
	real mehg_pop[n_pop];				 // MeHg contamination deviate by population
	vector[n_pop] mehg_pop_est;			 // mean MeHg contamination by population
	vector[n_pop] mehg_pop_est_scaled;		 // mean MeHg contamination by population scaled

  // bd
	real bd_ind[n_ind];				 // individual random effect deviates
	real X[ind_per_period_bd];		         // each individual's estimated bd per year 
	
  // population:year random effect deviates
	real bd_pop_year[n_pop_year];			 // pop-by-year variation in Bd level

  // Survival
	real inseason_pop[n_pop];			 // population specific survival within season (intercept)
	real offseason_pop[n_pop];			 // population specific survival between seasons (intercept)
	real offseason_pop_bd[n_pop];			 // population specific survival between seasons (slope over bd)
	real offseason_pop_len[n_pop];			 // population specific survival between seasons (slope over animal length)

  // Detection 
	real p_pop[n_pop];   				 // population-level detection deviates
	real p_day_dev[n_days];				 // day (nested in population) level detection deviates
	vector<lower=0,upper=1>[n_days] p_per_day;	 // detection per day at the level of the population

  // Long-form containers for estimates from t to t+1
	real<lower=0,upper=1> phi[ind_occ_min1];         // survival from t to t+1 for each individual (in long-form)
	real<lower=0,upper=1> p[ind_occ];                // detection at time t
	real<lower=0,upper=1> chi[ind_occ];              // Chi estimator -- probability an individual will never be seen again if they are alive at time t


// -----
// Imputed NA Length values
// -----

  // Population-specific length deviates
	for (pl in 1:n_pop) {
	  pop_len_dev[pl] = len_pop_sigma * len_pop_eps[pl];
	}

  // linear predictor for len regression on measured lengths
  	mu_len_have   = exp(
beta_len_sex[ind_sex[ind_len_which_have]] + 
beta_len_spec[ind_spec[ind_len_which_have]] +
pop_len_dev[ind_in_pop[ind_len_which_have]]
);  
	
  	rate_len_have = rep_vector(inverse_phi_len, n_ind_len_have) ./ mu_len_have;

  // linear predictor for len regression for imputing unknown lengths
  	mu_len_mis    = exp(
beta_len_sex[ind_sex[ind_len_which_mis]] + 
beta_len_spec[ind_spec[ind_len_which_mis]] +
pop_len_dev[ind_in_pop[ind_len_which_mis]]
);   	

  	rate_len_mis = rep_vector(inverse_phi_len, n_ind_len_mis) ./ mu_len_mis;

	ind_len[ind_len_which_have] = ind_len_have;	 // filling in the complete vector of ind_mehg with the data
	ind_len[ind_len_which_mis]  = ind_len_mis;       // filling in the complete vector of ind_mehg with the imputed values

  // Scaling the predicted lengths within-species
	for (ns in 1:n_spec) {

  // Jump through a hoop to select out all of the length values for a given species 
	  vector[ind_len_spec_size[ns]] temp_ind_len = segment(ind_len, ind_len_spec_first_index[ns], ind_len_spec_size[ns]);

	  ind_len_mean[ns] = mean(temp_ind_len);	 // take the mean of the lengths of all individuals of species ns
	  ind_len_sd[ns]   = sd(temp_ind_len);		 // take the sd of the lengths of all individuals of species ns

  // Scale the lengths of species ns and stick them in the complete long-form container
	  ind_len_scaled[ind_len_spec_first_index[ns]:(ind_len_spec_first_index[ns] + ind_len_spec_size[ns] - 1)] = (temp_ind_len - ind_len_mean[ns])/ind_len_sd[ns];

	}


// -----
// Estimated population-level mean MeHg 
// -----

	for (z in 1:n_pop) {
	  mehg_pop[z]     = mehg_pop_sigma * mehg_pop_eps[z];				// individual population deviates

  // calculate the mean at the population level	
	  mehg_pop_est[z] = exp(
beta_mehg_spec[spec_pop[z]]          + 							// species effect on MeHg			
beta_mehg_drawdown * pop_drawdown[z] + 							// drawdown effect on MeHg
mehg_pop[z]										// population deviate on MeHg
);
	} 

	mehg_pop_est_scaled = (mehg_pop_est - mean(mehg_pop_est))/sd(mehg_pop_est);     // scaled mehg population means

	for (i in 1:n_ind_mehg) {

  // mean estimate generated from each individuals measured bd
	  mu_mehg[i]  = exp(
beta_mehg_spec[ind_mehg_spec[i]]                   + 					// species effect on MeHg
beta_mehg_drawdown * pop_drawdown[ind_mehg_pop[i]] + 					// drawdown effect on MeHg
mehg_pop[ind_mehg_pop[i]]								// population deviate on MeHg
);

	}

	rate_mehg   = rep_vector(inverse_phi_mehg, n_ind_mehg) ./ mu_mehg;		// transformation to get one of the parameters for the gamma distribution


// -----
// bd submodel, contained to estimating within-season bd
// -----

	for (pp in 1:n_pop_year) {
	  bd_pop_year[pp] = bd_pop_sigma * bd_pop_eps[pp];				// pop-spec deviates
	} 

	for (i in 1:n_ind) {
	    
  // linear predictor for intercept for bd-response. Overall intercept + pop-specific intercept + individual random effect deviate
  	  bd_ind[i] = bd_ind_sigma[ind_in_pop[i]] * bd_ind_eps[i];  

	}

	for (t in 1:ind_per_period_bd) {

  // latent bd model before obs error (species effect + pop temp effect + ind deviate + pop*year deviate)
	  X[t] = beta_bd_spec[spec_for_bd[t]] + 
bd_ind[ind_bd_rep[t]] + 
bd_pop_year[ind_in_pop_year[t]] +
beta_bd_temp * pop_temp[pop_for_bd[t]] + 
beta_bd_len * ind_len_scaled[ind_bd_rep[t]];      

        }


// -----
// Survival probability over the whole period
// -----

  // pop-spec deviates
	for (pp in 1:n_pop) {
	 inseason_pop[pp]       = inseason_pop_sigma       * inseason_pop_eps[pp];
	 offseason_pop[pp]      = offseason_pop_sigma      * offseason_pop_eps[pp];
	 offseason_pop_bd[pp]   = offseason_pop_bd_sigma   * offseason_pop_bd_eps[pp];
	 offseason_pop_len[pp]  = offseason_pop_len_sigma  * offseason_pop_len_eps[pp];
	} 

	for (t in 1:ind_occ_min1) {

	 if (phi_zeros[t] == 1) {			// phi_zeros is 1 before an individual is caught for the first time
           phi[t] = 0;					// must be non-na values in stan, but the likelihood is only informed from first capture onward
	 } else {

	  if (offseason[t] == 0) {			// in season survival process

	   if (phi_ones[t] == 1) { 
	     phi[t] = 1;				// closed population assumption where survival is set to 1
           } else {

  // within-season survival given by a simple population-level intercept; could potentially consider increasing the complexity of this component eventually 
             phi[t] = inv_logit(beta_inseason[spec_phi[t]] + inseason_pop[pop_phi[t]]);

	   }

	  } else {

  // linear predictor for offseason survival -- given by a population-level intercept, a bd-effect, and a size effect 
	     phi[t] = inv_logit(
beta_offseason_int[spec_phi[t]]  + offseason_pop[pop_phi[t]]       + 
(beta_offseason_bd[spec_phi[t]]  + offseason_pop_bd[pop_phi[t]])   * X[phi_bd_index[t]] + 
(beta_offseason_len[spec_phi[t]] + offseason_pop_len[pop_phi[t]])  * ind_len_scaled[ind_occ_min1_rep[t]] +
beta_offseason_mehg[spec_phi[t]] * mehg_pop_est_scaled[pop_phi[t]] +
beta_offseason_sex[ind_sex[ind_occ_min1_rep[t]]]
);

	  }

	 }  

	}


// -----
// Detection probability over the whole period
// -----
	
  // population-level deviates
	for (pp in 1:n_pop) {
	  p_pop[pp]    = p_pop_sigma * p_pop_eps[pp];
	}

  // day-level deviates
	for (i in 1:n_days) {
  	  p_day_dev[i]  = p_day_delta_sigma[day_which_pop[i]] * p_day_delta_eps[i];  
	}

	for (t in 1:ind_occ) {

  // just a placeholder because stan can't have NA -- these are periods before an animal was captured that doesn't impact the likelihood
	 if (p_zeros[t] == 0) {       
	   p[t] = 0;
	 } else {
  // linear predictor for detection -- given by pop and day level deviates, drawdown, veg, and a length effect broken up by species
           p[t] = inv_logit(
p_pop[pop_p[t]]                        + 
p_day_dev[p_day[t]]                    + 
beta_p_drawdown * p_drawdown[p_day[t]] + 
beta_p_veg * p_veg[p_day[t]]           + 
beta_p_spec[spec_p[t]] * ind_len_scaled[ind_occ_rep[t]]
);
	 }
	}

	for (t in 1:n_days) {
  // calculating a daily estimate of detection to obtain a population size estimate
	  p_per_day[t] = inv_logit(
p_pop[day_which_pop[t]]         + 
p_day_dev[t]                    + 
beta_p_drawdown * p_drawdown[t] + 
beta_p_veg * p_veg[t]
);
	}	


// -----
// Probability of never detecting an individual again after time t
// -----

  // For each individual calculate the probability it won't be captured again
	for (i in 1:n_ind) {
  // Again, for chi have to use segment to jump through the hoop of extracting the correct 
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

  // fixed
	beta_bd_spec   ~ normal(0, 3);
	beta_bd_temp   ~ normal(0, 3);
	beta_bd_len    ~ normal(0, 3);

  // variances 
	bd_pop_sigma      ~ inv_gamma(8, 15);
	for (i in 1:n_pop) {
	  bd_ind_sigma[i] ~ inv_gamma(8, 15);			// Each population has its own amount of individual variation in bd
	}

	bd_obs            ~ inv_gamma(10, 4);

  // deviates
	for (i in 1:n_pop_year) {
	  bd_pop_eps[i]   ~ normal(0, 3);
	}

	for (i in 1:n_ind) {
	  bd_ind_eps[i] ~ normal(0, 3);
	}


  // Survival Priors

  // fixed
	beta_offseason_int  ~ normal(0, 0.65);
	beta_offseason_bd   ~ normal(0, 0.65);
	beta_offseason_len  ~ normal(0, 0.65);
	beta_offseason_mehg ~ normal(0, 0.65);
	beta_offseason_sex  ~ normal(0, 0.65);

	for (i in 1:n_spec) {
	  beta_inseason[i]  ~ normal(0, 1.45);
	}

  // variances
	offseason_pop_sigma     ~ inv_gamma(8, 15);
	offseason_pop_bd_sigma  ~ inv_gamma(8, 15);
	offseason_pop_len_sigma ~ inv_gamma(8, 15);
	inseason_pop_sigma      ~ inv_gamma(8, 15);

  // deviates
	for (i in 1:n_pop) {
	  inseason_pop_eps[i]      ~ normal(0, 1.45);
	  offseason_pop_eps[i]     ~ normal(0, 0.85);
	  offseason_pop_bd_eps[i]  ~ normal(0, 0.85);
	  offseason_pop_len_eps[i] ~ normal(0, 0.85);
	}


  // Detection Priors

  // fixed
	beta_p_spec     ~ normal(0, 0.65);
	beta_p_drawdown ~ normal(0, 0.65);
	beta_p_veg      ~ normal(0, 0.65);

  // variances and deviates
	p_pop_sigma    ~ inv_gamma(8, 15);

	for (i in 1:n_pop) {
	  p_day_delta_sigma[i] ~ inv_gamma(8, 15);
	  p_pop_eps[i]         ~ normal(0, 1.15);
	}

	for (i in 1:n_days) {
	  p_day_delta_eps[i]   ~ normal(0, 0.65);
	}


  // Imputed Covariates Priors: length

	inverse_phi_len  ~ inv_gamma(8, 15);	
	beta_len_sex     ~ normal(0, 3);
	beta_len_spec    ~ normal(0, 3);

	len_pop_sigma    ~ inv_gamma(8, 15);

	for (i in 1:n_pop) {
	  len_pop_eps[i] ~ normal(0, 3);
	}


  // Imputed Covariates Priors: mehg

  // fixed
	beta_mehg_spec     ~ normal(0, 1);		// Narrow to constrain crazy estimates given little information content
	beta_mehg_drawdown ~ normal(0, 1);

  // variance
	inverse_phi_mehg  ~ inv_gamma(8, 15);	
	mehg_pop_sigma    ~ inv_gamma(8, 15);

 // deviates
	for (i in 1:n_pop) {
	  mehg_pop_eps[i] ~ normal(0, 1);		// Narrow to constrain populations with no measures
	}


// -----
// Imputed NA Lengths
// -----

	ind_len_have ~ gamma(inverse_phi_len, rate_len_have);
	ind_len_mis  ~ gamma(inverse_phi_len, rate_len_mis);


// -----
// MeHg regression imputation
// -----

	ind_mehg ~ gamma(inverse_phi_mehg, rate_mehg);


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

  // for animals only captured once there are no periods where we know survival is one, so skip these calculations for these individuals
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
 
  // estimated population sizes based on daily detection and number of animals captured
	vector<lower=0>[n_days] pop_size;
 	pop_size = n_capt_per_day ./ p_per_day;
          
}


