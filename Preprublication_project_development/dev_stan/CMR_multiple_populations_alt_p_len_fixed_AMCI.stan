functions {
// ------------------------------ functions ------------------------------

// -----
// chi_t probability estimator
// -----

// Function to build the required chi_t probability estimator for the mark recapture model, which
 // is the probability of never recapturing an individual again after capturing them at time t
  // rewritten to run for an individual at a time given variable numbers of capture opportunities by individual (e.g. in different populations)
	
	real[] prob_uncaptured(int n_occ, vector p_sub, vector phi_sub) {

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

  // dimensional and bookkeeping data (non-vectors)
	int<lower=1> n_pop;				    // Number of distinct populations (sampling areas)
	int<lower=1> n_pop_year;			    // Index for pop*year
	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years and all populations)
	int<lower=1> ind_per_period_bd;	                    // n_ind * year (the unique periods of time in which each individual's bd is estimated)
	int<lower=1> ind_occ;			   	    // n_ind * all sampling periods (all events in which each individual could potentially have been captured)
	int<lower=1> ind_occ_min1;		 	    // n_ind * all sampling periods except the last 
	int<lower=1> n_days;				    // Number of sampling occasions
	int<lower=1> n_days_for_p;			    // Number of unique sampling day "types" informing the levels of the detection random effect
	int<lower=1> n_spec;				    // Total number of species
	int<lower=1> n_sex;				    // Number of unique entries for sex (M, F, U)
	int<lower=1> N_bd;				    // Number of defined values for bd
	int<lower=1> n_col_mm_int;			    // Number of unique intercepts for detection and survival models (i.e., number of columns in the model matrix)
	int<lower=1> n_u;				    // Number of random effects in between season survival portion
	
  // Index vectors with length ``n_ind'' (used in all model components)	
	int<lower=1> ind_occ_size[n_ind];		    // Number of sampling periods for all individuals
	int<lower=1> ind_occ_min1_size[n_ind];		    // Number of sampling periods -1 for all individuals
	int<lower=1> phi_first_index[n_ind];		    // The indexes of phi corresponding to the first entry for each individual
	int<lower=1> p_first_index[n_ind];	            // The indexes of p corresponding to the first entry for each individual
	int<lower=1> ind_in_pop[n_ind];		   	    // Which population each individual belongs to

  // Index vectors with length ``n_days'' or ``n_days_for_p''
	int<lower=1> day_which_pop[n_days_for_p];	    // Which pop is associated with each unique sampling day ''type'' (for day-level detection deviates)
	int<lower=1> p_rand_which_day[n_days];		    // Which random detection value is associated with each day
	int<lower=1> day_which_pop_rand[n_days];	    // Which pop is associated with each actual sampling day
		
  // Components for detection model (p)
	int<lower=1> n_p_zero;				    // Which entries of the longer vector p are set to zero
	int<lower=1> n_p_est;				    // Which entries of the longer vector p inform the likelihood

	int<lower=1> p_day[ind_occ];			    // Individual day identifier to try and estimate detection by day
	int<lower=1> pop_p[ind_occ];			    // Population index for detection predictors
	int<lower=1> ind_for_p[n_p_est];		    // Repeated individual for length for p predict
	matrix[n_p_est, n_col_mm_int] fe_mm_p_int;	    // Intercept component of the model matrix

  // Components for population size estimates
	int<lower=1> n_fe_mm_p_int_uni;				  // number of unique intercept combos
	matrix[n_fe_mm_p_int_uni, n_col_mm_int] fe_mm_p_int_uni;  // the unique model matrix entries (unique intercepts)
	matrix[n_days, 2] fe_mm_p_slope_uni;			  // entries of the slope covariates for detection for each day
	int<lower=1> spec_to_int[n_spec, n_sex];		  // which entries of the unique combos of the intercept are appropriate for the three sexes for each species
	int<lower=1> spec_pop_se[n_days];			  // the species present in each population, used to select the appropriate intercept for a given population
  
  // Components for survival model (phi) (Index vectors and model matrices)
	int<lower=1> n_phi_zero;  			    // Which entries of the longer vector phi are set to zero
	int<lower=1> n_phi_one;      			    // Which entries of the longer vector phi are set to one
	int<lower=1> n_phi_in;      			    // Which entries of the longer vector phi inform within season survival
	int<lower=1> n_phi_off;  			    // Which entries of the longer vector phi inform between season survival

	int<lower=1> ind_occ_min1_rep[ind_occ_min1];	    // Index vector of all individuals (each individual repeated the number of sampling occasions -1)
	int<lower=1> phi_bd_index[ind_occ_min1];	    // Which entries of latent bd correspond to each entry of phi
	int<lower=1> pop_phi[ind_occ_min1];		    // Population index for mortality predictors
	matrix[n_ind, n_spec] ind_spec;			    // The species of each individual	
	matrix[n_phi_off, n_col_mm_int] fe_mm_phi_int;      // Intercept component of the model matrix
	matrix[n_phi_off, n_spec] fe_mm_phi_slope;	    // Slope component of the model matrix
	  
  // Components for Bd model (Bd)
 	real X_bd[N_bd];			   	    // The bd values 
	int<lower=1> x_bd_index[N_bd];			    // Entries of X (latent bd) that have a corresponding real measure to inform likelihood with
	int<lower=1> ind_bd_rep[ind_per_period_bd];	    // Index of individual for individual bd estimates (as each individual gets one estimate per year)    
	int<lower=1> ind_in_pop_year[ind_per_period_bd];    // Index of pop*year for individual bd estimates
	int<lower=1> pop_bd[ind_per_period_bd];             // Index of population for individual bd estimates
	matrix[ind_per_period_bd, n_col_mm_int] spec_bd;    // Index of species identity for individual bd estimates

  // Components for length imputation (Dimensions, Index vectors, covariates, and model matrices)
	int<lower=1> n_ind_len_have;			    // Number of individuals that we have length data	  
	int<lower=1> n_ind_len_mis;			    // Number of individuals with missing length data
	int<lower=1> ind_len_which_have[n_ind_len_have];    // Index of individuals that we have length data
	int<lower=1> ind_len_which_mis[n_ind_len_mis];      // Index of individuals with missing length data
	vector[n_ind_len_have] ind_len_have;		    // The actual length values that we have
 	int<lower=1> ind_len_spec_first_index[n_spec]; 	    // First size index associated with each unique species (for species-specific scaling of length values)
 	int<lower=1> ind_len_spec_size[n_spec];      	    // Number of individuals of each species with lengths (for species-specific scaling of length values)
	matrix[n_ind, n_col_mm_int] ind_mm_len;		    // Intercept component of the model matrix

  // Site-level covariates
	real pop_drawdown[n_pop];	   		    // population specific covariate for proportion drawdown    
	real pop_temp[n_pop_year];		 	    // population*year specific covariate for temperature
	
  // Capture data
	int<lower=1> N_y;				    // Number of defined values for captures
  	int<lower=0, upper=1> y[N_y];		            // The capture values 
	matrix[n_days, n_sex] n_capt_per_day_sex;	    // Number of individuals of all sexes captured in a given day


  // Indices of phi, p, and chi that are 0, 1, or estimated, 
  // set up in R to avoid looping over the full length of phi and p here. See R code for details  
	int<lower=1> phi_zero_index[n_phi_zero];  
	int<lower=1> phi_one_index[n_phi_one];
	int<lower=1> phi_in_index[n_phi_in];
	int<lower=1> phi_off_index[n_phi_off]; 

	int<lower=1> p_zero_index[n_p_zero];
	int<lower=1> p_est_index[n_p_est];

  // Which entries of phi, p, and chi inform the likelihood.
  // set up in R to avoid looping over the full length of phi and p here. See R code for details
	int<lower=1> n_phi_ll;
	int<lower=1> n_p_ll;
	int<lower=1> which_phi_ll[n_phi_ll];
	int<lower=1> which_p_ll[n_p_ll];
	int<lower=1> which_chi_ll[n_ind];

}

parameters {
// ------------------------------ parameters ------------------------------

// -----
// bd submodel
// -----

  // fixed
	vector[n_col_mm_int] beta_bd_pop;		 // population-level average bd level
	real beta_bd_temp;				 // population-level temperature effect on bd levels
	real beta_bd_len;				 // individual-specific length effect on bd levels

  // error
	real<lower=0> bd_obs;    			 // observation noise for observed Bd compared to underlying state	
	
  // random: variance
	vector<lower=0>[n_pop] bd_ind_sigma;		 // variation in Bd among individuals (normal random effect variance) (nested in pops)		 
	real<lower=0> bd_pop_sigma;			 // variation in Bd among pop*year (normal random effect variance)

  // random: conditional modes (hereafter deviates)
	real bd_ind_eps[n_ind];         		 // individual specific bd load adjustment                
	real bd_pop_eps[n_pop_year];			 // pop-by-year average adjustment


// -----
// survival: offseason
// -----

  // fixed
	vector[n_col_mm_int] beta_offseason_int;	 // Intercept for between season survival
	vector[n_pop] beta_offseason_bd;		 // Bd effect on between season survival
	vector[n_pop] beta_offseason_len;		 // Length effect on between season survival 

// -----
// survival: inseason
// -----

  // fixed
	vector[n_pop] beta_inseason;			 // in season survival intercept

// -----
// detection
// -----

  // fixed
	vector[n_col_mm_int] beta_p_int;		 // species-level average detection
	vector[2] beta_p_slope; 			 // daily detection probably as a function of drawdown and vegetation
	real beta_p_len;	

  // random: variance
	real<lower=0> p_day_delta_sigma[n_pop];	         // variation in detection by day (nested within each population)

  // random: deviates
	real p_day_delta_eps[n_days_for_p];		 // day to day variation


// -----
// imputed covariates: length
// -----

	vector[n_col_mm_int] beta_len;
	vector[n_ind_len_mis] ind_len_mis;		 // the imputed values of len
  	real<lower=0> sd_len[n_pop]; 

}


transformed parameters {
// ------------------------------ transformed parameters ------------------------------

  // Individual lengths 
  	vector[n_ind_len_have] mu_len_have; 		  	 
  	vector[n_ind_len_mis] mu_len_mis; 
	vector[n_ind] ind_len;				 // all individual len (combining data and imputed values)
	vector[n_ind] ind_len_scaled;			 // all individual len scaled

  // bd
	real bd_ind[n_ind];				 // individual random effect deviates
	vector[ind_per_period_bd] X;		         // each individual's estimated bd per year
	vector[ind_per_period_bd] X_scaled;
	
  // population:year random effect deviates
	real bd_pop_year[n_pop_year];			 // pop-by-year variation in Bd level

  // Detection 
	vector[n_days_for_p] p_day_dev;			 // day (nested in population) level detection deviates

  // Long-form containers for estimates from t to t+1
	vector<lower=0,upper=1>[ind_occ_min1] phi;       // survival from t to t+1, each individual repeated the number of times its population was measured
	vector<lower=0,upper=1>[ind_occ] p;              // detection at time t
	real<lower=0,upper=1> chi[ind_occ];              // Chi estimator -- probability an individual will never be seen again if they are alive at time t


// -----
// Imputed NA Length values
// -----

  // linear predictor for len regression on measured lengths
	mu_len_have   = ind_mm_len[ind_len_which_have, ] * beta_len;  

  // linear predictor for len regression for imputing unknown lengths
  	mu_len_mis    = ind_mm_len[ind_len_which_mis, ] * beta_len;

  // filling in the complete vector of ind_len with the data and missing values
	ind_len[ind_len_which_have] = ind_len_have;
	ind_len[ind_len_which_mis]  = ind_len_mis;

  // Scaling the predicted lengths within-species (loop over species)
	for (ns in 1:n_spec) {

  // Jump through a hoop to select out all of the length values for a given species 
	  vector[ind_len_spec_size[ns]] temp_ind_len = segment(ind_len, ind_len_spec_first_index[ns], ind_len_spec_size[ns]);

  // Scale the lengths of species ns and stick them in the complete long-form container
	  ind_len_scaled[ind_len_spec_first_index[ns]:(ind_len_spec_first_index[ns] + ind_len_spec_size[ns] - 1)] = (temp_ind_len - mean(temp_ind_len))/sd(temp_ind_len);

	}


// -----
// bd submodel, contained to estimating within-season bd
// -----

  // pop-spec deviates
	for (pp in 1:n_pop_year) {
	  bd_pop_year[pp] = bd_pop_sigma * bd_pop_eps[pp];
	} 

  // linear predictor for intercept for bd-response. Overall intercept + pop-specific intercept + individual random effect deviate
	for (i in 1:n_ind) {
  	  bd_ind[i] = bd_ind_sigma[ind_in_pop[i]] * bd_ind_eps[i];  
	}

  // latent bd model before obs error (species effect + pop temp effect + ind deviate + pop*year deviate)
	for (t in 1:ind_per_period_bd) {
	  X[t] = spec_bd[t, ] * beta_bd_pop + bd_ind[ind_bd_rep[t]] + bd_pop_year[ind_in_pop_year[t]] + beta_bd_temp * pop_temp[pop_bd[t]] + beta_bd_len * ind_len_scaled[ind_bd_rep[t]];      
        }

	X_scaled = (X - mean(X)) / sd(X);


// -----
// Survival probability over the whole period
// -----

	phi[phi_zero_index] = rep_vector(0, n_phi_zero);
	phi[phi_one_index]  = rep_vector(1, n_phi_one);

	phi[phi_in_index]   = inv_logit(ind_spec[ind_occ_min1_rep[phi_in_index], ] * beta_inseason);

	phi[phi_off_index]  = inv_logit(
fe_mm_phi_int    * beta_offseason_int  +
(fe_mm_phi_slope * beta_offseason_bd) .* X_scaled[phi_bd_index[phi_off_index]] +
(fe_mm_phi_slope * beta_offseason_len) .* ind_len_scaled[ind_occ_min1_rep[phi_off_index]]
);


// -----
// Detection probability over the whole period
// -----

  // day-level deviates
	for (i in 1:n_days_for_p) {
  	  p_day_dev[i]  = p_day_delta_sigma[day_which_pop[i]] * p_day_delta_eps[i];  
	}

	p[p_zero_index] = rep_vector(0, n_p_zero);

	p[p_est_index]  = inv_logit(
fe_mm_p_int * beta_p_int      + 
p_day_dev[p_day[p_est_index]] +
fe_mm_p_slope * beta_p_slope  +
beta_p_len * ind_len_scaled[ind_for_p]
);


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
	beta_bd_pop   ~ normal(0, 3);
	beta_bd_temp  ~ normal(0, 3);
	beta_bd_len   ~ normal(0, 3);

  // variances 
	bd_pop_sigma  ~ inv_gamma(8, 15);
	bd_ind_sigma  ~ inv_gamma(8, 15);
	bd_obs        ~ inv_gamma(10, 4);

  // deviates
	bd_pop_eps ~ normal(0, 3);
	bd_ind_eps ~ normal(0, 3);

  // Survival Priors

  // fixed
	beta_offseason_int  ~ normal(0, 0.65);
	beta_offseason_bd   ~ normal(0, 0.65);
	beta_offseason_len  ~ normal(0, 0.65);
	beta_inseason       ~ normal(0, 1.75);

  // Detection Priors

  // fixed
	beta_p_int   ~ normal(0, 0.85);
	beta_p_slope ~ normal(0, 0.85);
	beta_p_len   ~ normal(0, 0.85);

  // variances and deviates
	p_day_delta_sigma ~ inv_gamma(8, 15);
	p_day_delta_eps   ~ normal(0, 0.85);


  // Imputed Covariates Priors: length

	sd_len    ~ inv_gamma(8, 15);	
	beta_len  ~ normal(0, 3);

// -----
// Imputed NA Lengths
// -----

	ind_len_have ~ normal(mu_len_have, sd_len[ind_in_pop[ind_len_which_have]]);
	ind_len_mis  ~ normal(mu_len_mis, sd_len[ind_in_pop[ind_len_which_mis]]);

// -----
// Bd Process and Data Model
// -----

	X_bd ~ normal(X[x_bd_index], bd_obs);
    

// -----
// Capture model
// -----

	1 ~ bernoulli(phi[which_phi_ll]);
	y[which_p_ll] ~ bernoulli(p[which_p_ll]);
	1 ~ bernoulli(chi[which_chi_ll]);

}

generated quantities {
// ------------------------------ generated quantities ------------------------------
          
  vector<lower=0>[n_days] pop_size;
  vector[n_fe_mm_p_int_uni] beta_p_each_int;
  matrix[n_days, n_sex] p_per_day;

  for (i in 1:n_fe_mm_p_int_uni) {
   beta_p_each_int[i] = fe_mm_p_int_uni[i, ] * beta_p_int;
  } 

  for (i in 1:n_days) {

   p_per_day[i] = to_row_vector(
    inv_logit(
     beta_p_each_int[spec_to_int[spec_pop_se[i], ]]     + 
     rep_vector(p_day_dev[p_rand_which_day[i]], n_sex)  + 
     rep_vector(fe_mm_p_slope_uni[i, ] * beta_p_slope, n_sex)
    )
   );

   pop_size[i]  = sum(n_capt_per_day_sex[i, ] ./ p_per_day[i]);

  }
 
}

