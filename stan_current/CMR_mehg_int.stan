
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
	int<lower=1> n_spec;				    // Total number of species
	int<lower=1> N_bd;				    // Number of defined values for bd
	int<lower=1> n_col_mm_int;			    // Number of unique intercepts for detection and survival models (i.e., number of columns in the model matrix)
	
  // Index vectors with length ``n_ind'' (used in all model components)	
	int<lower=1> ind_in_pop[n_ind];		   	    // Which population each individual belongs to
	  
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
	matrix[n_ind, n_col_mm_int] ind_mm_len;	  	    // Intercept component of the model matrix
	matrix[n_ind, n_spec] ind_spec;			    // The species of each individual

  // Components for MeHg model (Dimensions, Index vectors, covariates, and model matrices)
 	int<lower=0> n_ind_mehg_have;			    // Number of individuals with measured MeHg
	int<lower=0> n_ind_mehg_mis;			    // Number of individuals with missing MeHg data
	int<lower=0> ind_mehg_which_have[n_ind_mehg_have];  // Index of individuals that we have mehg data
	int<lower=0> ind_mehg_which_mis[n_ind_mehg_mis];    // Index of individuals with missing mehg data
	vector[n_ind_mehg_have] ind_mehg_have;		    // The actual mehg values that we have
	
 	int<lower=1> ind_mehg_spec_first_index[n_spec];     // First size index associated with each unique species (for species-specific scaling of length values)
 	int<lower=1> ind_mehg_spec_size[n_spec];      	    // Number of individuals of each species with lengths (for species-specific scaling of length values)

  // Site-level covariates
	vector[n_pop] pop_drawdown;	   		    // population specific covariate for proportion drawdown    
	real pop_temp[n_pop_year];		 	    // population*year specific covariate for temperature
	
}

parameters {
// ------------------------------ parameters ------------------------------

// -----
// bd submodel
// -----

  // fixed
	vector[n_col_mm_int] beta_bd_spec;		 // species-level average bd level
	real beta_bd_temp;				 // population-level temperature effect on bd levels
	real beta_bd_len;				 // individual-specific length effect on bd levels
	real beta_bd_mehg;				 // impact of individual-level MeHg on individual Bd loads

  // error
	real<lower=0> bd_obs;    			 // observation noise for observed Bd compared to underlying state	
	
  // random: variance
	vector<lower=0>[n_pop] bd_ind_sigma;		 // variation in Bd among individuals (normal random effect variance) (nested in pops)		 
	real<lower=0> bd_pop_sigma;			 // variation in Bd among pop*year (normal random effect variance)

  // random: conditional modes (hereafter deviates)
	real bd_ind_eps[n_ind];         		 // individual specific bd load adjustment                
	real bd_pop_eps[n_pop_year];			 // pop-by-year average adjustment


// -----
// imputed covariates: length
// -----

	vector[n_col_mm_int] beta_len;
	vector[n_ind_len_mis] ind_len_mis;		 // the imputed values of len
  	real<lower=0> sd_len[n_pop]; 


// -----
// imputed covariates: MeHg
// -----

  // fixed
	real<lower=0> inverse_phi_mehg;		         // variance parameter for gamma regression
	vector[n_col_mm_int] beta_mehg;			 // species-by-sex-specific MeHg means 
	real beta_mehg_drawdown;			 // effect of drawdown on MeHg
	real beta_mehg_len;				 // individual length effect on MeHg
	vector[n_ind_mehg_mis] ind_mehg_mis;		 // the imputed values of mehg

  // random: variance
	real<lower=0> mehg_pop_sigma;			 // variation in mean MeHg by population

  // random: deviates
	real mehg_pop_eps[n_pop];			 // variation among pops

}


transformed parameters {
// ------------------------------ transformed parameters ------------------------------

  // Individual lengths 
  	vector[n_ind_len_have] mu_len_have; 		  	 
  	vector[n_ind_len_mis] mu_len_mis; 
	vector[n_ind] ind_len;				 // all individual len (combining data and imputed values)
	vector[n_ind] ind_len_scaled;			 // all individual len scaled

  // MeHg

  	vector[n_ind_mehg_have] mu_mehg_have; 		 // the expected values for the gamma regression
  	vector[n_ind_mehg_have] rate_mehg_have; 	 // rate parameter for the gamma distribution

  	vector[n_ind_mehg_mis] mu_mehg_mis; 		 // the expected values (linear predictor) for the missing mehg values
  	vector[n_ind_mehg_mis] rate_mehg_mis; 		 // rate parameter for the gamma distribution for the missing mehg values

	vector[n_ind] ind_mehg;				 // all individual mehg (combining data and imputed values)
	vector[n_ind] ind_mehg_scaled;			 // all individual mehg scaled

	vector[n_pop] mehg_pop;				 // pop random deviates

  // bd
	real bd_ind[n_ind];				 // individual random effect deviates
	vector[ind_per_period_bd] X;		         // each individual's estimated bd per year
	vector[ind_per_period_bd] X_scaled;
	
  // population:year random effect deviates
	real bd_pop_year[n_pop_year];			 // pop-by-year variation in Bd level


// -----
// Imputed NA Length values
// -----

  // linear predictor for len regression on measured lengths
	mu_len_have   = ind_mm_len[ind_len_which_have, ] * beta_len;  

  // linear predictor for len regression for imputing unknown lengths
  	mu_len_mis    = ind_mm_len[ind_len_which_mis, ] * beta_len;

  // filling in the complete vector of ind_mehg with the data and missing values
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
// Imputed NA MeHg values
// -----

  // calculate the mean at the population level; species, drawdown, and pop deviates on MeHg
	for (z in 1:n_pop) {
	  mehg_pop[z]     = mehg_pop_sigma * mehg_pop_eps[z];								
	} 

// linear predictor for mehg regression	
  	mu_mehg_have   = exp(
ind_mm_len[ind_mehg_which_have, ] * beta_mehg +
mehg_pop[ind_in_pop[ind_mehg_which_have]] +
beta_mehg_drawdown * pop_drawdown[ind_in_pop[ind_mehg_which_have]] + 
beta_mehg_len * ind_len_scaled[ind_mehg_which_have]
);  
  	rate_mehg_have = rep_vector(inverse_phi_mehg, n_ind_mehg_have) ./ mu_mehg_have;

// linear predictor for mehg regression
  	mu_mehg_mis   = exp(
ind_mm_len[ind_mehg_which_mis, ] * beta_mehg +
mehg_pop[ind_in_pop[ind_mehg_which_mis]] +
beta_mehg_drawdown * pop_drawdown[ind_in_pop[ind_mehg_which_mis]] + 
beta_mehg_len * ind_len_scaled[ind_mehg_which_mis]
);    	
  	rate_mehg_mis = rep_vector(inverse_phi_mehg, n_ind_mehg_mis) ./ mu_mehg_mis;


// filling in the complete vector of ind_mehg with the data and imputed values
	ind_mehg[ind_mehg_which_have] = ind_mehg_have;	    						    
	ind_mehg[ind_mehg_which_mis]  = ind_mehg_mis;


  // Scaling the predicted lengths within-species (loop over species)
	for (ns in 1:n_spec) {

  // Jump through a hoop to select out all of the length values for a given species 
	  vector[ind_mehg_spec_size[ns]] temp_ind_mehg = segment(ind_mehg, ind_mehg_spec_first_index[ns], ind_mehg_spec_size[ns]);

  // Scale the lengths of species ns and stick them in the complete long-form container
	  ind_mehg_scaled[ind_mehg_spec_first_index[ns]:(ind_mehg_spec_first_index[ns] + ind_mehg_spec_size[ns] - 1)] = (temp_ind_mehg - mean(temp_ind_mehg))/sd(temp_ind_mehg);

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
	  X[t] = spec_bd[t, ] * beta_bd_spec + bd_ind[ind_bd_rep[t]] + bd_pop_year[ind_in_pop_year[t]] + beta_bd_temp * pop_temp[pop_bd[t]] + beta_bd_len * ind_len_scaled[ind_bd_rep[t]] + beta_bd_mehg * ind_mehg_scaled[ind_bd_rep[t]];      
        }

	X_scaled = (X - mean(X)) / sd(X);

}

model {
// ------------------------------ model ------------------------------

// -----
// Priors
// -----

  // Bd Model Priors

  // fixed
	beta_bd_spec  ~ normal(0, 3);
	beta_bd_temp  ~ normal(0, 3);
	beta_bd_len   ~ normal(0, 3);
	beta_bd_mehg  ~ normal(0, 3);

  // variances 
	bd_pop_sigma  ~ inv_gamma(8, 15);
	bd_ind_sigma  ~ inv_gamma(8, 15);
	bd_obs        ~ inv_gamma(10, 4);

  // deviates
	bd_pop_eps ~ normal(0, 3);
	bd_ind_eps ~ normal(0, 3);


  // Imputed Covariates Priors: length
	sd_len    ~ inv_gamma(8, 15);	
	beta_len  ~ normal(0, 3);


  // Imputed Covariates Priors: mehg

  // fixed
	beta_mehg	   ~ normal(0, 2);		// Narrow to constrain crazy estimates given little information content in some pops
	beta_mehg_drawdown ~ normal(0, 1);
	beta_mehg_len      ~ normal(0, 1);

  // variance
	inverse_phi_mehg  ~ inv_gamma(8, 15);	
	mehg_pop_sigma    ~ inv_gamma(8, 15);

 // deviates
	mehg_pop_eps ~ normal(0, 1);


// -----
// Imputed NA Lengths
// -----

	ind_len_have ~ normal(mu_len_have, sd_len[ind_in_pop[ind_len_which_have]]);
	ind_len_mis  ~ normal(mu_len_mis, sd_len[ind_in_pop[ind_len_which_mis]]);


// -----
// MeHg regression imputation
// -----

	ind_mehg_have ~ gamma(inverse_phi_mehg, rate_mehg_have);
	ind_mehg_mis  ~ gamma(inverse_phi_mehg, rate_mehg_mis);


// -----
// Bd Process and Data Model
// -----

	X_bd ~ normal(X[x_bd_index], bd_obs);
    

}
