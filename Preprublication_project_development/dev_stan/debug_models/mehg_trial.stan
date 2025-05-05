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


  // dimensional and bookkeeping params (single vals)
	int<lower=1> n_pop;				    // Number of distinct populations (sampling areas)
	int<lower=1> n_pop_year;			    // Index for pop*year
	int<lower=1> n_ind;				    // Total number of individuals caught (ever, over all years and all populations)
	int<lower=1> ind_per_period_bd;	                    // n_ind * year (the unique periods of time in which each individual's bd is estimated)
	int<lower=1> ind_occ;			   	    // n_ind * all sampling periods (all events in which each individual could potentially have been captured)
	int<lower=1> ind_occ_min1;		 	    // n_ind * all sampling periods except the last 
	int<lower=1> n_days;				    // Number of sampling occasions
	int<lower=1> n_spec;				    // Total number of species
	int<lower=1> n_sex;				    // Total number of entries given for Sex (should be 3 if data cleaning worked [F, M, U]
	
  // short vector indexes (length of n_ind)	
	int<lower=1> ind_occ_size[n_ind];		    // Number of sampling periods for all individuals
	int<lower=1> ind_occ_min1_size[n_ind];		    // Number of sampling periods -1 for all individuals
	int<lower=1> phi_first_index[n_ind];		    // The indexes of phi corresponding to the first entry for each individual
	int<lower=1> p_first_index[n_ind];	            // The indexes of p corresponding to the first entry for each individual
	int<lower=1> ind_in_pop[n_ind];		   	    // Which population each individual belongs to
	matrix[n_ind, n_sex] ind_sex;		  	    // Sex of each individual
	matrix[n_ind, n_spec] ind_spec;			    // The species of each individual

  // short vector indexes (length of n_pop)
	matrix[n_pop, n_spec] spec_pop;			    // Which species is found in each pop_spec

  // short vector indexes (length of n_days)
	int<lower=1> day_which_pop[n_days];		    // Which pop is associated with each unique sampling day (for day-level detection deviates)
	int<lower=1> spec_which_pop[n_days];		    // Which species is associated with each sampling day
		
  // long vector indices for observation model (p)
	int<lower=1> ind_occ_rep[ind_occ];		    // Index vector of all individuals (each individual repeated the number of sampling occasions)
	int<lower=1> p_month[ind_occ];		            // Vector designating shorter periods (here month) for observational model (all occasions)
	int<lower=0, upper=1> p_zeros[ind_occ];		    // Observation times for each individual in which we do not know if that individual is present
	int<lower=1> p_bd_index[ind_occ];		    // Which entries of latent bd correspond to each entry of p
	int<lower=1> p_day[ind_occ];			    // Individual day identifier to try and estimate detection by day
	int<lower=1> pop_p[ind_occ];			    // Population index for detection predictors

	matrix[ind_occ, n_spec] spec_p;		   	    // Species identity of each individual for phi
  
  // long vector indices for survival model (phi)
	int<lower=1> ind_occ_min1_rep[ind_occ_min1];	    // Index vector of all individuals (each individual repeated the number of sampling occasions -1)
	int<lower=0, upper=1> offseason[ind_occ_min1];	    // Vector indicating the last sampling periods of each season which gains offseason characteristics
	int<lower=1> phi_year[ind_occ_min1];		    // Same as p_year but without last sampling event
	int<lower=0, upper=1> phi_zeros[ind_occ_min1];	    // Observation times for each individual in advance of first detecting that individual
	int<lower=0, upper=1> phi_ones[ind_occ_min1];	    // Time periods where we force survival to be 1 (assuming a closed population)
	int<lower=1> phi_bd_index[ind_occ_min1];	    // Which entries of latent bd correspond to each entry of phi
	int<lower=1> pop_phi[ind_occ_min1];		    // Population index for mortality predictors

	matrix[ind_occ_min1, n_spec] spec_phi;		    // Species identity of each individual for phi
	matrix[ind_occ_min1, n_sex] sex_phi;		    // Sex of each individual for phi
	  
  // individual-level covariates
	int<lower=1> N_bd;				    // Number of defined values for bd
 	real X_bd[N_bd];			   	    // The bd values 
	int<lower=1> X_ind[N_bd];			    // Individual associated with each bd measure
	int<lower=1> bd_first_index[ind_per_period_bd];	    // First entry of latent bd associated with each individual 'by' period
	int<lower=1> bd_last_index[ind_per_period_bd];	    // Last entry of latent bd associated with each individual 'by' period
	int<lower=1> x_bd_index[N_bd];			    // entries of X (latent bd) that have a corresponding real measure to inform likelihood with

  // long vector indexes: bd stuff (bd)
	int<lower=1> ind_bd_rep[ind_per_period_bd];	    // Index of individual for individual bd estimates (as each individual gets one estimate per year)    
	int<lower=1> ind_in_pop_year[ind_per_period_bd];    // Index of pop*year for individual bd estimates
	int<lower=1> pop_bd[ind_per_period_bd];             // Index of population for individual bd estimates

	matrix[ind_per_period_bd, n_spec] spec_bd;	    // Index of species identity for individual bd estimates	    

  // covariates (length)
	int<lower=1> n_ind_len_have;			    // Number of individuals that we have length data	  
	int<lower=1> n_ind_len_mis;			    // Number of individuals with missing length data
	int<lower=1> ind_len_which_have[n_ind_len_have];    // Index of individuals that we have length data
	int<lower=1> ind_len_which_mis[n_ind_len_mis];      // Index of individuals with missing length data
	vector[n_ind_len_have] ind_len_have;		    // The actual length values that we have
 	int<lower=1> ind_len_spec_first_index[n_spec]; 	    // First size index associated with each unique species (for species-specific scaling of length values)
 	int<lower=1> ind_len_spec_size[n_spec];      	    // Number of individuals of each species with lengths (for species-specific scaling of length values)

	matrix[n_ind_len_have, n_spec] ind_len_spec_have;   // The species of all individuals that we have lengths for, in model matrix form
	matrix[n_ind_len_mis, n_spec] ind_len_spec_mis;	    // The species of all individuals that we don't have lengths for, in model matrix form
	matrix[n_ind_len_have, n_sex] ind_len_sex_have;	    // The sex of all individuals that we have lengths for, in model matrix form
	matrix[n_ind_len_mis, n_sex] ind_len_sex_mis;	    // The sex of all individuals that we don't have lengths for, in model matrix form

  // covariates (MeHg)
 	int<lower=0> n_ind_mehg;			    // Number of individuals with measured MeHg
 	vector<lower=0>[n_ind_mehg] ind_mehg;		    // Measured values of MeHg
 	int<lower=0> ind_mehg_pop[n_ind_mehg];	   	    // Populations associated with each measure of MeHg

	matrix[n_ind_mehg, n_spec] ind_mehg_spec;	    // Species associated with each measure of MeHg 	    

  // site-level covariates, converted continuous from cat
	real<lower=0, upper=1> pop_drawdown[n_pop];	    // population specific covariate for proportion drawdown    
  	
  // site-by-day covariates, converted continuous from cat
	vector<lower=0, upper=1>[n_days] p_drawdown;   	    // average value of drawdown among all of the sampled SubSites on a given day (for daily detection) 
  	vector<lower=0, upper=1>[n_days] p_veg;        	    // average value of vegetation among all of the sampled SubSites on a given day (for daily detection)

  // site-level covariates, taken continuous
	real pop_temp[n_pop_year];		 	    // population*year specific covariate for temperature
	
  // captures
	int<lower=1> N_y;				    // Number of defined values for captures
  	int<lower=0, upper=1> y[N_y];		            // The capture values 
  	int<lower=1> first[n_ind];         		    // Capture event in which each individual was first captured
  	int<lower=1> last[n_ind];         		    // Capture event in which each individual was last captured
	vector<lower=0>[n_days] n_capt_per_day;	   	    // Number of captures on each day 

  // indices of phi, p, and chi that are 0, 1, or estimated, and which entries inform the likelihood.
  // set up in R to avoid looping over the full length of phi and p here. See R code for details
	int<lower=1> n_phi_zero;  
	int<lower=1> n_phi_one;      
	int<lower=1> n_phi_in;      
	int<lower=1> n_phi_off;    
	int<lower=1> phi_zero_index[n_phi_zero];  
	int<lower=1> phi_one_index[n_phi_one];
	int<lower=1> phi_in_index[n_phi_in];
	int<lower=1> phi_off_index[n_phi_off]; 

	int<lower=1> n_p_zero;
	int<lower=1> n_p_est;
	int<lower=1> p_zero_index[n_p_zero];
	int<lower=1> p_est_index[n_p_est];

	int<lower=1> n_phi_ll;
	int<lower=1> n_p_ll;
	int<lower=1> which_phi_ll[n_phi_ll];
	int<lower=1> which_p_ll[n_p_ll];
	int<lower=1> which_chi_ll[n_ind];

}

parameters {
// ------------------------------ parameters ------------------------------

// -----
// imputed covariates: MeHg
// -----

  // fixed
	real<lower=0> inverse_phi_mehg;		         // variance parameter for gamma regression
	vector[n_spec] beta_mehg_spec;			 // species-specific MeHg means 
	real beta_mehg_drawdown;			 // effect of drawdown on MeHg

  // random: variance
	real<lower=0> mehg_pop_sigma;			 // variation in mean MeHg by population

  // random: deviates
	real mehg_pop_eps[n_pop];			 // variation among pops

}


transformed parameters {
// ------------------------------ transformed parameters ------------------------------

  // MeHg
  	vector[n_ind_mehg] mu_mehg;	 		 // the expected values for the gamma regression
  	vector[n_ind_mehg] rate_mehg;	 	         // rate parameter for the gamma distribution
	real mehg_pop[n_pop];				 // MeHg contamination deviate by population
	vector[n_pop] mehg_pop_est;			 // mean MeHg contamination by population
	vector[n_pop] mehg_pop_est_scaled;		 // mean MeHg contamination by population scaled

// -----
// Estimated population-level mean MeHg 
// -----

  // calculate the mean at the population level; species, drawdown, and pop deviates on MeHg
	for (z in 1:n_pop) {
	  mehg_pop[z]     = mehg_pop_sigma * mehg_pop_eps[z];								
	  mehg_pop_est[z] = exp(spec_pop[z, ] * beta_mehg_spec + beta_mehg_drawdown * pop_drawdown[z] + mehg_pop[z]); 
	} 

  // scaled mehg population means
	mehg_pop_est_scaled = (mehg_pop_est - mean(mehg_pop_est))/sd(mehg_pop_est);     

  // mean estimate generated from each individuals measured bd; species, drawdown, and pop deviates on MeHg
	for (i in 1:n_ind_mehg) {
	  mu_mehg[i]  = exp(ind_mehg_spec[i, ] * beta_mehg_spec + beta_mehg_drawdown * pop_drawdown[ind_mehg_pop[i]] + mehg_pop[ind_mehg_pop[i]]);	
	}

  // transformation to get one of the parameters for the gamma distribution
	rate_mehg   = rep_vector(inverse_phi_mehg, n_ind_mehg) ./ mu_mehg;		
}

model {
// ------------------------------ model ------------------------------

// -----
// Priors
// -----


  // Imputed Covariates Priors: mehg

  // fixed
	beta_mehg_spec     ~ normal(0, 1);		// Narrow to constrain crazy estimates given little information content in some pops
	beta_mehg_drawdown ~ normal(0, 1);

  // variance
	inverse_phi_mehg  ~ inv_gamma(8, 15);	
	mehg_pop_sigma    ~ inv_gamma(8, 15);

 // deviates
	mehg_pop_eps ~ normal(0, 1);


// -----
// MeHg regression imputation
// -----

	ind_mehg ~ gamma(inverse_phi_mehg, rate_mehg);


}


