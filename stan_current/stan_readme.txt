Current working stan models

-- Current primary working models --

SINGLE POPULATIONS
------------------

CMR_single_population_ind_rand_no_mehg
	-- Single population fit with length imputation and daily random effects for detection
CMR_single_population_ind_rand_mehg
	-- The above but also with MeHg imputation 



MULTIPLE POPULATIONS
--------------------

CMR_multiple_populations_full
	-- Model fitting to multiple species and populations
	-- As complicated as is probably sensible
	-- Uses categorical levels for a number of predictors

CMR_multiple_populations_reduced
	-- Reduces the number of predictors from ^^
	-- Converts the categorical factors to continuous (though a better strategy may be to use ordered factors)