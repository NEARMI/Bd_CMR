Current working stan models

-- Current primary working models --

SINGLE POPULATIONS
------------------

CMR_single_population.stan
	-- Older models (see dev_stan) now with things like sex in survival

CMR_single_population_gl.stan
	-- ^^ with expanded gamma regression for imputed lengths. 
	-- Leaving this and the above in stan dev for now because the length model hasn't been adequately debugged yet

CMR_single_population_mehg_gl.stan
	-- Single population with individual-level MeHg imputation and expanded Gamma length imputation
	-- Still need to debug the length imputation aspect



MULTIPLE POPULATIONS
--------------------

CMR_multiple_populations.stan
	-- Older models (see dev_stan) now with things like sex in survival

CMR_multiple_populations_gl.stan
	-- ^^ with expanded gamma regression for imputed lengths. 
	-- Leaving this and the above in stan dev for now because the length model hasn't been adequately debugged yet

