SINGLE POPULATIONS STAN MODELS 
------------------------------

-- Current primary working models --
  (for older models see dev_stan)
 



[1] CMR_single_population.stan
	-- Base single population model for most populations

CMR_single_population_nli.stan
	-- [1] without the need for length imputation because all animals have lengths measured

CMR_single_population_nl.stan
	-- [1] without length at all because no animals had lengths measured

CMR_single_population_no_gl.stan
	-- [1] without sex in length imputation because one sex was missing all lengths

CMR_single_population_no_in.stan
	-- [1] without any in-season for populations with none of the defined "in-season" time gaps


[2] CMR_single_population_mehg.stan
	-- [1] but also with individual-level MeHg imputation

CMR_single_population_mehg_nli.stan
	-- [2] with lo length imputation because all animals had lengths measured



[3] CMR_single_population_con.stan
	-- [1] but with seasonally varying bd loads at the individual level

