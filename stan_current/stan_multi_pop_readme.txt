MULTIPLE POPULATIONS STAN MODELS 
--------------------------------

-- Current primary working models --
  (for older models see dev_stan)
 


[1] CMR_multiple_populations.stan
	-- The top level model

[2] CMR_multiple_populations_mehg.stan
	-- [1] with individual-level MeHg imputation

CMR_multiple_populations_ssp.stan
	-- [1] for only the populations of one species 

CMR_multiple_populations_mehg_ssp.stan
	-- [2] for only the populations of one species (probably only fittable for RANA and maybe ANBO)

CMR_multiple_populations_ssp_mv.stan
	-- ^^ but with explicit Cholesky decomposition for the random effects to model the covariance among them. Tends to be a bit faster/more efficient which is also good. 
