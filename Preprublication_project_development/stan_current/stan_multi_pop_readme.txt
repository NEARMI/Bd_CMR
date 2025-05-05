MULTIPLE POPULATIONS STAN MODELS 
--------------------------------

-- Current primary working models --
  (all other models moved to /dev_stan)
 
[1] CMR_multiple_populations.stan
	-- The top level model

[2] CMR_multiple_populations_mehg.stan
	-- [1] with individual-level MeHg imputation

CMR_multiple_populations_ssp.stan
	-- [1] for only the populations of one species 

CMR_multiple_populations_mehg_ssp.stan
	-- [2] for only the populations of one species (probably only fittable for RANA and maybe ANBO)

CMR_multiple_populations_ssp_mv.stan, CMR_multiple_populations_ssp_fixed.stan moved to /dev_stan


CMR_multiple_populations_alt_p.stan
	-- [1] With a simplified model for p to allow for faster fits when the newt populations are included.
		-- Fits a random effect for detection by pop*primary period instead of per day

CMR_multiple_populations_alt_p_len.stan
	-- ^^ Above model but also with length in the detection part of the model



**** NOTE: 
	-- May need to create alt_p for the MeHg model, but may not since newts are not present and runs fast enough