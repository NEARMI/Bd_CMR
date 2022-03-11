Current working stan models


-- Current primary working model --

CMR_single_population_con_mi2
	-- One population
	-- continuous survival within season given number of days between sampling
	-- multiple imputation second attempt

[TO DO]
	-- Check multiple imputation
	-- clean up prim and sec periods with # days between samples for both phi and p
	-- figure out what to do with gamma

{IN DEV working through TO DO list}

CMR_single_population_con_mi2_p_phi_adj
	-- Adjusting phi and p 


-- Current primary working model --


CMR_multiple_populations_red_ni_cov
	-- Many populations
	-- reduced complexity (no Bd process in detection or within-season survival)
	-- no intercept
	-- using site-level covariates

[TO DO]
	-- Add commenting from the single population model
	-- Estimate site-level mean mercury from individual-level mercury
	-- update phi, p, and gamma according to choices in single population model
	-- try species as fixed effect location as random
	-- try fitting with just ANBO

{IN DEV working through TO DO list}

None yet
