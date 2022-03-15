Current working stan models

-- Current primary working models --

CMR_single_population_con_mi2_p_phi_adj
	-- One population
	-- continuous survival within season given number of days between sampling
	-- multiple imputation second attempt

CMR_single_population_con_mi2_p_phi_adj_p_per_day
	-- Trying a random effect for detection, allowing each day to have its own detection probability

CMR_single_population_con_mi2_p_phi_adj_p_per_day2
	-- Above but going back to periods with a closed population and periods with an open population

CMR_single_population_con_mi2_p_phi_adj_p_per_day2_ind_rand
	-- Above but also adding individual-level random effects for detection and within-season survival probabilities
	




[TO DO]
	-- Figure out what to do with Site, SubSite and sampling days









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
