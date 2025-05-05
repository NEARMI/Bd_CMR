This folder contains development and scratch stan models moving towards the working models in the main folder. Keeping these for my own notes and posterity sake 
-----------


Stan model names broken down into sections

e.g. : CMR_ind_all_no_XXXX.stan
        1   2   3   4   5

1) designates a CMR 

2) type of covariates
     ind = individual-based

3) how the covariates are collected
     all  = perfectly, no gaps -- no model for the covariates, just brought in as data
     some = always observed for a subset of those individuals captured  -- thus the model contains a model for the latent covariate
     pat  = sampling happens for all individuals captured but incompletely

4) random effects at the level of individual?
     no  = all individuals behave the same to their covariates (even if the covariates vary) 
     bd-p-phi = individuals vary in their Bd profiles and individuals vary in their capture probability and mortality probability

5) many possible extensions that designate:
     matrix or database form
     temperature as an added covariate
     adding recruitment between years
     playing with scaling for individuals not seen but which could be present in a previous year
     etc. (most models are commented reasonably well)


More extensive model notes: NOTE: I gave up keeping details after some time (see bullet 5 above)

Models build on each other moving down this list

CMR_ind_all_no    -- All individuals have the same Bd
CMR_ind_some_no   -- Individuals start in different places with their bd load but have the same slope
CMR_ind_some_no_2 -- ... and have different slopes

CMR_ind_some_bd-p-phi -- Individuals also vary in their response to detection and mortality (and not just to their variable bd load)
CMR_ind_some_bd-p-phi_all -- ^^ with estimating p and phi before first capture
CMR_ind_some_bd-p-phi_restrict -- ^^ with setting p and phi before first capture to 0

CMR_ind_pat_bd-p-phi     -- Modifying 'CMR_ind_some_bd-p-phi_restrict' for the actual structure of how
                              bd is sampled (from time to time for some individuals)

CMR_ind_pat_bd-p-phiXX   -- All of these are various attempts at debugging the next steps of uneven gaps in time between sampling events and other covariates 


CMR_ind_pat_bd_empirical -- ^^ but with random survival and detectability not as a function of bd to capture more of the individual variation
CMR_ind_pat_bd_empirical_red -- ^^ but with a completely different bd process model and no extra random effects to remove parameters
		^^ Important note about this model in the script
CMR_ind_pat_bd_empirical_red2 -- Extremely reduced form of ^^ to try and get it to fit to real data
CMR_ind_pat_bd_empirical_red3 -- ^^ but trying to adjust the n_occasions to restrict the likelihood to only get updated on days where sampling actually was taking place 
CMR_ind_pat_bd_empirical_red4 -- ^^ trying to debug. Problems with X as a transformed parameter and not parameter... 
	^^^^^ Will continue to just add a number to work on model development...

CMR_ind_pat_bd_empirical_01 -- Step back and try and fit just a binary of infected or not since load has so many issues

CMR_ind_pat_bd-p-phi_multi 
CMR_ind_pat_bd-p-phi_multi2 -- Converting the single season model into a multi-season model. Work in progress but seems to run


CMR_ind_pat_bd-p-phi_no_timegap_covariate
	- Random effect in bd intercept only
	- No individual variation in response to bd and no other random effects
	- Sampling occasions collapsed, meaning that the amount of time elapsed between sampling events is ignored
	- Survival from one sampling event to the next, even if there is a gap in time, is based on bd at the time of the first sampling event
CMR_ind_pat_bd-p-phi_no_average_bd
	- Random effect in bd intercept only
	- No individual variation in response to bd and no other random effects
	- Sampling occasions collapsed, but a covariate added to try and control for how much time elapsed between sampling events
	- Survival from one sampling event to the next, even if there is a gap in time, is based on bd at the time of the first sampling event
CMR_ind_pat_bd-p-phi_average_bd
	- Random effect in bd intercept only
	- No individual variation in response to bd and no other random effects
	- Sampling occasions collapsed, but a covariate added to try and control for how much time elapsed between sampling events
	- Survival from one sampling event to the next, is based on average estimated bd over the period between the sampling events


