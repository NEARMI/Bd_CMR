R scripts
---------

Individual_CMR.R              -- Only some individuals are swabbed for bd, but those that are are always swabbed when caught
Individual_CMR_low_samples.R  -- All individuals are occasionally swabbed when caught


Stan model names broken down into sections

e.g. : CMR_ind_all_no.stan
        1   2   3   4

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


Models
------

Models build on each other moving down this list

CMR_ind_all_no    -- All individuals have the same Bd
CMR_ind_some_no   -- Individuals start in different places with their bd load but have the same slope
CMR_ind_some_no_2 -- ... and have different slopes

CMR_ind_some_bd-p-phi -- Individuals also vary in their response to detection and mortality (and not just to their variable bd load)
CMR_ind_some_bd-p-phi_all -- ^^ with estimating p and phi before first capture
CMR_ind_some_bd-p-phi_restrict -- ^^ with setting p and phi before first capture to 0

CMR_ind_pat_bd-p-phi     -- Modifying 'CMR_ind_some_bd-p-phi_restrict' for the actual structure of how
                              bd is sampled (from time to time for some individuals)
CMR_ind_pat_bd_empirical -- ^^ but with random survival and detectability not as a function of bd to capture more of the individual variation
CMR_ind_pat_bd_empirical_red -- ^^ but with a completely different bd process model and no extra random effects to remove parameters
		^^ Important note about this model in the script
CMR_ind_pat_bd_empirical_red2 -- Extremely reduced form of ^^ to try and get it to fit to real data
CMR_ind_pat_bd_empirical_red3 -- ^^ but trying to adjust the n_occasions to restrict the likelihood to only get updated on days where sampling actually was taking place 
CMR_ind_pat_bd_empirical_red4 -- ^^ trying to debug. Problems with X as a transformed parameter and not parameter... 
	^^^^^ Will continue to just add a number to work on model development...

CMR_ind_pat_bd_empirical_01 -- Step back and try and fit just a binary of infected or not since load has so many issues


