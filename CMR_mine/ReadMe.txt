Stan model names broken down into sections

e.g. : CMR_ind_all_no.stan
        1   2   3   4

1) designates a CMR 

2) type of covariates
     ind = individual-based

3) how the covariates are collected
     all  = perfectly, no gaps -- no model for the covariates, just brought in as data
     some = observed for some captured individuals -- thus the model contains a model for the latent covariate

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

CMR_ind_some_bd_empirical -- Modifying 'CMR_ind_some_bd-p-phi_restrict' for the actual structure of how
                              bd is sampled (from time to time for some individuals)
 
