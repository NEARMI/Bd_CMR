####################################################################################################
## Processing of indices for the stan model to reduce looping for increasing computational speeds ##
####################################################################################################

## which individuals contribute to phi and p directly (individuals that were captured 2 or more times)
which_ind_ll      <- which(capture_range$first != capture_range$final)

## number of individuals that contribute to phi and p directly (used in the stan model to give length of which_ind_ll)
ind_ll_n <- which_ind_ll %>% length()

## which entries of the long-form phi and p are informed by these individuals? (these are the phi's and p's
 ## that will be used in the 1 ~ bernoulli() statements in the model block)
which_phi_ll <- numeric(0)
which_p_ll   <- numeric(0)

for (i in which_ind_ll) {
  ## extracting the phi's and p's for individual i (all of the times we know that individual i lives)
  these_phis   <- (phi_first_index[i] + (capture_range$first[i] + 1):capture_range$final[i] - 2)
  these_ps     <- (p_first_index[i] - 1 + ((capture_range$first[i] + 1):capture_range$final[i]))
  which_phi_ll <- c(which_phi_ll, these_phis)
  which_p_ll   <- c(which_p_ll, these_ps)
}

## do the same for chi, which requires a different lengthed loop because all individuals contribute one entry to chi
 ## (on the last day they were captured)
which_chi_ll    <- numeric(0)

for (i in 1:n_ind) {
  these_chis   <- p_first_index[i] - 1 + capture_range$final[i]
  which_chi_ll <- c(which_chi_ll, these_chis)
}


## Indices for which of all of the phi and p entries are constrained to be 0, 1, or estimated. Saves lots of time
 ## to calculate these indices here so that the stan model doesn't have to loop over some absurdly large number (phi and p are the
  ## longest of all of the entries as it is 1 per day per individual, so in the full model in the 400,000k range)

## phi's that must be 0 (phis for each individual before that individual was captured)
phi_zero_index  <- which(capt_history.phi$phi_zeros == 1) 

## phi's that must be 1 (phis for each individual inbetween captures)
phi_one_index   <- which(capt_history.phi$phi_ones == 1 & capt_history.phi$phi_zeros == 0) 

## phi's that are estimated ``within-season'' 
phi_in_index    <- which(capt_history.phi$offseason == 0 & capt_history.phi$phi_ones == 0 & capt_history.phi$phi_zeros == 0)

## phi's that are estimated ``between-season'' 
phi_off_index   <- which(capt_history.phi$offseason == 1 & capt_history.phi$phi_zeros == 0)  

## ps that must be 0 (ps for each individual before that individual was captured)
p_zero_index    <- which(capt_history.p$p_zeros == 0)

## all other ps
p_est_index     <- which(capt_history.p$p_zeros != 0) 
