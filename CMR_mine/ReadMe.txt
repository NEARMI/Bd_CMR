Stan model names broken down into sections

e.g. : CMR_ind_all_no.stan
        1   2   3   4

1) designates a CMR 

2) type of covariates
     ind = individual-based

3) how the covariates are collected
     all = perfectly, no gaps -- no model for the covariates, just brought in as data

4) random effects at the level of individual?
     no  = all individuals behave the same to their covariates (even if the covariates vary) 
