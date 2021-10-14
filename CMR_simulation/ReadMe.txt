Data Simulations and CMR model development

R scripts
---------

Individual_CMR_expanding.R
 -- long form R script to simulate data and fit stan model

Redoing_structure.R
 -- temporary script in progress as of OCT 14: working to convert stan model from matrix form to database form 

CMR_simulations.R
 -- script in progress as of OCT 14: working to convert Individual_CMR_expanding.R into functions for more convenient simulation of N populations with different structure


Stan Models
-----------

CMR_ind_pat_bd-p-phi_multi_recruit_free_temp.stan
 -- End step for the matrix model for a multi-year single population

CMR_ind_pat_bd-p-phi_multi_recruit_free_temp_db_simple_mp.stan
 -- Current in progress stan model that takes in long data to support multiple populations with different sampling schemes


Folders
-------

dev_stan
 -- dump of development and scratch stan models working towards the models in this main folder

dev_R
 -- dump of development and scratch R scripts working towards the models in this main folder

JAGS_attempts
 -- Briefly thought of moving to JAGS and then abandoned that idea

