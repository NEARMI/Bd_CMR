CMR models with a latent disease process. Only files relevant to the models fit for the publication in Wildlife Letters.
This is a child of the parent Bd_CMR and branch Scientific_Reports.

Some notes about this code and Stan CMR model:
	-- Most of this code was originally written for a larger CMR analysis of 20 populations of 10 species of amphibian (published in Scientific Reports, 2025). The single population analyzed here was a part of this larger analysis. The code has been adjusted for fitting this single population with the addition of CORT. 
	-- Some code and modeling choices seem overly complicated for what would be needed to analyze a single population, such as the indexing in the Stan model. These choices were to accommodate the 20 populations of variable amounts of sampling and would indeed be odd if this model was written from scratch for this single population. However, code cleaning and model fit are identical to a simpler structure and will give the same answers, but are, unfortunately, a bit more opaque.

To run code, open "top_level_script.R" and work from there. All other scripts are sourced from within this script. The use of these scripts is described in "top_level_script.R" and in each individual sourced script. 


-- R Scripts -- (alphabetically, for the order that scripts are used open and work through "top_level_script.R")

top_level_script.R         -- top level script from which all others are called and analysis is run

data_covariates_JP_cort.R  -- load and clean covariate data
data_manip.R               -- data manipulation 
data_stan.R                -- set up data in structure needed for stan model
explore_model_fit.R        -- explore model fit
JP_data.R                  -- load data
packages_functions.R       -- needed packages and functions
pseudo_growth_model_cort.R -- deal with missing cortstan_fit_single_JP_cort.R  -- send data to stan model and fit stan modelstan_indices.R             -- establish which indices of the data inform different pieces of the likelihood (to speed up computation) 	

-- Data -- (housed in data/final/)

pp_sp_JP.csv -- site visitsJP_RALU.csv  -- capture data


-- Stan models --

cmr_cort.stan -- Stan model adjusted from parent Stan models (see branch Scientific_Reports)



