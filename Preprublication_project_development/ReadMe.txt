CMR models with a latent disease process

Models and code for the empirical data found in the main repo

 -- Fits can be run on a local machine by opening ``top_level_script.R'' and working from there. That script will walk you through running anything. All notes needed to run everything are present that script.

 -- Script can also be streamlined a bit to run on a remote machine. For some older script examples see the Yeti folder in scr_debug_scripts


Some simulation data and model fits to these data (built very early in the project) and scratch development can be found in ``/CMR_simulation''


----------

* All needed scripts to run from "top_level_script.R" are in the main folder
* All of the most up to date Stan models are in stan_current/ (dev_stan/ holds many old and scratch Stan model iterations, many early upon which I built and many broken)

----------

-- Scripts -- (alphabetically, for the order that scripts are used open and work through "top_level_script.R")


top_level_script.R		-- top level script from which all others are called and analysis is run

capt_plot_multi.R		-- plot CMR data for multiple populationscapt_plot.R			-- plot CMR data for single populationcomplete_data.R			-- load data through 2022 data_covariates.R		-- load and clean site and individual-level covariatesdata_manip.R			-- organize data into the long form with zeros for individuals not captureddata_stan.R			-- organize data into form for the Stan modelsdata_temp_all.R			-- load and clean daily site-level temperature datadataset_notes.R			-- some info on the datadetermine_model.R		-- pick the Stan file to run based on user choicesestablishing_mm.R		-- build the [m]odel [m]atrices for the Stan model fitexplore_data.R			-- explore the data a bitexport_coefs.R			-- export coefficients from fits for plottingggplot_theme.R			-- ggplot theme for beautified plotsmanuscript_data_details.R	-- exports the info for the manuscript tablesmanuscript_figures.R		-- lots of figures for the manuscript (mostly supplemental)packages_functions.R		-- needed packages and functionsplotting_facet_names.R		-- expand shorthand species and site names for better looking plotsplotting_mehg.R			-- plot model for Bd-MeHg only modelplotting_multipop.R		-- plotting for multi-population fitplotting.R			-- plotting for single population fitProject_Notes.R			-- some project notesstan_fit_mehg_only.R		-- fit Bd-MeHg only model (no survival)stan_fit_mm.R			-- fit [m]odel [m]atrix style stan model for multiple populationsstan_fit_single.R		-- fit stan model for single populationstan_indices.R			-- establish which indices of the data inform different pieces of the likelihood (to speed up computation)
	

-- Data -- (housed in data/final/)

Hab_Cov.csv		-- categorical habitat variable information for each site
Temp_Precip.csv		-- temperature and precipitation summarized to yearly averages by site
all_cmr_years.csv	-- all CMR data for all cites
multi_swab.csv		-- data from individuals swabbed multiple times on a single capture
pp_sp_new.csv		-- all site visits in all populations 
raw_temp.csv		-- raw daily temperature in all sites


-- Stan models --

See readme files in the stan_current/ folder



