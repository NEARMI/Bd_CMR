CMR models with a latent disease process. Only files relevant to the models fit for the publication in Scientific Reports.

All project development scripts are found in the main branch of the repository.
The files in this branch include, but are not limited to, code for many alternative models and code for simulating CMR data and fitting to simulated data (for model dev)
ReadMe provided in the subfolder



 -- Fits can be run on a local machine by opening ``pipeline_for_publication.R'' and working from there. That script will walk you through running anything. All notes needed to run everything are present that script.

 -- Script can also be streamlined a bit to run on a remote machine. For some older script examples see the Yeti folder in scr_debug_scripts



-- Scripts -- (alphabetically, for the order that scripts are used open and work through "pipeline_for_publication.R")


pipeline_for_publication.R      -- top level script from which all others are called and analysis is run

complete_data.R			-- load data through 2022 data_covariates.R		-- load and clean site and individual-level covariatesdata_manip.R			-- organize data into the long form with zeros for individuals not captureddata_stan.R			-- organize data into form for the Stan modelsdata_temp_all.R			-- load and clean daily site-level temperature datadataset_notes.R			-- some info on the dataestablishing_mm.R		-- build the [m]odel [m]atrices for the Stan model fitggplot_theme.R			-- ggplot theme for beautified plotsmanuscript_data_details.R	-- exports the info for the manuscript tablesmanuscript_figures.R		-- lots of figures for the manuscript (mostly supplemental)packages_functions.R		-- needed packages and functionsstan_fit_mm.R			-- fit [m]odel [m]atrix style stan model for multiple populationsstan_indices.R			-- establish which indices of the data inform different pieces of the likelihood (to speed up computation)
	

-- Data -- (housed in data/final/)

Hab_Cov.csv		-- categorical habitat variable information for each site
Temp_Precip.csv		-- temperature and precipitation summarized to yearly averages by site
all_cmr_years.csv	-- all CMR data for all cites
multi_swab.csv		-- data from individuals swabbed multiple times on a single capture
pp_sp_new.csv		-- all site visits in all populations 
raw_temp.csv		-- raw daily temperature in all sites


-- Stan models --

See readme files in the stan_current/ folder



