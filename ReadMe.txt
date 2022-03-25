CMR models with a latent disease process

Models and code for the empirical data found in the main repo

 --- Open top_level_script.R and work from there. That script will walk you through running anything. All notes needed to run everything are present that script.

Simulation model and scratch development found in: CMR_simulation


---------
Some notes as of March 25
---------


1) "Full" multi-spec/multi-population model working, but yet to be fit on a large scale
----- 

A) Some notes about the model structure:
 -- Uses species-specific fixed effects wherever possible
  -- To do so collapses all RAXX to RALU (expect little success in predicting fixed effect level for this grouping)
 -- Includes pop-spec random effects for most parameters
 -- Imputes length at the species level
 -- Estimates population-average MeHg

B) Some thoughts about model fitting
 -- It appears as if many populations will provide little help for species fixed effects
   -- As in species fixed effects of various forms seem likely to basically always overlap 0 with wide-ish CI
   -- Because of this, it is likely that the pop-spec random effects will be heavily relied upon
 -- It seems that MeHg at the population level is going to provide little information
   -- Likely will have to rely on a second single-population analysis for the MeHg questions
   -- It appears that none of the covariates provide any resolution for MeHg


2) Some next steps 
------

A) A bit of cleanup still needed: 
	- expand.grid in data_manip.R
	- double check the single population model didn't get screwed up with all of the code cleaning for the multiple population model
	- make sure the model with categorical covariates is working correctly when all levels aren't represented
	- make sure the MeHg code in data_covariates.R still leads to appropriate model specification in both the single and multiple population models


3) Bigger Plan
--------------

A) Try and move forward with the species-level fixed effects and pop-spec random effects and step back and re-evaluate if that leads to troubles
B) Fit the "full" model for the main question
C) Fit two other models for sub-questions
 -- 1) Bd over time -- Newts
 -- 2) Relationship between Bd, MeHg and survival -- The well sampled MeHg populations


4) Potential adjustments / additions
------

A) Individual "injury" effect on survival


