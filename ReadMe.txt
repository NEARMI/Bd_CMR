CMR models with a latent disease process

Models and code for the empirical data found in the main repo

 -- Fits can be run on a local machine by opening ``top_level_script.R'' and working from there. That script will walk you through running anything. All notes needed to run everything are present that script.

 -- Fits can be run on a remote machine with slurm using scripts that read ``XXXX_yeti.R''. These wrapper scripts have less commenting than ``top_level_script.R''
   -- Planned fits can be found in model_runs.numbers


Simulation model and scratch development found in ``/CMR_simulation''


---------
Some notes as of May 6
---------


1) Current status of work
----- 

A) Single and Multiple populations finished for now. Many run on Yeti, others waiting to be run. Where I am leaving this project for about 10 days:
	-- Most individual models were run, but FL populations failed for some reason. Need to determine why
	-- Plotting still not working for a number of populations. Need to retrieve the Rds files from PC and work through the plotting
	-- Single NOVI continuous models running, but didn't check if they were going to finish or fail with timeout
	-- Multi-population jobs running, but didn't check if they were going to finish or fail with timeout
	-- Overleaf methods a bit behind the current fitting strategy (population groupings)
	-- After all of this need to make some final decisions on what to run

B) Some notes on fitting process (3500 samples and three chains on Yeti)
	-- All individual populations apart from PA and MA NOVI run within 6h
	-- ANBO populations job running (but expected to finish within 12h)
	-- RANA populations job running (with fitting MeHg) expected to finish 24-36h
	-- All pops together without the few NOVI pops running, unclear about time expectation (maybe 36h)

C) Notes on fits
	-- Small populations consistently are predicted to have increasing survival with increasing Bd
		-- adding Bd to detection doesn't seem to help these small populations
	-- Most RANA populations have a very strong positive relationship between length and survival
		-- These has too narrow of CIs to be reasonable. Need to try length in detection

D) To Do
	1) Plot, debug, and upload some fits to Overleaf
	2) Clean up overleaf language, upload some new figures from fits, and share
	3) Another new model needed for SMNWR_E AMCI because no individuals were found infected
	4) Try length and Bd in detection for these RANA populations (with the MeHg single species multipop model)
	