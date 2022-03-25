##########################################################
## Build capture history data frame for all populations ##
##########################################################

## number of unique sites, species and pop-spec (my treatment of a "population") 
u_sites <- unique(sampling$pop_spec)
n_sites <- u_sites %>% length()
n_spec  <- unique(sampling$Species) %>% length()

## Find the first and last time period (here months) ever sampled in each population
period_range <- sampling %>% 
  group_by(pop_spec) %>% 
  summarize(
    min_period = min(Month)
  , max_period = max(Month)
  )

## Also, total number of years each population was sampled
year_range <- sampling %>% 
  group_by(pop_spec) %>% 
  summarize(
    n_years = length(unique(Year))
  )

## sampled_periods created from the master sampling data frame 
sampled_periods <- sampling %>% 
  group_by(Site, Year, Month) %>%
  mutate(rep_sec = seq(n())) %>%
  mutate(sampled = 1) %>%
  arrange(Year, pop_spec, Month, SecNumConsec) 

## Just the sampled years in each population as well to be able to grab the correct continuous covariate values
sampled_years <- sampling %>% 
  group_by(Site, pop_spec, Year) %>%
  slice(1) %>%
  dplyr::select(Site, Species, pop_spec, Year)

## Create a data frame of the capture history of each unique individual in each population as well as extract individual covariates

all_ind <- 0 ## establish a counter that is used at the bottom of the loop to count all individuals across all populations
for (i in 1:n_sites) {
  
## Extract a given site and sampling characteristics of that site
period_range.i    <- period_range %>% filter(pop_spec == u_sites[i])
data.i            <- data.all %>% filter(pop_spec == u_sites[i]) %>% mutate(Mark = as.factor(Mark)) %>% mutate(Mark = as.numeric(Mark))
sampled_periods.i <- sampled_periods %>% filter(pop_spec == u_sites[i])

## *** Pretty sure this can be simplified dropping Month and Year and just using unique(Mark) and unique(CaptureDate)
 ## Will come back to this as it is just cleanup and wont change anything
capt_history.t <- 
  ## First create that "all possible combinations" data frame (all secondary periods in which each animal
   ## could possibly have been caught)
expand.grid(
  Month    = seq(from = period_range.i$min_period, to = period_range.i$max_period, by = 1)
, Year     = unique(sampled_periods.i$Year)
, pop_spec = u_sites[i]
, Mark     = unique(data.i$Mark)) %>% 
  ## Add species to not screw up multi-species sites
  left_join(., data.i %>% dplyr::select(Mark, Species) %>% distinct()) %>%
  ## Add in which periods were sampled and which individuals were sampled
  left_join(., sampled_periods.i) %>% 
  left_join(., (data.i %>% dplyr::select(Month, Year,  Mark, SecNumConsec, Species, bd_load, BdSample))) %>% 
  ## drop the times that were never sampled
  filter(!is.na(SecNumConsec)) %>% 
  mutate(sampled = ifelse(is.na(sampled), 0, 1)) %>%
  ## just using a random non-na column to find captures (convenient given how left-join works)
  rename(captured = BdSample) %>% 
  ## convert other NAs to 0s now that we have only sampled days left
  mutate(
    captured = ifelse(is.na(captured), 0, 1)
  , swabbed  = ifelse(is.na(bd_load), 0, 1)) %>%
  mutate(
  ## *** eeek! Need to do better here
    log_bd_load = log(bd_load + 1)
  ## Just because NA are not allowed. This doesn't actually do anything because swabbed == 1 is a logical gate
  , log_bd_load = ifelse(is.na(log_bd_load), 0, log_bd_load)
  )

## And calculate the duration between each sampling day for each individual
capt_history.t %<>% group_by(Mark, Year) %>% 
  rename(capture_date = CaptureDate) %>%
    ## In reality the 0 should be an NA, but with the way the model is set up the 0 isn't used (and stan cant have NA)
  mutate(capture_gap = c(as.numeric((capture_date - lag(capture_date, 1))[-1]), 0)) 

## *** There are some duplicate entries. Rare but does happen. Probably either poor QA/QC or possibly individuals caught more than once in a day.
 ## Mostly these are identical entries duplicated so just take the first
capt_history.t %<>% group_by(Mark, SecNumConsec) %>% slice(1)

## If using less than the full population for debug purposes take those random individuals here
if (red_ind) {
  which_ind      <- sample(unique(capt_history.t$Mark), min(length(unique(capt_history.t$Mark)), num_ind))
  capt_history.t %<>% filter(Mark %in% which_ind) %>% droplevels()
}

## Before converting Mark to a numeric, find the individual specific covariates to be used later
ind_cov <- data.i %>% group_by(Mark) %>% 
  summarize(
    merc = mean(MeHg_conc_ppb, na.rm = T)
    ## returns the one value or a mean if caught multiple times in one primary period
  , len  = mean(as.numeric(SVLmm), na.rm = T)
  , inj  = ifelse(any(potential_injury_effect == 1), 1, 0)
  ) %>% distinct()

capt_history.t %<>% left_join(., ind_cov)

## Drop all individuals that were only ever captured in the final sampling period as these individuals cannot contribute anything to the model
 ## Have to also then rename individuals so their numbers are consecutive
ind_at_end <- (capt_history.t %>% 
  group_by(Mark) %>% 
  mutate(total_caps = cumsum(captured)) %>%
  filter(total_caps == 1) %>% 
  filter(captured == 1) %>% 
  ungroup() %>%
  filter(capture_date == max(capture_date)) %>%
  summarize(ind_at_end = unique(Mark)))$ind_at_end

capt_history.t %<>% filter(Mark %notin% ind_at_end) %>% droplevels()

## For speed for the very largest (mutate and as.factor is really slow)
new_marks           <- capt_history.t$Mark %>% as.factor() %>% as.numeric()
capt_history.t$Mark <- new_marks

## Jump through a minor hoop to name unique individuals
n_inds         <- max(capt_history.t$Mark)
capt_history.t %<>% mutate(Mark = Mark + all_ind)
all_ind        <- all_ind + n_inds

## Add in other needed indexing columns
capt_history.t %<>% ungroup() %>% 
  arrange(Year, Month, SecNumConsec, Mark, Site) %>% 
  mutate(
## periods counted from the first week in each year
   month_year  = interaction(Month, Year)
 , month_year  = as.factor(month_year)
 , month_year  = as.numeric(month_year)
) %>% relocate(month_year, .after = Year) %>% 
  dplyr::select(-Notes)

if (i == 1) {
capt_history <- capt_history.t
} else {
capt_history <- rbind(capt_history, capt_history.t)
}

print("--------------")
print(paste("Through", i, "of", n_sites, "sites", sep = " "))
print(paste(all_ind, "Total Individuals Organized", sep = " "))
print(paste("Uses", object.size(capt_history)/1E9, "Gb", sep = " "))
print("--------------")

}

## Sort data frame in the appropriate order (counting through consecutive weeks one individual at a time)
capt_history %<>% arrange(Mark, Year, Site, Month, SecNumConsec)

## individuals' measured bd 
capt_history.bd_load <- capt_history %>% 
  ungroup() %>%
  arrange(Mark, month_year) %>%
  filter(swabbed == 1)

## first and last _OF THE CAPTURE EVENTS_ in which each individual was captured
 ## (possible min and max will vary by which population individuals are in)
capture_range  <- capt_history %>% 
  group_by(Mark) %>% 
  filter(sampled == 1) %>%  
  summarize(
    first = min(which(captured == 1))
  , final = max(which(captured == 1))) %>% 
  dplyr::select(first, final) %>% 
  ## Remove all individuals in the data set that were never captured (in case there are any left over for w/e reason,
   ## but they should all have been removed already)
  filter(!is.infinite(first) | !is.infinite(final))
