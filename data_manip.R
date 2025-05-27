##########################################################
## Build capture history data frame for all populations ##
##########################################################

## number of unique sites, species and pop-spec (my treatment of a "population") 
u_sites <- unique(sampling$Site)
n_sites <- u_sites %>% length()
n_spec  <- unique(sampling$Species) %>% length()

## Find the first and last time period (here months) ever sampled in each population
period_range <- sampling %>% 
  group_by(Site) %>% 
  summarize(
    min_period = min(Month)
  , max_period = max(Month)
  )

## Also, total number of years each population was sampled
year_range <- sampling %>% 
  group_by(Site) %>% 
  summarize(
    n_years = length(unique(Year))
  )

## sampled_periods created from the master sampling data frame 
sampled_periods <- sampling %>% 
  group_by(Site, Year, Month) %>%
  mutate(rep_sec = seq(n())) %>%
  mutate(sampled = 1) %>%
  arrange(Year, Month, SecNumConsec) 

## Just the sampled years as well to be able to grab the correct continuous covariate values
sampled_years <- sampling %>% 
  group_by(Site, Year) %>%
  slice(1) %>%
  dplyr::select(Site, Species, Year)

## Create a data frame of the capture history of each unique individual 
 ## as well as extract individual covariates

## NOTE: only 1 site here, but keeping this loop (which was used for our other CMR
 ## work because it is working correctly as is)
for (i in 1:n_sites) {
  
## Extract a given site and sampling characteristics of that site
period_range.i    <- period_range %>% filter(Site == u_sites[i])
data.i            <- data.all %>% filter(Site == u_sites[i]) %>% 
  mutate(Mark = as.factor(Mark)) %>% mutate(Mark = as.numeric(Mark)) %>%
  mutate(index = seq(n()))
data.i            %<>% group_by(Mark, CaptureDate) %>% slice(1) %>% ungroup()
sampled_periods.i <- sampled_periods %>% filter(Site == u_sites[i])

capt_history.t <- 
expand.grid(
  Month    = seq(from = period_range.i$min_period, to = period_range.i$max_period, by = 1)
, Year     = unique(sampled_periods.i$Year)
, Site = u_sites[i]
, Mark     = unique(data.i$Mark)) %>%
  ## Add species to not screw up multi-species sites
  left_join(., data.i %>% dplyr::select(Mark, Species) %>% distinct()) %>%
  ## Add in which periods were sampled and which individuals were sampled
  left_join(., sampled_periods.i) %>% 
  left_join(., (data.i %>% dplyr::select(
    Month, Year,  Mark, SecNumConsec, Species, bd_load, BdSample, index
  ))) %>% 
  ## drop the times that were never sampled
  filter(!is.na(SecNumConsec)) %>% 
  mutate(sampled = ifelse(is.na(sampled), 0, 1)) %>%
  ## just using a random non-na column to find captures (convenient given how left-join works)
  rename(captured = index) %>% 
  ## convert other NAs to 0s now that we have only sampled days left
  mutate(
    captured = ifelse(is.na(captured), 0, 1)
  , swabbed  = ifelse(is.na(bd_load), 0, 1)) %>%
  mutate(
  ## *** eeek! Need to do better here
    log_bd_load = log(bd_load + 1)
  ## Just because NA are not allowed. This doesn't actually do anything because swabbed == 1 is a logical gate
  , log_bd_load = ifelse(is.na(log_bd_load), 0, log_bd_load)
  ) %>% arrange(Mark, CaptureDate)

## And calculate the duration between each sampling day for each individual
capt_history.t %<>% 
  group_by(Mark) %>% 
  rename(capture_date = CaptureDate) %>%
    ## In reality the 0 should be an NA, but with the way the model is set up the 0 isn't used (and stan cant have NA)
  mutate(capture_gap = c(as.numeric((capture_date - lag(capture_date, 1))[-1]), 0)) 

## Impute missing measured SVL within year
data.i %<>% 
  group_by(Mark, Year) %>% 
  mutate(SVLmm = ifelse(is.na(SVLmm), mean(SVLmm, na.rm = T), SVLmm)) %>% 
  ungroup()

capt_history.t %<>% 
  left_join(., data.i %>% dplyr::select(Mark, Year, CaptureDate, SVLmm) %>%
  rename(capture_date = CaptureDate, len = SVLmm))

## Fill in SVL for all possible capture days
source("pseudo_growth_model_cort.R")

## Average cort and merc if measured multiple times for the same individual. See paper
 ## and supplemental figure for more info on this choice
ind_cov <- data.i %>% 
  group_by(Mark) %>% 
  summarize(
    merc = mean(merc, na.rm = T)
  , cort_base_conc   = mean(cort_base_conc, na.rm = T)
  , cort_stress_conc = mean(cort_stress_conc, na.rm = T)
  ) %>% distinct()

capt_history.t %<>% left_join(., ind_cov)

## Check for multiple different or NA sex entries for cleanup
ind_sex    <- data.i %>% group_by(Mark) %>% summarize(all_sex = unique(Sex))

## only enter in the loop to assign a single sex for those individuals with more than one sex given
sex_to_fix <- ind_sex %>% group_by(Mark) %>% summarize(nsex = n()) %>% filter(nsex > 1) %>% left_join(., ind_sex)
ind_sex    <- ind_sex %>% filter(Mark %notin% unique(sex_to_fix$Mark))

uni_ind_sex_mark <- unique(sex_to_fix$Mark)

if (length(uni_ind_sex_mark) > 0) {

for (ii in 1:length(unique(uni_ind_sex_mark))) {
  ind_sex.temp <- sex_to_fix %>% filter(Mark == uni_ind_sex_mark[ii])
  num_sex      <- apply(as.matrix(sex_entries), 1, FUN = function(x) length(which(ind_sex.temp$all_sex == x))) %>% unlist()
  fill.sex     <- ind_sex.temp[1, ]
  num_sex      <- c(
    f = num_sex[1]
  , m = num_sex[3]
  , u = num_sex[2]
  )
   ## if more than one sex entry is given for an individual, need to assign a sex
  if (sum(num_sex != 0) > 1) {
   ## there are 4 possibilities for greater than 1 entry
     ## 1, 1, 0
     ## 1, 0, 1
     ## 0, 1, 1
     ## 1, 1, 1
   ## Need to assign a choice to each option
    if (num_sex["f"] > 0 & num_sex["m"] > 0 & num_sex["u"] == 0) {
      ## for any animal assigned m, f, and u, it is pretty clear we don't know so give them u
      fill.sex$all_sex <- "U"
    } else if (num_sex["f"] > 0 & num_sex["m"] == 0 & num_sex["u"] > 0) {
      ## for any animal assigned f and u give the benefit of the doubt and go with f
      fill.sex$all_sex <- "F" 
    } else if (num_sex["f"] == 0 & num_sex["m"] > 0 & num_sex["u"] > 0) {
      ## for any animal assigned m and u give the benefit of the doubt and go with m
      fill.sex$all_sex <- "M" 
    } else if (num_sex["f"] > 0 & num_sex["m"] > 0 & num_sex["u"] > 0) {
      ## for any animal assigned m and f go with u
      fill.sex$all_sex <- "U" 
    }
  } 
  if (ii == 1) {
    ind_sex.a <- fill.sex
  } else {
    ind_sex.a <- rbind(ind_sex.a, fill.sex)
  }
}
  
ind_sex <- rbind(ind_sex, ind_sex.a) %>% dplyr::select(-nsex) %>% arrange(Mark) %>% rename(Sex = all_sex)
  
} else {
  
ind_sex %<>% rename(Sex = all_sex)
  
}

capt_history.t %<>% left_join(., ind_sex)

## Drop all individuals that were only ever captured in the final sampling period as these individuals
 ## cannot contribute anything to the model. Rename individuals so their numbers are consecutive
ind_at_end <- (capt_history.t %>% 
  group_by(Mark) %>% 
  mutate(total_caps = cumsum(captured)) %>%
  filter(total_caps == 1) %>% 
  filter(captured == 1) %>% 
  ungroup() %>%
  filter(capture_date == max(capt_history.t$capture_date)) %>%
  summarize(ind_at_end = unique(Mark)))$ind_at_end

if (length(ind_at_end) > 0) {
  capt_history.t %<>% filter(Mark %notin% ind_at_end) %>% droplevels()  
}

## Holdover from larger dataset -- For speed mutate and as.factor is slow
new_marks           <- capt_history.t$Mark %>% as.factor() %>% as.numeric()
capt_history.t$Mark <- new_marks

## Add in other needed indexing columns
capt_history.t %<>% ungroup() %>% 
  arrange(Year, Month, SecNumConsec, Mark, Site) %>% 
  mutate(
## periods counted from the first week in each year
   month_year  = interaction(Month, Year)
 , month_year  = as.factor(month_year)
 , month_year  = as.numeric(month_year)
) %>% relocate(month_year, .after = Year) 

}

capt_history <- capt_history.t; rm(capt_history.t)

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

## Scale length here for cort analysis
capt_history %<>% mutate(len_raw = len, len = scale(len)[, 1])

