##########################################################
## Build capture history data frame for all populations ##
##########################################################

## number of unique sites (populations) in the combined dataset
n_sites <- unique(data.all$Site) %>% length()
u_sites <- unique(data.all$Site)

## Find the first and last time period ever sampled in each population
period_range <- data.all %>% 
  group_by(Site) %>% 
  summarize(
    min_period = min(Month)
  , max_period = max(Month)
  )

## Also, each year each population was sampled
year_range <- data.all %>% 
  group_by(Site) %>% 
  summarize(
    n_years = length(unique(Year))
  )

## Find the number of secondary periods in each of the primary periods
sampled_periods <- data.all %>% 
  group_by(Site, Year) %>%
  summarize(Month = unique(Month)) %>%
  left_join(.
    , data.all %>% 
  group_by(Site, Year, Month) %>%
  summarize(SecNumConsec = unique(SecNumConsec))
    ) %>% 
  group_by(Site, Year, Month) %>%
  mutate(rep_sec = seq(n())) %>%
  mutate(sampled = 1) %>%
  arrange(Year, Site, Month, SecNumConsec) 

## Create a data frame of the caputre history of each unique individual in each population as
 ## well as extract individual covariates
all_ind <- 0
for (i in 1:n_sites) {
  
## Extract a given site and sampling characteristics of that site
period_range.i    <- period_range %>% filter(Site == u_sites[i])
data.i            <- data.all %>% filter(Site == u_sites[i])
sampled_periods.i <- sampled_periods %>% filter(Site == u_sites[i])
  
capt_history.t <- 
  ## First create that "all possible combinations" data frame (all secondary periods in which each animal
   ## could possibly have been caught)
expand.grid(
  Month = seq(from = period_range.i$min_period, to = period_range.i$max_period, by = 1)
, Year  = unique(data.i$Year)
, Site  = u_sites[i]
, Mark  = unique(data.i$Mark)) %>% 
  ## Add in which periods were sampled and which individuals were sampled
  left_join(., sampled_periods.i) %>% 
  left_join(., (data.i %>% dplyr::select(Month, Year,  Mark, SecNumConsec, Species, bd_load))) %>% 
  ## drop the times that were never sampled
  filter(!is.na(SecNumConsec)) %>% 
  mutate(sampled = ifelse(is.na(sampled), 0, 1)) %>%
  rename(
  ## just using a random non-na column to find captures (convenient given how left-join works)
    captured = Species) %>% 
  ## convert other NAs to 0s now that we have only sampled days left
  mutate(
    captured = ifelse(is.na(captured), 0, 1)
  , swabbed  = ifelse(is.na(bd_load), 0, 1)) %>%
  mutate(
  ## eeek! Need to do better here
    log_bd_load = log(bd_load + 1)
  ## Just because NA are not allowed. This doesn't actually do anything because swabbed == 1 is a logical gate
  , log_bd_load = ifelse(is.na(log_bd_load), 0, log_bd_load)
  )

## There are some duplicate entries? (Individuals caught more than once in a day??)
# capt_history.t %>% group_by(Mark, SecNumConsec) %>% summarize(n_entry = n()) %>% arrange(desc(n_entry))
capt_history.t %<>% group_by(Mark, SecNumConsec) %>% slice(1)

## Before converting Mark to a numeric, find the individual specific covariates to be used later
if ("MeHgConc" %in% names(data.all)) {

ind_cov <- data.i %>% group_by(Mark, Year, Month) %>% 
  summarize(
    merc = mean(MeHgConc, na.rm = T)
    ## returns the one value or a mean if caught multiple times in one primary period
  , size = mean(as.numeric(MassG), na.rm = T)  
    ) %>% distinct()

} else {
  
ind_cov <- data.i %>% group_by(Mark, Year, Month) %>% 
  summarize(
   size = mean(as.numeric(MassG), na.rm = T)  
    ) %>% distinct()
  
}

capt_history.t %<>% left_join(., ind_cov)

## Jump through a few hoops to name unique individuals
 ## NOTE: this is an issue if individuals move populations (as that individual in each population will be
  ## labeled as a differnet individual). Not sure what to do about it though
capt_history.t %<>% mutate(Mark = as.factor(Mark)) %>%
  mutate(Mark = as.numeric(Mark))
n_inds <- max(capt_history.t$Mark)
capt_history.t %<>% mutate(Mark = Mark + all_ind)
all_ind <- all_ind + n_inds

## Add in other needed indexing columns
capt_history.t %<>% ungroup() %>% 
  arrange(Year, Month, SecNumConsec, Mark, Site) %>% 
  mutate(
## periods counted from the first week in each year
   month_year  = interaction(Month, Year)
 , month_year  = as.factor(month_year)
 , month_year  = as.numeric(month_year)
)

if (i == 1) {
capt_history <- capt_history.t
} else {
capt_history <- rbind(capt_history, capt_history.t)
}

print("")
print("")
print(paste("Through", i, "of", n_sites, "sites", sep = " "))

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
  ## Remove all individuals in the data set that were never captured (in case there are any for w/e reason)
  filter(!is.infinite(first) | !is.infinite(final))
