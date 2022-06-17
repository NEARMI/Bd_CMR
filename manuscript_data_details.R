##########################################################
## Various details about the dataset for the manuscript ##
##########################################################

library(xtable)

####
## Details on sampling dates
####

sampling   <- read.csv("data/cleaned_cov_csv/PP_SP.csv") 

date_convert <- apply(matrix(sampling$CaptureDate), 1, FUN = function (x) {
  a <- strsplit(x, split = "/")[[1]]
  b <- a[3]
  b <- strsplit(b, "")[[1]]
  if (length(b) > 2) {
  b <- b[c(3, 4)] %>% paste(collapse = "")
  } else {
  b <- paste(b, collapse = "")
  }
  paste(c(a[c(1, 2)], b), collapse = "/")
})

sampling$Year  <- apply(matrix(sampling$CaptureDate), 1, FUN = function (x) strsplit(x, "/")[[1]][3]) %>% as.numeric()
sampling$Month <- apply(matrix(sampling$CaptureDate), 1, FUN = function (x) strsplit(x, "/")[[1]][1]) %>% as.numeric()

sampling %<>% mutate(CaptureDate = as.Date(date_convert, "%m/%d/%y"))

sampling %<>% mutate(pop_spec = interaction(Species, Site)) %>% droplevels() %>%
  mutate(pop_spec = as.character(pop_spec)) %>%
  arrange(pop_spec) %>% 
  mutate(pop_spec = factor(pop_spec, levels = unique(pop_spec)))

sampling %<>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species)
    , to = c(
  "Ambystoma cingulatum"
, "Anaxyrus boreas"
, "Pseudacris maculata"
, "Notophthalmus viridescens"
, "Rana boylii"
, "Rana cascadae"
, "Rana draytonii"
, "Rana luteiventris"
, "Rana pretiosa"
, "Rana sierrae"
    )))

sampling.dates <- sampling %>% group_by(Species, State, Site, Year) %>% 
  summarize(Sampling_Occasions = n_distinct(CaptureDate)) %>% 
  ungroup() %>%
  pivot_wider(c(Species, State, Site), names_from = Year, values_from = Sampling_Occasions)
  
print(xtable(sampling.dates), include.rownames = FALSE)

sampling %>% group_by(Species, State, Site, Year) %>% 
  summarize(Sampling_Occasions = n_distinct(Month)) %>% 
  ungroup() %>%
  pivot_wider(c(Species, Site), names_from = Year, values_from = Sampling_Occasions) %>%
  as.data.frame()

####
## Run all of the data cleaning, but without removing various pieces for the analysis
####

## Modified data loading to not drop things
{
## Loop over all files 
data.files <- list.files("data/cleaned_cmr_csv")
data.files <- paste("data/cleaned_cmr_csv/", data.files, sep = "")

## Load the detailed sampling histor for each location
sampling   <- read.csv("data/cleaned_cov_csv/PP_SP.csv") %>% 
  filter(Notes != "Opportunistic") %>% 
  droplevels()

## Deal with non-date dates 
date_convert <- apply(matrix(sampling$CaptureDate), 1, FUN = function (x) {
  a <- strsplit(x, split = "/")[[1]]
  b <- a[3]
  b <- strsplit(b, "")[[1]]
  if (length(b) > 2) {
  b <- b[c(3, 4)] %>% paste(collapse = "")
  } else {
  b <- paste(b, collapse = "")
  }
  paste(c(a[c(1, 2)], b), collapse = "/")
})

sampling$Year  <- apply(matrix(sampling$CaptureDate), 1, FUN = function (x) strsplit(x, "/")[[1]][3]) %>% as.numeric()
sampling$Month <- apply(matrix(sampling$CaptureDate), 1, FUN = function (x) strsplit(x, "/")[[1]][1]) %>% as.numeric()

sampling %<>% mutate(CaptureDate = as.Date(date_convert, "%m/%d/%y"))

## Loop over all cleaned files
for (i in seq_along(data.files)) {

## load the file
data.temp  <- read.csv(data.files[i])

## deal with dates
if (strsplit(data.temp$CaptureDate[1], split = " ")[[1]] %>% length() > 1) {
  ddate <- apply(matrix(data.temp$CaptureDate), 1, FUN = function(x) strsplit(x, split = " ")[[1]][1])
  data.temp$CaptureDate <- ddate
}
  
if (length(grep("/", data.temp$CaptureDate[1])) > 0) {

date_convert <- apply(matrix(data.temp$CaptureDate), 1, FUN = function (x) {
  a <- strsplit(x, split = "/")[[1]]
  b <- a[3]
  b <- strsplit(b, "")[[1]]
  if (length(b) > 2) {
  b <- b[c(3, 4)] %>% paste(collapse = "")
  } else {
  b <- paste(b, collapse = "")
  }
  paste(c(a[c(1, 2)], b), collapse = "/")
})

Year  <- apply(matrix(date_convert), 1, FUN = function (x) strsplit(x, "/")[[1]][3])
Year  <- paste("20", Year, sep = "") %>% as.numeric()
Month <- apply(matrix(date_convert), 1, FUN = function (x) strsplit(x, "/")[[1]][1]) %>% as.numeric()

data.temp$CaptureDate <- date_convert
data.temp      %<>% mutate(
  CaptureDate = as.Date(CaptureDate, "%m/%d/%y")
, Year        = Year
, Month       = Month)

} else {
  
Year  <- apply(matrix(data.temp$CaptureDate), 1, FUN = function (x) strsplit(x, "/")[[1]][3])
Year  <- paste("20", Year, sep = "") %>% as.numeric()
Month <- apply(matrix(data.temp$CaptureDate), 1, FUN = function (x) strsplit(x, "/")[[1]][1]) %>% as.numeric()
  
data.temp      %<>% mutate(
  CaptureDate = as.Date(CaptureDate)
, Year        = Year
, Month       = Month)
  
}
  
####
## A few other modifications
####

## name changes for convenience
if ("TrgetCopies.swb" %in% names(data.temp)) {
  data.temp %<>% rename(bd_load = TrgetCopies.swb) 
} else if ("TargetCopies.swab" %in% names(data.temp)) {
  data.temp %<>% rename(bd_load = TargetCopies.swab) 
} else {
  
}
  
## Convert bd load to a numeric and drop all individuals that were never marked as they can't be used to inform the model
 ## *** For now also drop individuals found dead. Could add a flag in the model to indicate that we know that an individual
  ## died at some point between the last capture and the current capature of the dead individual, but this is kind of a hassle
   ## for just a couple of individuals and definitely will not have an impact
data.temp %<>% mutate(bd_load = as.numeric(bd_load)) %>%
  rename(Mark = IndividualID) %>%
  filter(!is.na(Mark), Mark != "") %>%
  filter(dead != 1 | is.na(dead))

## *** check which years have no swabbing and remove them as this is just going to complicate the model structure and
 ## isn't going to help resolve bd vs survival (and is rare--maybe just one or two years total across all populations)
no.swabyear <- data.temp %>% group_by(Year) %>% 
  summarize(tot_ss = length(which(!is.na(bd_load)))) %>% 
  filter(tot_ss == 0)

data.temp %<>% filter(Year %notin% no.swabyear$Year) %>% droplevels()

## Create a subset of sampling to link up with the data to get the correct prim num and sec num
sampling.temp <- sampling %>% filter(Site %in% data.temp$Site)

## Deal with the mix of numeric and non-numeric subsites
if ((class(data.temp$SubSite) == "integer" | class(data.temp$SubSite) == "numeric") & 
    class(sampling.temp$SubSite) == "character") {
  sampling.temp %<>% mutate(SubSite = as.numeric(SubSite))
}

if (class(data.temp$SubSite) == "character" & 
    (class(sampling.temp$SubSite) == "integer" | class(sampling.temp$SubSite) == "numeric")) {
  sampling.temp %<>% mutate(SubSite = as.character(SubSite))
}

## Link up the PrimNum and SecNumConsec from "sampling" with data.temp
data.temp %<>% left_join(., sampling.temp %>% dplyr::select(Site, SubSite, Species, CaptureDate, PrimNum, SecNumConsec))

## Some datasets have extra (and different) data columns than other datasets. Extract the necessary 
 ## columns from each (already processed so that all files contain at least these)
data.temp %<>% dplyr::select(
    Site, SubSite, Species, CaptureDate, Year, Month, PrimNum, SecNumConsec
  , Mark, BdSample, BdResult, SwabLost, SVLmm, MassG, Sex, Age, potential_injury_effect, bd_load, HgSampleID
  , flagged, reason, Notes 
  ) %>%
  mutate(dataset = i)

## *** The current strategy is to collapse all SubSites and just use the main site. 
 ## Which SubSites are sampled is still used to calculate each animals detection probability
  ## Recalculate SecNumConsec as well just as a count from 1:n() of days each site is sampled
sampling.temp %<>% 
  rename(SecNumConsec_corrected = SecNumConsec) %>%
  arrange(Site, CaptureDate) %>% 
  dplyr::select(-SubSite, -SecNumConsec_corrected, -Notes, -Region, -State) %>%
  distinct() %>%
  arrange(Site, CaptureDate, desc(Captures)) %>%
  group_by(Site, Species, CaptureDate) %>%
  slice(1) %>% 
  ungroup(CaptureDate) %>%
  mutate(SecNumConsec_corrected = seq(n())) 

data.temp %<>% left_join(., sampling.temp) %>% dplyr::select(-SecNumConsec) %>% rename(SecNumConsec = SecNumConsec_corrected)

## *** There are a few "opportunistic" sampling events in JonesPond when other work was being done. Drop these as the
 ## population is not the same as it is on-season which adds a wrinkle to the model that isn't worth dealing with
if (data.temp$Site[1] == "JonesPond" & data.temp$Species[1] == "ANBO") {
  data.temp %<>% filter(Month < 6)
}

## *** Some locations (for just a few animals) have "duplicate" swabs (the same individual swabbed multiple times)
 ## For now collapsing these measures to a single measure per animal in a sensible way 
  ## (details given in this ~50 line section)
if (data.temp %>% filter(reason == "duplicate") %>% nrow() > 0) {

data.temp.d <- data.temp %>% filter(reason == "duplicate")
data.temp %<>% filter(reason != "duplicate")
  
## Probably can do this in a more efficient way, but not clearly aware of one
data.temp.d.ind <- unique(data.temp.d$Mark)

for (ind_i in seq_along(data.temp.d.ind)) {
  data.temp.d.t <- data.temp.d %>% filter(Mark == data.temp.d.ind[ind_i])
  
  ## Take the positive swab if negative and positive
  if ("neg" %in% data.temp.d.t$BdResult & "pos" %in% data.temp.d.t$BdResult) {
    data.temp.d.t %<>% filter(BdResult == "pos")
  } else if (all(data.temp.d.t$BdResult == "neg")) {
  ## Take negative if both negative
    data.temp.d.t <- data.temp.d.t[1, ]
  } else if (all(data.temp.d.t$BdResult == "pos")) {
  ## Take the average of the positives if both positive
    temp_bd_load          <- mean(data.temp.d.t$bd_load)
    data.temp.d.t         <- data.temp.d.t[1, ]
    data.temp.d.t$bd_load <- temp_bd_load
  } else {
  ## Shouldn't be a 4th option so break
    print("Check entries for the duplicates, there appears to be a problem")
    break
  }
  
  if (ind_i == 1) {
    data.temp.d.a <- data.temp.d.t
  } else {
    data.temp.d.a <- rbind(data.temp.d.a, data.temp.d.t)
  }
  
}

data.temp <- rbind(data.temp, data.temp.d.a)
  
}

if (i == 1) {
  data.all <- data.temp
} else {
  data.all <- rbind(data.all, data.temp)
}
  
}

daily_hab_covar <- sampling %>% dplyr::select(Region, Site, SubSite, Species, CaptureDate)

## *** Adjust the master "sampling" file to match the choices for the data (i.e., dropping subsites as I am doing for now)
sampling %<>% 
  group_by(Site) %>% 
  arrange(Site, CaptureDate) %>% 
  dplyr::select(-SubSite, -SecNumConsec) %>%
  distinct() %>%
  arrange(Site, CaptureDate, desc(Captures)) %>%
  group_by(Site, Species, CaptureDate) %>%
  slice(1) %>% 
  ungroup(CaptureDate) %>%
  mutate(SecNumConsec = seq(n()))

## population defined as a species in a site
data.all %<>% mutate(pop_spec = interaction(Species, Site)) %>% droplevels() %>%
  mutate(pop_spec = as.character(pop_spec)) %>%
  arrange(pop_spec) %>% 
  mutate(pop_spec = factor(pop_spec, levels = unique(pop_spec)))
sampling %<>% mutate(pop_spec = interaction(Species, Site)) %>% droplevels() %>%
  mutate(pop_spec = as.character(pop_spec)) %>%
  arrange(pop_spec) %>% 
  mutate(pop_spec = factor(pop_spec, levels = unique(pop_spec)))

####
## Habitat characteristics and MeHg
####

## Load and subset to the sites that I have data for so far
Oth_hab_cov     <- read.csv("data/cleaned_cov_csv/ARMI_CMR_OtherHabitatCovariates.csv") %>% filter(Site %in% unique(data.all$Site))
temp_precip_hab <- read.csv("data/cleaned_cov_csv/ARMI_CMR_TempPrecip.csv") %>% filter(Site %in% unique(data.all$Site))
MeHG            <- read.csv("data/cleaned_cov_csv/CMR_MeHg.csv") 

## Add some columns converting binned quantiative values into continuous metrics
Oth_hab_cov %<>% mutate(
  drawdown_cont = plyr::mapvalues(DRAWDOWN, from = c(1, 2, 3, 4), to = c(12.5, 37.5, 62.5, 87.5))
, canopy_cont   = plyr::mapvalues(CANOPY, from = c(1, 2, 3, 4), to = c(12.5, 37.5, 62.5, 87.5))
, veg_cont      = plyr::mapvalues(VEG, from = c(1, 2, 3, 4), to = c(12.5, 37.5, 62.5, 87.5))
)

## Edit daily_hab_covar dataframe to get covaraties on each day for detection
daily_hab_covar %<>% left_join(., Oth_hab_cov) %>% 
  group_by(Site, CaptureDate) %>%
  summarize(
    drawdown      = mean(DRAWDOWN, na.rm = T) %>% round()
  , veg           = mean(VEG, na.rm = T) %>% round()
  , drawdown_cont = mean(drawdown_cont, na.rm = T)
  , veg_cont      = mean(veg_cont, na.rm = T)
  )

## 
  
## Data frame to look for some MeHg errors
Me_Hg_error_check <- data.all %>% 
  ungroup() %>% 
  group_by(Site) %>% 
  summarize(num_with_MeHgID = length(which(HgSampleID != ""))) %>%
  left_join(.
  , data.all %>% dplyr::select(Site, Mark, HgSampleID) %>%
  left_join(
    .
  , (MeHG %>% rename(HgSampleID = Incoming_IDCode) %>% dplyr::select(HgSampleID, MeHg_conc_ppb))
  ) %>% filter(!is.na(MeHg_conc_ppb)) %>% group_by(Site) %>% 
  droplevels() %>%
  summarize(num_with_measured_MeHg = n()))

## Add the Mercury data to the main data frame (rely on the barcode)
data.all %<>% left_join(.
  , MeHG %>% rename(HgSampleID = Incoming_IDCode) %>% dplyr::select(HgSampleID, MeHg_conc_ppb)
  ) %>% dplyr::select(-HgSampleID)

data.all %<>% relocate(c(flagged, reason, Notes), .after = MeHg_conc_ppb)
data.all %<>% relocate(pop_spec, .after = Species)

## Quick check how often each covariate is measured in each population was measured 
cov_meas <- data.all %>%
  group_by(Site) %>%
  summarize(
    prop_merc   = length(which(!is.na(MeHg_conc_ppb))) / n()
  , prop_weight = length(which(!is.na(MassG))) / n()
  , prop_length = length(which(!is.na(SVLmm))) / n()
  , prop_sex    = length(which(!is.na(Sex))) / n()
  , prop_age    = length(which(!is.na(Age))) / n()
  )

## store all of the types of entries for sex for use in cleaning later
data.all    %<>% mutate(Sex = plyr::mapvalues(Sex, from = "", to = "U"))

sex_entries <- unique(data.all$Sex)
}
red_ind    <- FALSE
## Modified data manipulation to not drop things
{
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

capt_history.t <- 
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
capt_history.t %<>% group_by(Mark) %>% 
  rename(capture_date = CaptureDate) %>%
    ## In reality the 0 should be an NA, but with the way the model is set up the 0 isn't used (and stan cant have NA)
  mutate(capture_gap = c(as.numeric((capture_date - lag(capture_date, 1))[-1]), 0)) 

## *** There are some duplicate entries. Rare but does happen. Probably either poor QA/QC or possibly individuals caught more than once in a day.
 ## Mostly these are identical entries duplicated so just take the first
capt_history.t %<>% group_by(Mark, SecNumConsec) %>% slice(1)

## If using less than the full population for debug purposes take those random individuals here
#if (red_ind_PA_debug) {
#  which_ind      <- (capt_history.t %>% group_by(Mark) %>% filter(swabbed == 1) %>% 
#      summarize(nswabs = n()) %>% arrange(desc(nswabs)) %>% slice(1:num_ind))$Mark
#  capt_history.t %<>% filter(Mark %in% which_ind) %>% droplevels()
#} else if (red_ind & !red_ind_PA_debug) {
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

## *** Also deal with individual sex, for which there can be lots of possibilities because of misidentification. Don't want to deal with
 ## a whole extra complication in the model itself (misidentification) so doing some pre-processing here to settle on a sex for each individual
ind_sex          <- data.i %>% group_by(Mark) %>% summarize(all_sex = unique(Sex))

## only enter in the loop to assing a single sex for those individuals with more than one sex given
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
}

capt_history %<>%  mutate(
  Species = plyr::mapvalues(Species, from = unique(Species)
    , to = c(
  "Ambystoma cingulatum"
, "Anaxyrus boreas"
, "Pseudacris maculata"
, "Notophthalmus viridescens"
, "Rana boylii"
, "Rana cascadae"
, "Rana draytonii"
, "Rana luteiventris"
, "Rana pretiosa"
, "Rana sierrae"
    ))) 

recaps <- capt_history %>% group_by(Species, State, Site, Year, Mark) %>% 
  filter(captured == 1) %>%
  summarize(recaps = sum(captured)) %>% 
  mutate(recaptured = ifelse(recaps > 1, 1, 0)) %>% 
  ungroup(Mark) %>%
  summarize(
    inds_capt  = n_distinct(Mark)
  , recapt_ind = sum(recaptured)) %>% 
  ungroup() 

recaps %<>% dplyr::select(-inds_capt) %>% pivot_wider(c(Species, State, Site), names_from = Year, values_from = recapt_ind)
  
print(xtable(recaps), include.rownames = FALSE)



