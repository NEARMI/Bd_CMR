##################################################
## Load each of the data files in the directory ##
##################################################

## Loop over all files 
data.files <- list.files("data/cleaned_cmr_csv")
data.files <- paste("data/cleaned_cmr_csv/", data.files, sep = "")

## Load the detailed sampling histor for each location
sampling   <- read.csv("data/cleaned_cov_csv/PP_SP.csv") %>% 
## remove the 2 individual SMNWR_E-NOVI site. It is not helpful to have these data in the model
  filter(!(Site == "SMNWR_E" & Species == "NOVI")) %>%
## remove the minimal sampling of the Blackrock population in the late season
  filter(Notes != "Opportunistic") %>% 
  droplevels()

## *** For now (as of March 25) convert all RAXX to RANA to help [the jury is still out] fixed effect estimates
RANA_spec <- unique(sampling$Species)[grep("RA", unique(sampling$Species))]
sampling  %<>% mutate(Species = plyr::mapvalues(Species, from = RANA_spec, to = rep("RANA", length(RANA_spec))))

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
check.rana <- unique(data.temp$Species)[grep("RA", unique(data.temp$Species))]

## *** if RAXX convert to RANA
if (length(check.rana) > 0) {
 data.temp  %<>% mutate(Species = plyr::mapvalues(Species, from = check.rana, to = rep("RANA", length(check.rana))))   
}
  
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

## Build a data frame that contains a description of the habitat structure in each
 ## SubSite sampled on that day (for detection in the model)
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
