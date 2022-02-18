##################################################
## Load each of the data files in the directory ##
##################################################

data.files <- list.files("data")
data.files <- data.files[-which(data.files == "xlsx")]
data.files <- paste("data/", data.files, sep = "")

sampling   <- read.csv("data/xlsx/PP_SP.csv")

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

sampling$Year  <- apply(matrix(sampling$CaptureDate), 1, FUN = function (x) strsplit(x, "/")[[1]][3]) %>% 
  paste("20", ., sep = "") %>% as.numeric()

sampling$Month <- apply(matrix(sampling$CaptureDate), 1, FUN = function (x) strsplit(x, "/")[[1]][1]) %>% 
  as.numeric()

sampling %<>% mutate(CaptureDate = as.Date(date_convert, "%m/%d/%y"))

## can collapse all SubSites and rewrite SecNumConsec
sampling %<>% 
  group_by(Site) %>% 
  arrange(Site, CaptureDate) %>% 
  dplyr::select(-SubSite, -SecNumConsec) %>%
  distinct() %>%
  arrange(Site, CaptureDate, desc(Captures)) %>%
  group_by(Site, Species, CaptureDate) %>%
  slice(1) %>% 
  ungroup(CaptureDate) %>%
  mutate(
    SecNumConsec = seq(n())
  , SubSite      = 1)

for (i in seq_along(data.files)) {

  data.temp <- read.csv(data.files[i])
  
if (strsplit(data.temp$CaptureDate[1], split = " ")[[1]] %>% length() > 1) {
  ddate <- apply(matrix(data.temp$CaptureDate), 1, FUN = function(x) strsplit(x, split = " ")[[1]][1])
  data.temp$CaptureDate <- ddate
}
  
####
## Deal with annoying date formats
####
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

Year <- apply(matrix(date_convert), 1, FUN = function (x) strsplit(x, "/")[[1]][3])
Year <- paste("20", Year, sep = "") %>% as.numeric()

Month <- apply(matrix(date_convert), 1, FUN = function (x) strsplit(x, "/")[[1]][1]) %>% as.numeric()

data.temp$CaptureDate <- date_convert
data.temp      %<>% mutate(
  CaptureDate = as.Date(CaptureDate, "%m/%d/%y")
, Year        = Year
, Month       = Month)

} else {
  
Year <- apply(matrix(data.temp$CaptureDate), 1, FUN = function (x) strsplit(x, "/")[[1]][3])
Year <- paste("20", Year, sep = "") %>% as.numeric()

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
  
data.temp %<>% mutate(bd_load = as.numeric(bd_load)) %>%
  rename(Mark = PitTagCode) %>%
  filter(!is.na(Mark))

## check which years have no swabbing and remove them
no.swabyear <- data.temp %>% group_by(Year) %>% 
  summarize(tot_ss = length(which(!is.na(bd_load)))) %>% 
  filter(tot_ss == 0)

data.temp %<>% filter(Year %notin% no.swabyear$Year) %>% droplevels()

data.temp %<>% dplyr::select(
    Site, SubSite, Species, CaptureDate, Year, Month, PrimNum, SecNumConsec
  , Mark, BdSample, BdResult, SVLmm, MassG, bd_load, HgSampleID) %>%
  mutate(dataset = i)

## Different data sets define SecNumConsec in different ways. Homogenize the choice by making this 
 ## variable a count from 1 to n() and make sure that dates with no captures are still part of the data
sampling.temp <- sampling %>% filter(Site %in% data.temp$Site) %>% 
  dplyr::select(-Notes) %>%
  rename(SecNumConsec_corrected = SecNumConsec)

if ((class(data.temp$SubSite) == "integer" | class(data.temp$SubSite) == "numeric") & 
    class(sampling.temp$SubSite) == "character") {
  sampling.temp %<>% mutate(SubSite = as.numeric(SubSite))
}

if (class(data.temp$SubSite) == "character" & 
    (class(sampling.temp$SubSite) == "integer" | class(sampling.temp$SubSite) == "numeric")) {
  sampling.temp %<>% mutate(SubSite = as.character(SubSite))
}

sampling.temp %<>% dplyr::select(-SubSite)

data.temp %<>% left_join(., sampling.temp)

data.temp %<>% dplyr::select(-SecNumConsec) %>% rename(SecNumConsec = SecNumConsec_corrected)

#adj.SecNumConsec <- data.temp %>% group_by(CaptureDate) %>% 
#  summarize(SecNumConsec = unique(SecNumConsec)) %>% 
#  ungroup() %>%
#  mutate(SecNumConsec_corrected = seq(n()))

#data.temp %<>% left_join(., adj.SecNumConsec) %>% 
#  dplyr::select(-SecNumConsec) %>% 
#  rename(SecNumConsec = SecNumConsec_corrected)

## drop all entries with no mark
data.temp %<>% filter(Mark != "")

## -- this will have to be non-dynamic?? -- ##
## certain species disappear and capture events become mostly opportunistic. Need to figure out what to 
 ## do with these data, but for now just drop them
if (data.temp$Site[1] == "JonesPond" & data.temp$Species[1] == "ANBO") {
  data.temp %<>% filter(Month < 6)
}

if (i == 1) {
  data.all <- data.temp
} else {
  data.all <- rbind(data.all, data.temp)
}
  
}

## If a swab was taken and the result was negative, convert the copies/swab to
 ## zero and not na
data.all %<>% mutate(
  BdResult = ifelse(
    is.na(BdResult) | BdResult == ""
    , "no_samp"
    , BdResult))

data.all[
  data.all$BdSample == "Y"   & 
  data.all$BdResult == "neg" 
    , ]$bd_load <- 0

sampling %<>% dplyr::select(-SubSite, -Notes)

## drop the opportunistic sampling of ANBO from JonesPond
sampling %<>% filter(!(Site == "JonesPond" & Species == "ANBO" & SecNumConsec > 13))

## population defined as a species in a site
data.all %<>% mutate(pop_spec = interaction(Site, Species)) %>% droplevels() %>%
  mutate(pop_spec = factor(pop_spec, levels = unique(pop_spec)))
sampling %<>% mutate(pop_spec = interaction(Site, Species)) %>% droplevels() %>%
  mutate(pop_spec = factor(pop_spec, levels = unique(pop_spec)))


####
## Mercury and Habitat characteristics
####

## Load and subset to the sites that I have data for so far
Oth_hab_cov     <- read.csv("data/xlsx/ARMI_CMR_OtherHabitatCovariates.csv") %>% filter(Site %in% unique(data.all$Site))
temp_precip_hab <- read.csv("data/xlsx/ARMI_CMR_TempPrecip.csv") %>% filter(Site %in% unique(data.all$Site))

## Habitat covaraite details:
# HYDRO	     -- Hydroperiod: Temporary (T) or Permanent (P)
# DRAWDOWN	 -- Drawdown % disappeared, 1 (0-25), 2 (26-50), 3 (51-75), 4 (76-100)
# CANOPY	   -- % Woody canopy cover over wetland 1 (0-25), 2 (26-50), 3 (51-75), 4 (76-100)
# VEG	       -- % sub and emergent vegetation in wetland 1 (0-25), 2 (26-50), 3 (51-75), 4 (76-100)
# SUB	       -- substrate, organic (O) or inorganic (I)
# WCOL	     -- water color: clear (C), brown (B), clear but tannin browned (TB)
# SULF	     -- smell of sulfur: Yes (Y), Sometimes (S), No (N)
# PRODUC	   -- Productivity: Oligo (O), Meso (M), Eutrophic (E)
# Temp and Precip available 2017-2020, named sensibly

## Mercury load 
MeHG            <- read.csv("data/xlsx/CMR_MeHg.csv")

## Checking for errors
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

## Add the Mercury data to the main data frame
 ## For mercury will simply rely on the barcode and not site -- which can take on different names :(
data.all %<>% left_join(.
  , MeHG %>% rename(HgSampleID = Incoming_IDCode) %>% dplyr::select(HgSampleID, MeHg_conc_ppb)
  ) %>% dplyr::select(-HgSampleID)

