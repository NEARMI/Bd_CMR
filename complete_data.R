####
## Load and clean the compiled data
####

## Load required data and make a few minor adjustments
cmr        <- read.csv("data/final/all_cmr_years.csv") %>% mutate(CaptureDate = as.Date(SurveyDate, "%m/%d/%y")) %>%
  dplyr::select(-SurveyDate)
sampling   <- read.csv("data/final/pp_sp_new.csv") %>% filter(Model == 1) %>% droplevels() %>% mutate(CaptureDate = as.Date(SurveyDate, "%m/%d/%y")) %>%
  dplyr::select(-SurveyDate)
sites      <- sampling %>% dplyr::select(MasterSite, SiteName, FriendlySiteName, SPEC, CaptureDate, Model) %>% distinct()

## Make a few changes to sites that were lost or merged or changed name
 ## SubSite Adjustments
  ## SMNWR0103        -- DROP (abandoned)
  ## SMNWR1002        -- DROP (abandoned)
  ## SMNWR0145        -- Roll into SMNWR1601 (merged)
  ## SMNWR0124        -- Roll into SMNWR3001 (merged)
  ## San Francisquito -- SFC (to match previous)
cmr        %<>% mutate(
  FriendlySiteName = plyr::mapvalues(FriendlySiteName
    , from = c("SMNWR0145", "SMNWR0124", "San Francisquito"), to = c("SMNWR1601", "SMNWR3001", "SFC"))
, SiteName = plyr::mapvalues(SiteName
    , from = c("SMNWR1601", "SMNWR3001", "SMNWR0145", "SMNWR0124"), to = c("SMNWR_E", "SMNWR_W", "SMNWR_E", "SMNWR_W"))
  ) %>% filter(FriendlySiteName %notin% c("SMNWR0103", "SMNWR1002"))

## Drop all abandoned sites and those with far too few captures as well as any dates in any sites where ALL CAPTURED animals were not marked
 ## (one other dates SOME animals were not marked, these removed separately)
cmr %<>% left_join(., sites) %>% relocate(c(MasterSite, SPEC), .before = SiteName) %>% filter(!is.na(MasterSite)) %>% dplyr::select(-Model)

## Drop all animals that were not marked
cmr %<>% filter(AnimalMark.ID %notin% c("notMarked")) %>% filter(!is.na(AnimalMark.ID)) %>% droplevels()

## Remove columns not needed for analysis and rename columns to match what I used to develop the code and model
sampling %<>% dplyr::select(-Species) %>% 
  rename(Site = MasterSite, SubSite = FriendlySiteName, Species = SPEC) %>%
  dplyr::select(Site, SubSite, Species, CaptureDate, SecNumConsec, Captures)

cmr %<>% dplyr::select(-Species) %>% 
  rename(Site = MasterSite, SubSite = FriendlySiteName, Species = SPEC
    , Mark = AnimalMark.ID, BdSample = Swab_barcode, BdResult = Bd
    , SVLmm = BodyLength, bd_load = Ct_Bd, merc = TissueMeHg) %>%
  dplyr::select(Site, SubSite, CaptureDate, Species, Mark, BdSample, BdResult, Sex, SVLmm, bd_load, merc)

## Break up date into component pieces
cmr %<>% mutate(
  Year   = apply(matrix(as.character(CaptureDate)), 1, FUN = function (x) strsplit(x, "[-]")[[1]][1]) %>% as.numeric()
, Month  = apply(matrix(as.character(CaptureDate)), 1, FUN = function (x) strsplit(x, "[-]")[[1]][2]) %>% as.numeric()
, julian = as.POSIXlt(CaptureDate)$yday
) %>% relocate(c(Year, Month, julian), .after = CaptureDate)

sampling %<>% mutate(
  Year   = apply(matrix(as.character(CaptureDate)), 1, FUN = function (x) strsplit(x, "[-]")[[1]][1]) %>% as.numeric()
, Month  = apply(matrix(as.character(CaptureDate)), 1, FUN = function (x) strsplit(x, "[-]")[[1]][2]) %>% as.numeric()
, julian = as.POSIXlt(CaptureDate)$yday
) %>% relocate(c(Year, Month, julian), .after = CaptureDate)

## Tidy up a number of data columns
cmr %<>% mutate(
  ## Convert the Bd Ct to load
  bd_load = ct_to_load(as.numeric(bd_load))
  ## Fix placeholder lengths
, SVLmm = ifelse(SVLmm == -99.0, NA, SVLmm)
  ## Make sure there are only 3 entries for Sex
, Sex = plyr::mapvalues(Sex, from = c("Unknown", "Male", "Female", ""), to = c("U", "M", "F", "U"))
  ## Convert Tissue MeHg to full body MeHg
, merc = ifelse(Species %in% c("NOVI", "AMCI"), scale_newt(merc), scale_frog(merc))
  ) %>% mutate(
  ## if negative add 0 load
  bd_load = ifelse(BdResult == "neg", 0, bd_load)
  )

sex_entries <- unique(cmr$Sex)

## Convert all RAXX to RANA
cmr      %<>% mutate(Species = plyr::mapvalues(Species, from = c("RALU", "RACA", "RAPR", "RADR", "RASI", "RABO"), to = rep("RANA", 6)))
sampling %<>% mutate(Species = plyr::mapvalues(Species, from = c("RALU", "RACA", "RAPR", "RADR", "RASI", "RABO"), to = rep("RANA", 6)))

## Save a data frame for which subsites were visited on each day for the detection model
sampling_for_p <- sampling

## Save a data frame that contains a description of the habitat structure in each SubSite sampled on that day (for detection in the model)
daily_hab_covar <- sampling %>% dplyr::select(Site, SubSite, Species, CaptureDate)

####
## Collapse SubSites 
####

## Collapses so that every visit to the main site gets sorted in order regardless of what subsites are visited
 ## ^^ data frame above stores which subsites are visited on each of these days to inform detection
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

## Define the pop * spec interactions. Used as the levels of the random effects
cmr %<>% 
  mutate(pop_spec = interaction(Species, Site)) %>% droplevels() %>%
  mutate(pop_spec = as.character(pop_spec)) %>%
  arrange(pop_spec) %>% 
  mutate(pop_spec = factor(pop_spec, levels = unique(pop_spec))) %>%
  relocate(pop_spec, .before = Site)

sampling %<>% 
  mutate(pop_spec = interaction(Species, Site)) %>% droplevels() %>%
  mutate(pop_spec = as.character(pop_spec)) %>%
  arrange(pop_spec) %>% 
  ## make sure sampling pop_spec levels match the cmr data
  mutate(pop_spec = factor(pop_spec, levels = unique(cmr$pop_spec))) %>%
  relocate(pop_spec, .before = Site)

## Add the sampling event number to the main CMR data 
cmr %<>% left_join(., sampling %>% dplyr::select(-Captures)) %>% 
  relocate(SecNumConsec, .after = CaptureDate) %>%
  arrange(pop_spec, CaptureDate)

## rename to match code written with 2018-2021 data
data.all <- cmr; rm(cmr)

## Double check (and correct any minor errors) for sites with days with no captures
print(
  paste(
  "All Captures/No Captures checks out?"
, 
sampling %>% left_join(
 . , data.all %>% group_by(pop_spec, SubSite, CaptureDate) %>% slice(1) %>% dplyr::select(pop_spec, CaptureDate) %>%
  mutate(some_capt = 1) 
) %>% mutate(some_capt = ifelse(is.na(some_capt), 0, 1)) %>%
  mutate(capt_miss = Captures - some_capt) %>% 
  filter(capt_miss != 0) %>% nrow() == 0
, sep = " ---- ")
)

####
## Habitat characteristics and MeHg
####

## Load and subset to the sites that I have data for so far
Oth_hab_cov     <- read.csv("data/final/Hab_Cov.csv") %>% filter(Site %in% unique(data.all$Site))
temp_precip_hab <- read.csv("data/final/Temp_Precip.csv") %>% filter(Site %in% unique(data.all$Site))

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

