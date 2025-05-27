####
## Load and clean data
####

## Load required data and make a few minor adjustments
cmr        <- read_excel("data/final/JP_RALU.xlsx") %>% 
  mutate(CaptureDate = as.Date(SurveyDate)) %>%
  dplyr::select(-SurveyDate, -SubSite)
sampling   <- read.csv("data/final/pp_sp_JP.csv") %>% 
  filter(Model == 1) %>% droplevels() %>% mutate(CaptureDate = as.Date(SurveyDate, "%m/%d/%y")) %>%
  dplyr::select(-SurveyDate)
sites      <- sampling %>% dplyr::select(MasterSite, SiteName, FriendlySiteName, SPEC, CaptureDate, Model) %>% distinct()

## Bit of site cleanup
cmr %<>% left_join(., sites) %>% 
  relocate(c(MasterSite, SPEC), .before = SiteName) %>% 
  filter(!is.na(MasterSite)) %>% dplyr::select(-Model)

## Double check to make sure that all animals were marked 
cmr %<>% filter(AnimalMark.ID %notin% c("notMarked")) %>% 
  filter(!is.na(AnimalMark.ID)) %>% droplevels()

## Remove columns not needed for analysis and rename columns to match what we used
 ## to develop the code and model for the original CMR paper  
sampling %<>% dplyr::select(-Species) %>% 
  rename(Site = MasterSite, SubSite = FriendlySiteName, Species = SPEC) %>%
  dplyr::select(Site, SubSite, Species, CaptureDate, SecNumConsec, Captures)

cmr %<>% dplyr::select(-Species) %>% 
  rename(Site = MasterSite, SubSite = FriendlySiteName, Species = SPEC
         , Mark = AnimalMark.ID, BdSample = Swab_barcode, BdResult = Bd
         , SVLmm = BodyLength, merc = MeHg_conc_ppb, bd_load = TargetCopies) %>%
  dplyr::select(
    Site, SubSite, CaptureDate, Species, Mark, BdSample, BdResult, Sex, SVLmm
  , bd_load, merc, cort_base_conc, cort_stress_conc)

## Break date up into component pieces
cmr %<>% mutate(
    Year   = apply(matrix(as.character(CaptureDate)), 1, FUN = function (x) strsplit(x, "[-]")[[1]][1]) %>% as.numeric()
  , Month  = apply(matrix(as.character(CaptureDate)), 1, FUN = function (x) strsplit(x, "[-]")[[1]][2]) %>% as.numeric()
  , julian = as.POSIXlt(CaptureDate)$yday
) %>% mutate(
  Week = ceiling(julian/7)
) %>%
  relocate(c(Year, Month, Week, julian), .after = CaptureDate)

sampling %<>% mutate(
    Year   = apply(matrix(as.character(CaptureDate)), 1, FUN = function (x) strsplit(x, "[-]")[[1]][1]) %>% as.numeric()
  , Month  = apply(matrix(as.character(CaptureDate)), 1, FUN = function (x) strsplit(x, "[-]")[[1]][2]) %>% as.numeric()
  , julian = as.POSIXlt(CaptureDate)$yday
) %>% mutate(
  Week = ceiling(julian/7)
) %>%
  relocate(c(Year, Month, Week, julian), .after = CaptureDate)

## Tidy up a number of data columns
cmr %<>% mutate(
  ## Fix placeholder lengths if there are any
    SVLmm = ifelse(SVLmm == -99.0, NA, SVLmm)
  ## Make sure there are only 3 entries for Sex
  , Sex = plyr::mapvalues(Sex, from = c("Unknown", "Male", "Female", ""), to = c("U", "M", "F", "U"))
  ## Convert Tissue MeHg to full body MeHg
  , merc = ifelse(Species %in% c("NOVI", "AMCI"), scale_newt(merc), scale_frog(merc))
) %>% mutate(
  ## if negative add 0 load
  bd_load = ifelse(BdResult == "neg", 0, bd_load)
)

## Double check number of unique sex vals
sex_entries <- unique(cmr$Sex)
if (n_distinct(cmr$Sex) > 3) {
  stop("Check Sex Labels")
}

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

## Add the sampling event number to the main CMR data 
cmr %<>% left_join(., sampling %>% dplyr::select(-Captures)) %>% 
  relocate(SecNumConsec, .after = CaptureDate) %>%
  arrange(Site, CaptureDate)

## rename to match code written with 2018-2021 data
data.all <- cmr; rm(cmr)
