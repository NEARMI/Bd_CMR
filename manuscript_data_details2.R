##########################################################
## Various details about the dataset for the manuscript ##
##########################################################

####
## Details on sampling dates
####

## Table S1 ------

sampling   <- read.csv("data/final/pp_sp_new.csv") %>% filter(Model == 1) %>% droplevels() %>% mutate(CaptureDate = as.Date(SurveyDate, "%m/%d/%y")) %>%
  dplyr::select(-SurveyDate)

sampling %<>% dplyr::select(-Species) %>% 
  rename(Site = MasterSite, SubSite = FriendlySiteName, Species = SPEC) %>%
  dplyr::select(Site, SubSite, Species, CaptureDate, SecNumConsec, Captures)

sampling %<>% mutate(
  Year   = apply(matrix(as.character(CaptureDate)), 1, FUN = function (x) strsplit(x, "[-]")[[1]][1]) %>% as.numeric()
, Month  = apply(matrix(as.character(CaptureDate)), 1, FUN = function (x) strsplit(x, "[-]")[[1]][2]) %>% as.numeric()
, julian = as.POSIXlt(CaptureDate)$yday
) %>% mutate(
  Week = ceiling(julian/7)
) %>%
  relocate(c(Year, Month, Week, julian), .after = CaptureDate)

sampling %>% group_by(Species, Site, Year) %>% 
  summarize(Sampling_Occasions = n_distinct(CaptureDate)) %>% 
  ungroup() %>%
  pivot_wider(c(Species, Site), names_from = Year, values_from = Sampling_Occasions)

sampling %>% group_by(Species, Site, Year) %>% 
  summarize(Sampling_Occasions = n_distinct(Month)) %>% 
  ungroup() %>%
  pivot_wider(c(Species, Site), names_from = Year, values_from = Sampling_Occasions) %>%
  as.data.frame()

## Table S2 ------

recaps <- capt_history %>% group_by(Species, Site, Year, Mark) %>% 
  filter(captured == 1) %>%
  summarize(recaps = sum(captured)) %>% 
  mutate(recaptured = ifelse(recaps > 1, 1, 0)) %>% 
  ungroup(Mark) %>%
  summarize(
    inds_capt  = n_distinct(Mark)
  , recapt_ind = sum(recaptured)) %>% 
  ungroup() 

recaps %>% dplyr::select(-recapt_ind) %>% 
  pivot_wider(c(Species, Site), names_from = Year, values_from = inds_capt)

recaps %>% dplyr::select(-inds_capt) %>% 
  pivot_wider(c(Species, Site), names_from = Year, values_from = recapt_ind)

## Table 1 ------

# 1 -- Total Individuals Captured
# 2 -- Total Captures
# 3 -- Percent of Individuals Recaptured
# 4 -- Total Bd Swabs
# 5 -- Percent of Captured Individuals Swabbed for Bd
# 6 -- Percent of Captured Individuals Recaptured After Being Swabbed
# 7 -- Percent of Captured Individuals Measured for MeHg

# 1, 2, 4
capt_history %>% 
  group_by(pop_spec) %>%
  summarize(
    tot_ind = n_distinct(Mark)
  , tot_cap = sum(captured)
  , tot_swa = sum(swabbed)
  )

# 3, 5, 7
capt_history %>% group_by(pop_spec, Mark) %>% 
  filter(captured == 1) %>%
  summarize(
    recaps = sum(captured)
  , swabs  = sum(swabbed)
  , mehg   = sum(merc, na.rm = T)) %>% 
  mutate(
    recaptured = ifelse(recaps > 1, 1, 0)
  , swabbs     = ifelse(swabs > 0, 1, 0)) %>% 
  ungroup(Mark) %>%
  summarize(
    inds_capt   = n_distinct(Mark)
  , recapt_ind  = sum(recaptured)
  , swabbed_ind = sum(swabbs)
  , merced_ind  = length(which(mehg != 0))) %>% 
  ungroup() %>% 
  mutate(
    prop_recap = recapt_ind / inds_capt
  , prop_bd    = swabbed_ind / inds_capt
  , prop_merc  = merced_ind / inds_capt
    ) %>% 
  left_join(
    .
    #6
  , capt_history %>% 
  filter(captured == 1) %>%
  group_by(pop_spec, Mark) %>%
  summarize(
    capt_then_swabbed = ifelse(any(which(captured == 1) > min(which(swabbed == 1))), 1, 0)
  ) %>% filter(
    capt_then_swabbed == 1
  ) %>% ungroup() %>%
  group_by(pop_spec) %>%
  summarize(tot_capt_then_swabbed = sum(capt_then_swabbed))
  ) %>% mutate(
    prop_capt_then_swabbed = tot_capt_then_swabbed / inds_capt
  )


  
