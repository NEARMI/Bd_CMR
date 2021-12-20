##############################################
## Examine the data from across the country ##
##############################################

####
## Notes December 20:
####

## Some initial thoughts from data exploration
 ## 1) A bit unclear what needs to happen with subsites where animals move between subsites (but maybe only between a few)
  ## -- combine all subsites or treat them as unique entities
 ## 2) A majority of individuals in some sites do test negative, which could throw a wrench into the problem of no discrete states
 ## 3) Still need to figure out how to pull in all of the sites with the no captures
 ## 4) OVERALL there really isn't much trend over month and if there is (if you squint a bit) it is somewhat variable by
  ## location, so definitely will need to rely on temp, the data of which is yet to be seen

### packages and functions

needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan")
lapply(needed_packages, require, character.only = TRUE)
source("../ggplot_theme.R")
set.seed(10002)
'%notin%' <- Negate('%in%')


### load and combine data

data.files <- list.files("data")
data.files <- data.files[-which(data.files == "xlsx")]
data.files <- paste("data/", data.files, sep = "")

sampling   <- read.csv("data/xlsx/PP_SP.csv")

for (i in seq_along(data.files)) {

  data.temp <- read.csv(data.files[i])
  
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
  , Mark, BdSample, BdResult, SVLmm, MassG, bd_load) %>%
  mutate(dataset = i)

## Different data sets define SecNumConsec in different ways. Homogenize the choice by just making this 
 ## variable a count from 1 to n()

adj.SecNumConsec <- data.temp %>% group_by(CaptureDate) %>% 
  summarize(SecNumConsec = unique(SecNumConsec)) %>% 
  ungroup() %>%
  mutate(SecNumConsec_corrected = seq(n()))

data.temp %<>% left_join(., adj.SecNumConsec) %>% 
  dplyr::select(-SecNumConsec) %>% 
  rename(SecNumConsec = SecNumConsec_corrected)

## drop all entries with no mark
data.temp %<>% filter(Mark != "")

## -- this will have to be non-dynamic?? -- ##
## certain species disappear and capture events become mostly opportunisitc. Need to figure out what to 
 ## do with these data, but for now just drop them
#if (i == 1) {
#  data.temp %<>% filter(Month < 7)
#}

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
    , BdResult)
)

data.all[
  data.all$BdSample == "Y"   & 
  data.all$BdResult == "neg" 
    , ]$bd_load <- 0


data.all %<>% mutate(Species = plyr::mapvalues(Species,
  from = "BUBO", to = "ANBO"))

### Some exploration of sampling effort, distribution etc.

data.all %>% group_by(Site) %>%
  summarize(length(unique((SubSite))))

## Individuals seen across 2 or 3 sites at the most... Treat these
 ## as a single site is one extreme
data.all %>% 
  filter(Site == "Blackrock") %>% 
  group_by(Mark) %>%
  summarize(nsites = length(unique(SubSite))) %>%
  arrange(desc(nsites)) %>% head(20)

data.all %>% 
  group_by(Site, SubSite, Mark) %>% 
  summarize(num_capt = n()) %>% 
  ungroup(Mark, SubSite) %>% 
  summarize(
    one_capt  = length(which(num_capt == 1)) / n()
  , two_capt  = length(which(num_capt == 2)) / n()
  , more_capt = length(which(num_capt > 2)) / n())

data.all %>% 
  group_by(Species, Mark) %>% 
  summarize(num_capt = n()) %>% 
  ungroup(Mark) %>% 
  summarize(
    one_capt  = length(which(num_capt == 1)) / n()
  , two_capt  = length(which(num_capt == 2)) / n()
  , more_capt = length(which(num_capt > 2)) / n())

data.all %>% 
  group_by(Site, SubSite, Mark) %>% 
  summarize(num_capt = n()) %>% 
  ungroup(Mark) %>% {
    ggplot(., aes(x = num_capt)) + 
      geom_histogram(bins = 20) +
      facet_wrap(~Site, scales = "free")
  }

data.all %>% 
  group_by(Site, SubSite, Mark) %>% 
  summarize(num_swab = length(which(BdSample == "Y"))) %>% 
  ungroup(Mark, SubSite) %>% 
  summarize(
    no_swab    = length(which(num_swab == 0)) / n()
  , one_swab   = length(which(num_swab == 1)) / n()
  , two_swab   = length(which(num_swab == 2)) / n()
  , more_swabs = length(which(num_swab > 2)) / n())

data.all %>% 
  mutate(
    Month = as.numeric(Month)
  , Year  = as.factor(Year)) %>%
  group_by(Site, Month, Year) %>% 
  summarize(num_swab = length(which(BdSample == "Y"))) %>% {
    ggplot(., aes(Month, num_swab)) + 
      geom_line(aes(colour = Year)) +
      geom_point(aes(colour = Year)) +
      facet_wrap(~Site) +
      ylab("Total Swabs")
  }


### Bd and other plots

data.all %>% filter(bd_load > 0) %>%
  mutate(Month = as.factor(Month)) %>% {
    ggplot(., aes(Month, bd_load)) + 
    geom_boxplot(width = 0.3) +
    scale_y_log10() +
    facet_wrap(~Site, scales = "free")
  }

data.all %>%
  mutate(Month = as.factor(Month)) %>% {
    ggplot(., aes(Month, bd_load)) + 
    geom_violin() +
    scale_y_log10() +
    facet_wrap(~Site, scales = "free")
  }

data.all %>%
  group_by(Month, Year, Site) %>%
  filter(BdSample == "Y") %>% 
  summarize(
    test_pos = length(which(bd_load > 0)) / n()
  , nswabs  = n()) %>% {
    ggplot(., aes(Month, test_pos)) + 
    geom_line(aes(colour = as.factor(Year))) +
    geom_point(aes(colour = as.factor(Year), size = nswabs)) +
    scale_size_continuous(breaks = c(10, 30, 50, 80, 110, 140, 170)) +
    facet_wrap(~Site)
}

data.all %>% filter(bd_load > 0) %>% {
    ggplot(., aes(Month, bd_load)) + 
    geom_point(aes(colour = as.factor(Year))) +
    geom_line(aes(group = Mark)) + 
    scale_y_log10() +
    facet_wrap(~Site, scales = "free")
}

data.all %>% filter(Site == "Y029") %>% {
    ggplot(., aes(Month, bd_load)) + 
    geom_point(aes(colour = as.factor(Year))) +
    geom_line(aes(group = Mark)) + 
    scale_y_continuous(
      trans  = "pseudo_log"
    , breaks = c(1, 10, 100, 1000, 10000, 1E6, 1E7)) +
    facet_wrap(~Year, scales = "free")
}

data.all %>% {
    ggplot(., aes(Month, bd_load)) + 
    geom_point(aes(colour = as.factor(Year))) +
    geom_line(aes(group = Mark)) + 
    scale_y_continuous(
      trans  = "pseudo_log"
    , breaks = c(1, 10, 100, 1000, 10000, 1E6, 1E7)) +
    facet_wrap(~Species, scales = "free")
}

data.all %>% filter(Species == "RALU") %>% {
    ggplot(., aes(Month, bd_load)) + 
    geom_point(aes(colour = as.factor(Year))) +
    geom_line(aes(group = Mark)) + 
    scale_y_continuous(
      trans  = "pseudo_log"
    , breaks = c(1, 10, 100, 1000, 10000, 1E6, 1E7)) +
    facet_wrap(~Year, scales = "free")
}



