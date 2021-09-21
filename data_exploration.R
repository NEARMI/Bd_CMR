##############################
## Explore the NE newt data ##
##############################

## !! Note, data sets not included in the public repo

####
## Data Columns in Bd_Newts_AllSites.csv
####

## SA           - overall study area (pennsylvania, Massachusetts, Wisconsin)
## Site         - name of the pond where the newt was caught
## PondLocation - This is a location within a pond from Evan's site (we have something similar for ours but I don't see it in this file)
## Recap        - was it previously caught and marked
## Mark         - the individual ID for each newt (it is a color pattern of marks - WI may not have individual marks).
## SVL/TL       - length measures SVL is snout-vent-length and TL for total length. Typically use SVL for herps since tail damage can make total length noisy (we also took tail samples for mercury, which affect total length)
## Swab         - was the animal sampled for chytrid - for PA we have a lot of captures without swabs
## SwabID/SwabBatch
  ##  /SwabCode - these are IDs used on the sample tubes and by the lab when the samples are processed
## Tail         - was the animal sampled for mercury (a much smaller subset of animals)
## TailCode     - the ID for the tube with the mercury sample used by the lab
## Status       - was Bd detected by the lab 
## status.2     - I think this is just status converted to 0 and 1
## copies       - the relative amount of Bd detected (I forget the units here - other studies will use zoospore equivalents - either way it is a relative amount of dna on the swab sample)

## HP           - (hydroperiod) relative number of days a pond is inundated with water in a year (short versus long)
## Abundance    - estimated number of newts in the pond from the mark-recapture data
## Area         - estimated water surface area of the pond
## Density      - abundance/area

## The rest of the fields are different measures of water or air temperature (the fields differ in how sensors were averaged and the number of days that are averaged to get the reading)

####
## Data notes
####

## From Jill on Sep 10

 ## 1) Part of the problem for 2019 may be that we swabbed all amphibians captured (lots of tadpoles, etc) but only marked newts, so some of the swab codes I see in the spreadsheet are actually data from spring peepers, etc.
 ## 2) Some of the ‘NA’s in 2020 and 2021 were because we captured more than we had the ability to process on a few days, so there are instances where an animal was released without any processing, an animal was swabbed but not marked, and/or was marked but not swabbed. We still included those in our data entry so that we would know the # and locations of captures even if we didn’t have the other information.
 ## 3) Also, there was some confusion with the marks since I last sent Jim our data I believe. There were some double marks (2 individuals with the same mark) that we’ve since fixed either in the field or with database clues (1 was a male, 1 female. etc). 

####
## Data notes as of Sep 20
####

## 1) In the cleaned up data set some individuals seem to be unidentified ? (Mark == NA) and still have Bd taken,
  ## especially in WI. Also in WI, most individuals were never recaptured but bd swabs are almost 100%. 
   ## So maybe WI will just be good to capture the distribution of Bd throughout the year?

## 2) Load follows fraction positive quite closely. Ramps up from April through June and then starts falling
 ## Is this because of temperature?
 ## Are individuals infected with very low doses early in the season that isn't caught OR
  ## are individuals actually free of infection and become infected and then load increases 
  ## i.e. increasing load leads to infections in other individuals?

## 3) Load within specific individuals---if you squint pretty hard---could maybe vaguely follow
 ## a quadratic, with a peak around June 1. 
  ## However, individuals are super variable and really go in all which ways
 ## It may be sensible to try and model a general trend in Bd following day, and have individuals deviate
  ## randomly from the pop mean at any given day. To do this could maybe have a latent process of
   ## the mean, and then have each individual have a random deviate from which they differ from the mean?
    ## Coded something like observation error, but don't allow each individual to have a latent bd process?

## 4A) Some individuals are clearly being mis-marked as not infected when they are infected given 
 ## abrupt jumps from sick to not sick back to sick again (with reasonably high copy numbers)

## 4B) Bd load (copy number) also has lots of observation noise. Definitely need an observation error
 ## model for Bd... e.g., see individual Mark ROXXXR and BRBYXX and Red tail VIE

## 5) Quite a bit of variation in Bd between sites in the same state (e.g., see EC and ML in WI)

## 6) Seems that MA and WI could sample later in the season

## 7) Can we use known relationships between Bd and temperature and the temperature in these
 ## locations to constrain our fits?

####
## Data issues (as of Sep 20)
####

## 1) Five entries with Status == "Negative" but copies.swab > 0
## 2) Unclear about newts vs other species in MA
## 3) 

####
## Summary notes about patterns thinking about a model
####

## 1) It is not very helpful that the place with good bd sampling doesn't track individuals well
 ## -- In this location it doesn't seem that Bd load changes much over the sampling periods
 ## -- It also seems in WI there is more site variation than year variation or season variation
##    OR that the locations that sample individuals well don't sample Bd all that often (MA)
 ## -- In this location there is maybe a small quadratic relationship in bd load over season
## 2) Overall, bd dynamics seems to vary quite a lot among states and among sites within those states
 ## -- In theory one overarching capture model could contain all locations and the bd submodel could vary by location:site?
 ##    However, this seems like quite a pain / impossible [?] to fit
## 3) There does seem to be some small quadratic-like pattern in changes in bd load overall at the level of sites
 ## -- Though individuals are highly variable in their trajectories (when their bd load increases, falls, or if it mostly stays constant)
 ## -- This points to what will probably be a big difficulty modeling latent individual trajectories in Bd
 ##    what will probably need to happen is a seasonal model for bd with individual-specific deviates that are constant over the season
 ##    (i.e., individual specific slopes are probably quite impractical)
## 4) Much of the increasing pattern in Bd load early in the season is individuals becoming infected and thus
 ##   bd load increasing from 0 to some non-zero value (instead of always being infected with a low detectable level in the early season)
 ## -- This may point to the need for some sort of joint model of infection and then bd growth, but that
 ##    may be unreasonable to fit
## 5) Both positive infection and copy number need an observation model because both are measured
 ##   with a decent amount of error (in some individuals copy plummets to 0 after a high swab then returns to a high value in a few days)
 ## -- May need to be a bit complicated: i.e., detecting positivity depends on load. That is,
 ##    a small load may lead to a zero estimate 

####
## For most model thoughts see data_model_notes.numbers (or .csv), but...
####

## 1) How to write the model by SA:site over many years? 
 ## random effects by location?
 ## for time, how to deal with the fact that sampling occurs in uneven intervals?
  ## how often is a specific site revisited? Can we write down every date and just fill in 
  ## the 1s for when the individuals are caught? That is, if we didn't have a capture of an 
  ## individual at a specific date we can stick in a zero? 
   ## Some individuals move sites, how do we deal with this in a small/simple model?
    ## For example can we cluster sampling by the week across all ponds (or something like that)?

## 2) For a model of load there will need to also be a model of infected/not (I think)
 ## to differentiate changes in load from the process of becoming infected

####
## Packages and functions
####
needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan")
lapply(needed_packages, require, character.only = TRUE)
source("../ggplot_theme.R")
set.seed(10002)
'%notin%' <- Negate('%in%')

## Data from Montana? -- put this aside for now
# RALUCaptDat            <- readxl::read_excel("fwcmrcode/RALUCaptDat.xlsx")
# names(RALUCaptDat)[28] <- "bd_copies"
# RALUCaptDat            %<>% mutate(bd_copies = as.numeric(bd_copies))

# RALUCaptDat %>% filter(Year == 2018, Month == 7) %>% {
#   ggplot(., aes(x = bd_copies)) +
#     geom_histogram(bins = 50)
# }

## Data from WI, MA, and PA
# Bd_Newts_AllSites <- read.csv("Bd_Newts_AllSites.csv")
Bd_Newts_AllSites   <- read.csv("Bd_Newts_AllSites_9.12.21.csv")

## Stupid dates in R
if (length(grep("/", Bd_Newts_AllSites$Date[1])) > 0) {

date_convert <- apply(matrix(Bd_Newts_AllSites$Date), 1, FUN = function (x) {
  a <- strsplit(x, split = "/")[[1]]
  b <- a[3]
  b <- strsplit(b, "")[[1]][c(3, 4)] %>% paste(collapse = "")
  paste(c(a[c(1, 2)], b), collapse = "/")
})

Bd_Newts_AllSites$Date <- date_convert
Bd_Newts_AllSites      %<>% mutate(Date = as.Date(Date, "%m/%d/%y"))

} else {
  
Bd_Newts_AllSites      %<>% mutate(Date = as.Date(Date))
  
}

## data characteristics
head(Bd_Newts_AllSites)

Bd_Newts_AllSites %>% summarize_all(n_distinct)

## Capture distribution
cap_dist <- Bd_Newts_AllSites %>% 
  group_by(Mark, Site, SA) %>%
  summarize(
    num_captures = n()
  , num_swabs    = length(which(Swab == TRUE))) %>% 
  arrange(desc(num_captures)) %>% 
  filter(!is.na(Mark))                 ## There seems to be a lot of individuals with no mark notated 

cap_dist %>% ungroup() %>% summarize_all(n_distinct)
cap_dist %>% group_by(SA) %>% summarize_all(n_distinct)

## What does Mark look like by state?
Bd_Newts_AllSites %>% 
  group_by(SA) %>% 
  filter(!is.na(Mark)) %>%
  slice(1:5) %>% 
  summarize(Mark = unique(Mark))

## SEP 10: Ok, so some individuals without a mark are able to be known as a recapture? How?
 ## UPDATE SEPT 20: These have been cleared up with the new data set
Bd_Newts_AllSites %>% filter(is.na(Mark)) %>% 
  summarize(no_recap = length(which(Recap == TRUE)) / n())

## Some individuals do move sites, but not all that many
cap_dist %>% ungroup(Site) %>% summarize(num_sites = n()) %>% {
    ggplot(., aes(x = num_sites)) + 
      geom_histogram(bins = 100) +
      xlab("Number of Sites where an Individual was caught")
  }
  
## Clearly more prolific sites or more sampled sites. Distributions don't seem that different by site?
cap_dist %>% filter(!is.na(Site)) %>% {
    ggplot(., aes(x = num_captures)) + 
      geom_histogram(bins = 100) + 
      facet_wrap(~Site, scales = "free")
}

## Look at time sampling of bd in all individuals 
 ## Unfortunately most individuals are never swabbed or swabbed only once
  ## Some reasonable[?] number of individuals are swabbed twice or three times but very few
   ## are swabbed more times than that. It could be hard with these data to model changes in bd load
    ## (or even infection status) over time
Bd_Newts_AllSites %>% 
  group_by(Mark) %>% 
  summarize(bd_prop = length(which(Swab == TRUE))) %>%
  arrange(desc(bd_prop)) %>% 
  filter(!is.na(Mark)) %>% {
    ggplot(., aes(x = bd_prop)) +
    geom_histogram(bins = 100) +
      xlab("bd swabs by individual") 
  }

Bd_Newts_AllSites %>% 
  group_by(Mark) %>% 
  summarize(bd_prop = length(which(Swab == TRUE)) / n()) %>% {
    ggplot(., aes(x = bd_prop)) +
    geom_histogram(bins = 100) +
      xlab("Proportion of individual captures with bd swabs") 
  }

## Of all of the captures, bd was sampled extremely well in WI and quite well in MA, but
 ## somewhat poorly in PA. Unfortunately, in WI few of the captures are marked and recaptures are low
Bd_Newts_AllSites %>% 
  group_by(SA) %>%
  summarize(
    mark        = length(which(!is.na(Mark)))
  , no_mark     = length(which(is.na(Mark)))
  , recap       = length(which(Recap == TRUE))
  , no_recap    = length(which(Recap != TRUE))
  , prop_mark   = 1 - (length(which(is.na(Mark))) / n())
  , prop_recap  = length(which(Recap == TRUE)) / n()
  , bd_prop = length(which(Swab == TRUE)) / n())

## And in the individuals captured many times?
Bd_Newts_AllSites %>% filter(Mark %in% {
  (cap_dist %>% filter(num_captures > 8))$Mark
}) %>% group_by(Mark) %>% 
  summarize(bd_prop = length(which(Swab == TRUE)) / n()) %>% {
    ggplot(., aes(x = bd_prop)) +
    geom_histogram(bins = 100) +
      xlab("Proportion of samples with bd swabs")
  }

## What does Bd positivity look like by location?
Bd_Newts_AllSites %>% group_by(SA, month, year, Site) %>%
  filter(Swab == TRUE) %>% 
  summarize(
    bd_positive = length(which(Status == "Positive")) / n()
  , num_samples = n()
  ) %>% mutate(Date = as.Date(paste(year, month, "01", sep = "-"))) %>% {
      ggplot(., aes(Date, bd_positive)) + 
      geom_point(aes(size = num_samples, colour = Site)) +
      geom_line(aes(colour = Site)) +
      scale_size_continuous(name = "Samples") +
      ylab("Fraction Bd Positive") +
      xlab("Date") +
      facet_wrap(~SA) +
      theme(axis.text.x = element_text(size = 10))
    }

## what about BD load?
Bd_Newts_AllSites %>% group_by(SA, month, year) %>%
  filter(Swab == TRUE) %>% 
  summarize(
    mid = quantile(copies.swab, 0.50, na.rm = T)
  , lwr = quantile(copies.swab, 0.10, na.rm = T)
  , upr = quantile(copies.swab, 0.90, na.rm = T)
  , bd_positive = length(which(Status == "Positive")) / n()
  , num_samples = n()
  ) %>% mutate(Date = as.Date(paste(year, month, "01", sep = "-"))) %>% {
      ggplot(., aes(Date, mid)) + 
      geom_point(aes(size = num_samples, colour = bd_positive)) +
      geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
      scale_size_continuous(name = "Samples") +
      scale_color_gradient(low = "dodgerblue3", high = "firebrick3"
        , name = "Bd Positive") +
      ylab("Bd 'load'") +
      xlab("Date") +
      facet_wrap(~SA) +
      theme(axis.text.x = element_text(size = 10)) +
      scale_y_log10()
  }

## Pretty strong positive correlation between fraction positive and average load
 ## BUT this average takes into account 0s (non-infected) so it _should_ be the case
  ## that these are positively correlated
Bd_Newts_AllSites %>% group_by(SA, month, year) %>%
  filter(Swab == TRUE) %>% 
  summarize(
    mid = quantile(copies.swab, 0.50, na.rm = T)
  , lwr = quantile(copies.swab, 0.10, na.rm = T)
  , upr = quantile(copies.swab, 0.90, na.rm = T)
  , bd_positive = length(which(Status == "Positive")) / n()
  , num_samples = n()
  ) %>% mutate(Date = as.Date(paste(year, month, "01", sep = "-"))) %>% {
      ggplot(., aes(bd_positive, mid)) + 
      geom_point(aes(size = num_samples)) +
      geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
      scale_size_continuous(name = "Samples") +
      ylab("Bd 'load'") +
      xlab("Bd Positive Fraction") +
      facet_wrap(~SA) +
      theme(axis.text.x = element_text(size = 10)) +
      scale_y_log10()
  }

## A slightly different question is whether they remain positively
 ## correlated subsetting to infected individuals (e.g., because of progression 
  ## at the level of the population and progression within individuals over time)
Bd_Newts_AllSites %>% group_by(SA, month, year) %>%
  filter(Swab == TRUE, Status == "Positive") %>% 
  summarize(
    mid = quantile(copies.swab, 0.50, na.rm = T)
  , lwr = quantile(copies.swab, 0.10, na.rm = T)
  , upr = quantile(copies.swab, 0.90, na.rm = T)
  ) %>% left_join(.
    , 
    {
Bd_Newts_AllSites %>% group_by(SA, month, year) %>%
  filter(Swab == TRUE) %>% 
  summarize(
    bd_positive = length(which(Status == "Positive")) / n()
  , num_samples = n()
  )
    }) %>% mutate(Date = as.Date(paste(year, month, "01", sep = "-"))) %>% {
      ggplot(., aes(bd_positive, mid)) + 
      geom_point(aes(size = num_samples)) +
      geom_errorbar(aes(ymin = lwr, ymax = upr)) + 
      scale_size_continuous(name = "Samples") +
      ylab("Bd 'load'") +
      xlab("Bd Positive Fraction") +
      facet_wrap(~SA) +
      theme(axis.text.x = element_text(size = 10)) +
      scale_y_log10()
  }

## Just a few data issues given that copies/swab can be non-zero for
 ## Status == Negative entries. Not many, but need to be dealt with
Bd_Newts_AllSites %>% filter(Status == "Negative") %>% 
  arrange(desc(copies.swab)) %>% dplyr::select(copies.swab) %>% head(10)

## And within a year?
year.gg1 <- Bd_Newts_AllSites %>% group_by(SA, month, year, Site) %>%
  filter(Swab == TRUE) %>% 
  summarize(
    mid = quantile(copies.swab, 0.50, na.rm = T)
  , lwr = quantile(copies.swab, 0.10, na.rm = T)
  , upr = quantile(copies.swab, 0.90, na.rm = T)
  , bd_positive = length(which(Status == "Positive")) / n()
  , num_samples = n()
  ) %>% mutate(Date = as.Date(paste(year, month, "01", sep = "-"))) %>% {
      ggplot(., aes(month, bd_positive)) + 
      geom_point(aes(size = num_samples, colour = SA, group = Site), alpha = 0.6) +
      geom_line(aes(colour = SA, group = Site), alpha = 0.6) +
      scale_colour_brewer(palette = "Dark2") +
      scale_size_continuous(name = "Samples") +
      ylab("Bd Fraction Positive") +
      theme(
        axis.text.x = element_blank()
      , axis.title.x = element_blank()) +
      xlab("Month") +
      facet_wrap(~year) +
      scale_x_continuous(breaks = c(3, 4, 5, 6, 7), lim = c(3, 7))
  }

year.gg2 <- Bd_Newts_AllSites %>% group_by(SA, month, year, Site) %>%
  filter(Swab == TRUE, Status == "Positive") %>% 
  summarize(
    mid = quantile(copies.swab, 0.50, na.rm = T)
  , lwr = quantile(copies.swab, 0.10, na.rm = T)
  , upr = quantile(copies.swab, 0.90, na.rm = T)
  , bd_positive = length(which(Status == "Positive")) / n()
  , num_samples = n()
  ) %>% mutate(Date = as.Date(paste(year, month, "01", sep = "-"))) %>% {
      ggplot(., aes(month, mid)) + 
      geom_point(aes(size = num_samples, colour = SA, group = Site), alpha = 0.6) +
      geom_line(aes(colour = SA, group = Site), alpha = 0.6) +
      scale_colour_brewer(palette = "Dark2") +
      scale_size_continuous(name = "Samples") +
      ylab("Bd 'load'") +
      xlab("Month") +
      facet_wrap(~year) +
      theme(axis.text.x = element_text(size = 10)) +
      scale_y_log10(
        breaks = c(1E2, 1E3, 1E4, 1E5), 
        labels = c("1E2", "1E3", "1E4", "1E5")) +
      scale_x_continuous(breaks = c(3, 4, 5, 6, 7), lim = c(3, 7))
  }

## It is pretty convincing that Bd load is increasing through the season
 ## but to model this it will be important to separate out if this is because
  ## more individuals are becoming infected or if each individual is experiencing
   ## an increasing load throughout the season
gridExtra::grid.arrange(year.gg1, year.gg2, ncol = 1)

## So... the amount of bd at the individual level?

## First, really well captured individuals. It is pretty hard to look at anything
 ## but faceted by individual, but this makes it difficult to see general patterns
  ## over days
Bd_Newts_AllSites %>% filter(Mark %in% {
  (cap_dist %>% filter(num_swabs > 2))$Mark
}, Site == "P1", Swab == TRUE) %>% {
  ggplot(., aes(julian, copies.swab)) + 
    geom_line(aes(colour = as.factor(year))) +
    geom_point(aes(colour = as.factor(year))) +
    xlab("Date") + 
    ylab("Bd Load") + 
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    scale_y_continuous(trans = "pseudo_log"
      , breaks = c(1, 10, 1E3, 1E4, 1E5, 1E6)) +
    theme(axis.text.x = element_text(size = 10)) +
    facet_wrap(~Mark) +
    theme(axis.text.y = element_text(size = 8))
}

## Don't facet by individual to try and get a better look at temporal patterns
Bd_Newts_AllSites %>% filter(
  Site == "P1"
  , Swab == TRUE) %>% {
  ggplot(., aes(julian, copies.swab)) + 
    geom_line(aes(colour = as.factor(year), group = Mark)) +
    geom_point(aes(colour = as.factor(year), group = Mark)) +
    xlab("Date") + 
    ylab("Bd Load") + 
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    scale_y_continuous(trans = "pseudo_log"
      , breaks = c(1, 10, 1E3, 1E4, 1E5, 1E6)) +
      facet_wrap(~year)
}

## Or all individuals captured that have marks -- 
 ## note that the best bd sampling is in WI which doesn't have many listed
Bd_Newts_AllSites %>% filter(!is.na(year), copies.swab < 2E6) %>% {
  ggplot(., aes(julian, copies.swab)) + 
    geom_line(aes(group = Mark, colour = SA), alpha = 0.5) + 
    xlab("Day") + 
    ylab("Bd Load") + 
    facet_wrap(~year) +
    scale_y_continuous(trans = "pseudo_log"
      , breaks = c(1, 10, 1E3, 1E4, 1E5, 1E6, 2E6)) +
    theme(axis.text.x = element_text(size = 10)) +
    scale_colour_brewer(palette = "Dark2", name = "Year")
}

## Lots of these individuals seem to be swabbed but have no load. But to model
 ## these data need to get to the bottom of new infections vs changes in load
Bd_Newts_AllSites %>% filter(!is.na(year), copies.swab < 2E6
  , Status == "Positive") %>% {
  ggplot(., aes(julian, copies.swab)) + 
    geom_line(aes(group = Mark, colour = SA), alpha = 0.5) + 
    xlab("Day") + 
    ylab("Bd Load") + 
    facet_wrap(~year) +
    scale_y_continuous(trans = "pseudo_log"
      , breaks = c(1, 10, 1E3, 1E4, 1E5, 1E6, 2E6)) +
    theme(axis.text.x = element_text(size = 10)) +
    scale_colour_brewer(palette = "Dark2", name = "Year")
}

## A bit hard to look at due to the sampling and the lack of recaptures; need some violins.
 ## Still a bit deceiving given the variation by site (e.g. see WI)
Bd_Newts_AllSites %>% filter(Status == "Positive") %>% 
  mutate(year_state = interaction(year, SA)) %>%
  mutate(julian = as.factor(julian)) %>% {
  ggplot(., aes(julian, copies.swab)) + 
    geom_violin() + 
    xlab("Day") + 
    ylab("Bd Load") + 
    facet_wrap(~year_state) +
    scale_y_continuous(trans = "pseudo_log", breaks = c(1, 10, 1E3, 1E4, 1E5, 1E6, 2E6)) +
    theme(
      axis.text.x = element_text(size = 8, angle = 300)) +
    ggtitle("Only infected individuals")
  }

## looks like 2020 is the best sampled, just look at that for some better resolution
## Quite noisy. I guess there is some reasonable curvature here[?] with
 ## a general peak around day 150?
## It would be nice to sample later in the season, for example to see if
 ## load in SPR1 falls
Bd_Newts_AllSites %>% filter(!is.na(Mark)) %>% 
  filter(year == 2020) %>% {
  ggplot(., aes(julian, copies.swab)) + 
    geom_line(aes(group = Mark)) + 
    xlab("Day") + 
    ylab("Bd Load") + 
    facet_wrap(~Site) +
    scale_y_continuous(trans = "pseudo_log"
      , breaks = c(1, 10, 1E3, 1E4, 1E5, 1E6, 2E6)) +
    theme(axis.text.x = element_text(size = 10)) 
}

## And by site for just 2020 using violins
Bd_Newts_AllSites %>% filter(!is.na(Mark)) %>% 
  filter(year == 2020) %>% mutate(julian = as.factor(julian)) %>% {
  ggplot(., aes(julian, copies.swab)) + 
    geom_violin() +
    xlab("Day") + 
    ylab("Bd Load") + 
    facet_wrap(~Site) +
    scale_y_continuous(trans = "pseudo_log"
      , breaks = c(1, 10, 1E3, 1E4, 1E5, 1E6, 2E6)) +
    theme(axis.text.x = element_text(size = 8)) 
}

## And with violins. Less clear than with all individuals, but there
 ## does still seem to be some quadratic like pattern
Bd_Newts_AllSites %>% filter(!is.na(Mark)) %>% 
  filter(year == 2020, Status == "Positive") %>% mutate(julian = as.factor(julian)) %>% {
  ggplot(., aes(julian, copies.swab)) + 
    geom_violin() +
    xlab("Day") + 
    ylab("Bd Load") + 
    facet_wrap(~Site) +
    scale_y_continuous(trans = "pseudo_log"
      , breaks = c(1, 10, 1E3, 1E4, 1E5, 1E6, 2E6)) +
    theme(axis.text.x = element_text(size = 8)) 
}

## It looks to me like it may be hard to model individuals' loads as unique
 ## trajectories instead of deviates from some population average over time 
  ## (which could follow a quadratic over the season?)

## What about just positivity?
Bd_Newts_AllSites %>% filter(!is.na(Mark)) %>%
  filter(Mark %in% {
  (cap_dist %>% filter(num_swabs > 3))$Mark
}) %>% 
  filter(Swab == TRUE) %>%
  filter(Status == "Negative" | Status == "Positive") %>% 
  mutate(bd_positive = as.numeric(as.factor(Status)) - 1) %>% {
  ggplot(., aes(julian, bd_positive)) + 
    geom_point(aes(group = Mark)) + 
    xlab("Date") + 
    ylab("Bd Load") + 
    facet_wrap(~Mark) +
    scale_y_continuous(breaks = c(0, 1), lim = c(-0.5, 1.5)) +
    theme(
      axis.text.x  = element_text(size = 12)
    , strip.text.x = element_text(size = 12)
    , axis.text.y  = element_text(size = 12))
  }

## Take a look at a few specific individuals to get a sense of their loads
Bd_Newts_AllSites %>% filter(Mark %in% 
    c("OOYRBB", "RBBBRR", "OxOxxx", "YORRxx", "BxRxRR", "OOxBRR")) %>% 
  mutate(year_mark = interaction(Mark, year)) %>% {
  ggplot(., aes(julian, copies.swab)) + 
    geom_point() +
    scale_y_continuous(trans = "pseudo_log"
      , breaks = c(1, 10, 1E3, 1E4, 1E5, 1E6, 2E6)) +
    xlab("Date") +
    ylab("bd 'load'") +
    facet_wrap(~year_mark)
}

Bd_Newts_AllSites %>% filter(!is.na(Mark)) %>%
  filter(Mark %in% {
  (cap_dist %>% filter(num_swabs > 4))$Mark
}) %>% {
  ggplot(., aes(julian, copies.swab)) + 
    geom_point(aes(colour = as.factor(year))) +
    scale_y_continuous(trans = "pseudo_log"
      , breaks = c(1, 10, 1E3, 1E4, 1E5, 1E6, 2E6)) +
    scale_color_brewer(palette = "Dark2", name = "Year") +
    xlab("Date") +
    ylab("bd 'load'") +
    facet_wrap(~Mark)
}

## given that that is hard to look at, lets facet by individual for just 2020.
 ## Still really hard to look at and very noisy...
  ## One thing that jumps out is that load is clearly measured with
   ## lots of uncertainty -- see ROXXXR for example
Bd_Newts_AllSites %>% filter(!is.na(Mark)) %>%
  filter(Mark %in% {
  (cap_dist %>% filter(num_swabs > 3))$Mark
}) %>% filter(year == 2020) %>% {
  ggplot(., aes(julian, copies.swab)) + 
    geom_point(aes(group = Mark)) +
    geom_line(aes(group = Mark)) + 
    xlab("Date") + 
    ylab("Bd Load") + 
    facet_wrap(~Mark) +
    scale_y_continuous(trans = "pseudo_log"
      , breaks = c(1, 10, 1E3, 1E4, 1E5, 1E6, 2E6)) +
    theme(
      axis.text.x  = element_text(size = 8)
    , strip.text.x = element_text(size = 8)
    , axis.text.y  = element_text(size = 8))
}

## What about a more detailed look at the WI data

## Patterns are very different by Site
 ## Extremely noisy in PA, cleaner in WI (and maybe MA)
gg.j.1 <- Bd_Newts_AllSites %>% 
  filter(SA == "WI") %>% 
  group_by(Date, julian, month, year, Site) %>% 
  summarize(prop_positive = length(which(Status == "Positive")) / n()) %>% ungroup() %>% {
  ggplot(., aes(julian, prop_positive)) +
    geom_point(aes(colour = as.factor(year))) +
    geom_line(aes(colour = as.factor(year))) +
    scale_color_brewer(palette = "Dark2", name = "Year") +
    facet_wrap(~Site) +
      theme(legend.position = c(0.7, 0.7))
}

## Load is pretty consistent
gg.j.2 <- Bd_Newts_AllSites %>% 
  filter(SA == "WI") %>% 
  group_by(Date, julian, month, year, Site) %>% {
  ggplot(., aes(as.factor(julian), copies.swab)) +
    geom_violin() +
    scale_color_brewer(palette = "Dark2", name = "Year") +
      scale_y_log10() +
    facet_wrap(~Site) + 
    xlab("Day") +
    ylab("Bd 'Load'") +
    theme(axis.text.x = element_text(size = 9, angle = 300))
}

## There doesn't seem to be much pattern in Bd load through the year in WI
gridExtra::grid.arrange(gg.j.1, gg.j.2, ncol = 1)

## Lets see an example recapture matrix for a given location (I guess PA for now
 ## as this is the only state with marks). Try for a single site first
Bd_Newts_AllSites %>% 
  filter(SA == "PA") %>% 
  group_by(Site) %>% 
  summarize(num_ind = n_distinct(Mark))

A11 <- Bd_Newts_AllSites %>% 
  filter(
    Site == "A11"
  , year == 2020
  )

capt_history <- expand.grid(
  Date = unique(A11$Date)
, Mark = unique(A11$Mark)
)

## just finding the captures for all unique individuals on each day sampling occurred at A11
capt_history %<>% 
  left_join(.
    , {A11 %>% dplyr::select(Date, Mark, julian, copies.swab)}
    ) %>% rename(
      captured = julian
    , bd_load  = copies.swab) %>% 
  mutate(
    captured = ifelse(is.na(captured), 0, 1)
  , swabbed  = ifelse(is.na(bd_load), 0, 1)) %>%
  mutate(log_bd_load = log(bd_load + 1))  ### eeek!

## On the actual dates
capt_history %>% {
  ggplot(., aes(Date, Mark, fill = as.factor(captured))) + 
    geom_tile(alpha = 0.8) +
    xlab("Sampling Event") + ylab("Individual") +
    scale_fill_manual(
        values = c("dodgerblue4", "firebrick4")
      , name   = "Detected?"
      , labels = c("No", "Yes")) +
    theme(
      axis.text.y = element_text(size = 6)
    , axis.text.x = element_text(size = 6)
    , legend.text = element_text(size = 12)
    , legend.key.size = unit(.55, "cm")
    ) 
}

## collapsed
capt_history %>% mutate(event = as.factor(Date)) %>% {
  ggplot(., aes(event, Mark, fill = as.factor(captured))) + 
    geom_tile(alpha = 0.8) +
    xlab("Sampling Event") + ylab("Individual") +
    scale_fill_manual(
        values = c("dodgerblue4", "firebrick4")
      , name   = "Detected?"
      , labels = c("No", "Yes")) +
    theme(
      axis.text.y = element_text(size = 6)
    , axis.text.x = element_text(size = 8, angle = 300)
    , legend.text = element_text(size = 12)
    , legend.key.size = unit(.55, "cm")
    ) 
}
