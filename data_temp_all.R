####
## Use daily local temperature 
####

temp.all   <- read.csv("data/final/raw_temp.csv") %>% rename(Year = year) %>% dplyr::select(-dayl, -prcp)
temp_sites <- unique(temp.all$Site)

for (i in seq_along(temp_sites)) {

data.temp <- temp.all %>% filter(Site == temp_sites[i]) 
  
data.temp <- rbind(
  data.temp
,  data.frame(
  Site  = rep(temp_sites[i], unique(data.temp$yday) %>% length())
, Year  = rep(2022, unique(data.temp$yday) %>% length())
, yday  = unique(data.temp$yday)
, tmean = (data.temp %>% group_by(yday) %>% summarize(tmean = mean(tmean)))$tmean
, tmin  = (data.temp %>% group_by(yday) %>% summarize(tmin  = mean(tmin)))$tmin
, tmax  = (data.temp %>% group_by(yday) %>% summarize(tmax  = mean(tmax)))$tmax
)
)

if (i == 1) {
  data.temp.all <- data.temp
} else {
  data.temp.all <- rbind(data.temp.all, data.temp)
}

}

data.temp.all       %<>% mutate(Date = as.Date(as.character(paste(Year, 01, 01, sep = "-"))) + yday - 1)
data.month          <- apply(matrix(as.character(data.temp.all$Date)), 1, FUN = function(x) strsplit(x, "-")[[1]][2] %>% as.numeric())
data.temp.all$Month <- data.month

## *** Trying a few different possible covariates for temp here.
 ## Merits to different choices, but all have their issues.
data.temp.all %<>% group_by(Site, Year) %>% 
  mutate(
    cumtemp_m = cumsum(tmean)
  , cumtemp_l = cumsum(tmin)
  , cumtemp_u = cumsum(tmax)
    )

## Thermal optimum for Bd growth between 17 and 25 C
data.temp.all %<>% group_by(Site, Year) %>% 
  mutate(
    in_opt_m = ifelse(tmean >= 17 & tmean <= 25, 1, 0)
  , in_opt_l = ifelse(tmin >= 17  & tmin <= 25, 1, 0)
  , in_opt_u = ifelse(tmax >= 17  & tmax <= 25, 1, 0)
    ) %>%
  mutate(
    days_in_opt_m = cumsum(in_opt_m)
  , days_in_opt_l = cumsum(in_opt_l)
  , days_in_opt_u = cumsum(in_opt_u)
    )

## *** Need to figure out the correct scaling for temp here...
data.temp.all %<>% ungroup() %>% 
  group_by(Site) %>% 
  mutate(
    cumtemp_m_s = scale(cumtemp_m)[, 1]
  , cumtemp_l_s = scale(cumtemp_l)[, 1]
  , cumtemp_u_s = scale(cumtemp_u)[, 1]
    )

## *** Though scaling the number of days in the thermal optima seems like an easy scale
 ## across all years, days, and pops
data.temp.all %<>% ungroup() %>% 
  mutate(
    days_in_opt_m_s = scale(days_in_opt_m)[, 1]
  , days_in_opt_l_s = scale(days_in_opt_l)[, 1]
  , days_in_opt_u_s = scale(days_in_opt_u)[, 1]
    )

