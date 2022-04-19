#########################################
## Clean and combine the raw temp data ##
#########################################

data.files <- list.files("data/raw_temp")
data.locs  <- apply(as.matrix(data.files), 1, FUN = function(x) strsplit(x, ".csv")[[1]])
data.files <- paste("data/raw_temp/", data.files, sep = "")

which_sites <- unique(capt_history$Site) %>% as.character()
which_sites <- which(data.locs %in% which_sites)

data.locs   <- data.locs[which_sites]
data.files  <- data.files[which_sites]

for (i in seq_along(data.files)) {
  
data.temp  <- read.csv(data.files[i]) %>% mutate(
  Site = data.locs[i]
) %>% dplyr::select(year, yday, Site, tmean) %>%
  rename(Year = year) %>%
  filter(Year %in% (capt_history %>% filter(Site == data.locs[i]) %>%
      summarize(Year = unique(Year)))$Year)

data.temp %<>% mutate(Date = as.Date(as.character(paste(Year, 01, 01, sep = "-"))) + yday - 1)

data.month <- apply(matrix(as.character(data.temp$Date)), 1
  , FUN = function(x) strsplit(x, "-")[[1]][2] %>% as.numeric())
 
data.temp$Month <- data.month

## *** Trying a few different possible covariates for temp here.
 ## Merits to different choices, but all have their issues.
data.temp %<>% group_by(Year) %>% mutate(cumtemp = cumsum(tmean))

## Thermal optimum for Bd growth between 17 and 25 C
data.temp %<>% group_by(Year) %>% 
  mutate(in_opt = ifelse(tmean >= 17 & tmean <= 25, 1, 0)) %>%
  mutate(days_in_opt = cumsum(in_opt))


if (i == 1) {
temp_data.all <- data.temp
} else {
temp_data.all <- rbind(temp_data.all, data.temp)
}

}

## *** Need to figure out the correct scaling for temp here...
temp_data.all %<>% ungroup() %>% 
  group_by(Site) %>% mutate(cumtemp_s = scale(cumtemp)[, 1])

## *** Though scaling the number of days in the thermal optima seems like an easy scale
 ## across all years, days, and pops
temp_data.all %<>% ungroup() %>% mutate(days_in_opt_s = scale(days_in_opt)[, 1])

gg.temp_test1 <- temp_data.all %>% {ggplot(., aes(yday, cumtemp_s)) + 
    geom_line(aes(colour = as.factor(Year))) + 
    facet_wrap(~Site, scales = "free")}

gg.temp_test2 <- temp_data.all %>% {ggplot(., aes(yday, days_in_opt_s)) + 
    geom_line(aes(colour = as.factor(Year))) + 
    facet_wrap(~Site, scales = "free")}

