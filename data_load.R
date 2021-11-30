##################################################
## Load each of the data files in the directory ##
##################################################

data.all <- list.files("data")
data.all <- data.all[-which(data.all == "xlsx")]
data.all <- paste("data/", data.all, sep = "")

for (i in seq_along(data.all)) {

  data.temp <- read.csv(data.all[i])
  
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

data.temp$CaptureDate <- date_convert
data.temp      %<>% mutate(CaptureDate = as.Date(CaptureDate, "%m/%d/%y"))

} else {
data.temp      %<>% mutate(CaptureDate = as.Date(CaptureDate))
}
  
####
## A few other modifications
####

## name changes for convenience
data.temp %<>% rename(bd_load = TargetCopies.swab) %>%
  mutate(bd_load = as.numeric(bd_load)) %>%
  rename(Mark = PitTagCode) %>%
  filter(!is.na(Mark))

## check which years have no swabbing and remove them
no.swabyear <- data.temp %>% group_by(Year) %>% 
  summarize(tot_ss = length(which(!is.na(bd_load)))) %>% filter(tot_ss == 0)

data.temp %<>% filter(Year %notin% no.swabyear$Year) %>% droplevels()
  
if (i == 1) {
  data.all <- data.temp
} else {
  data.all <- rbind(data.all, data.temp)
}
  
}

