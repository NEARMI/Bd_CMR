##################################################
## Load each of the data files in the directory ##
##################################################

data.files <- list.files("data")
data.files <- data.files[-which(data.files == "xlsx")]
data.files <- paste("data/", data.files, sep = "")

for (i in seq_along(data.files)) {

  data.temp <- read.csv(data.files[i])
  
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
  summarize(tot_ss = length(which(!is.na(bd_load)))) %>% filter(tot_ss == 0)

data.temp %<>% filter(Year %notin% no.swabyear$Year) %>% droplevels()

data.temp %<>% dplyr::select(
    Site, Species, CaptureDate, Year, Month, PrimNum, SecNumConsec
  , Mark, BdSample, SVLmm, MassG, bd_load) %>%
  mutate(dataset = i)

## Different data sets define SecNumConsec in different ways. Homogenize the choice by just making this 
 ## variable a count from 1 to n()

adj.SecNumConsec <- data.temp %>% group_by(CaptureDate) %>% 
  summarize(SecNumConsec = unique(SecNumConsec)) %>% 
  ungroup() %>%
  mutate(SecNumConsec_corrected = seq(n()))

data.temp %<>% left_join(., adj.SecNumConsec) %>% 
  dplyr::select(-SecNumConsec) %>% rename(SecNumConsec = SecNumConsec_corrected)

## drop all entries with no mark
data.temp %<>% filter(Mark != "")

## -- this will have to be non-dynamic?? -- ##
## certain species disappear and capture events become mostly opportunisitc. Need to figure out what to 
 ## do with these data, but for now just drop them
if (i == 1) {
  data.temp %<>% filter(Month < 7)
}

if (i == 1) {
  data.all <- data.temp
} else {
  data.all <- rbind(data.all, data.temp)
}
  
}

