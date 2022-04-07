######################################################
## A place to store code for all manuscript figures ##
######################################################

## April 5:
 ## For now just a place to drop some initial supp stuff while working on some writing


#### Supplemental Newt Quadratic Figure ----

Bd_Newts_AllSites   <- read.csv("data/cleaned_cmr_csv/SB_NOVI.csv")

## Stupid dates in R
if (length(grep("/", Bd_Newts_AllSites$CaptureDate[1])) > 0) {

date_convert <- apply(matrix(Bd_Newts_AllSites$CaptureDate), 1, FUN = function (x) {
  a <- strsplit(x, split = "/")[[1]]
  b <- a[3]
  b <- strsplit(b, "")[[1]][c(3, 4)] %>% paste(collapse = "")
  paste(c(a[c(1, 2)], b), collapse = "/")
})

Bd_Newts_AllSites$CaptureDate <- date_convert
Bd_Newts_AllSites      %<>% mutate(CaptureDate = as.Date(CaptureDate, "%m/%d/%y"))

} else {
  
Bd_Newts_AllSites      %<>% mutate(CaptureDate = as.Date(CaptureDate))
  
}

Bd_Newts_AllSites$Year   <- apply(matrix(as.character(Bd_Newts_AllSites$CaptureDate)), 1, FUN = function (x) strsplit(x, "-")[[1]][1]) %>% as.numeric()
Bd_Newts_AllSites$Month  <- apply(matrix(as.character(Bd_Newts_AllSites$CaptureDate)), 1, FUN = function (x) strsplit(x, "-")[[1]][2]) %>% as.numeric()
Bd_Newts_AllSites$Julian <- as.POSIXlt(Bd_Newts_AllSites$CaptureDate)$yday

ggplot(
  Bd_Newts_AllSites %>% filter(BdSample == "Y") %>% droplevels()
  , aes(Julian, TargetCopies.swab)) + 
  geom_line(aes(group = IndividualID, colour = as.factor(Year))) +
  facet_wrap(~Year) +
  scale_y_log10() +
  scale_colour_brewer(palette = "Dark2", name = "Year") +
  xlab("Julian Day") +
  ylab("Bd Load (Copies/Swab)") + 
  theme(legend.key.size = unit(0.8, "cm"))
