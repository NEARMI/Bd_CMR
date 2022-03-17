##############################################################################
## Print some characteristics of the data for the loop over all populations ##
##############################################################################

data.sum <- capt_history %>% group_by(Mark) %>% filter(captured == 1) %>% 
  summarize(
    mehg = mean(merc, na.rm = T)
  , len  = mean(len, na.rm = T)
  , swab = sum(swabbed, na.rm = T)
    )

print("-------------------------")
print(paste("Dataset has -", nrow(data.sum), "- total captured individuals", sep = " "))
print(paste("Dataset has -", sum(data.sum$swab), "- total -Bd swabs- across all individuals", sep = " "))
print(paste("Dataset has -", length(which(data.sum$swab > 1)), "- total individuals that were -swabbed- at least once", sep = " "))
print(paste("Dataset has -", length(which(!is.na(data.sum$len))), "- total individuals with -length- measrued", sep = " "))
print(paste("Dataset has -", length(which(!is.na(data.sum$mehg))), "- total individuals with -MeHg- measrued", sep = " "))
print("-------------------------")
