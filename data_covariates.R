##############################################################
## Deal with any other needed covariates for the stan model ##
##############################################################

## -- Individual specific covariates -- ##

ind.size <- (capt_history %>% group_by(Mark) %>%
  summarize(size = mean(size, na.rm = T)))$size
ind.size[which(is.na(ind.size))] <- mean(ind.size[-which(is.na(ind.size))])
ind.size <- scale(ind.size)[, 1]

ind.hg <- (capt_history %>% group_by(Mark) %>%
  summarize(merc = mean(merc, na.rm = T)))$merc
ind.hg[which(is.na(ind.hg))] <- mean(ind.hg[-which(is.na(ind.hg))])
ind.hg <- scale(ind.hg)[, 1]

## -- Population level covariates -- ##

