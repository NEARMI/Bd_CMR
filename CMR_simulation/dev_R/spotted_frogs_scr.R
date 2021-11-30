test_capt <- capt_history %>% 
  filter(Year == 2020) %>% 
  droplevels() %>%
  group_by(Mark, Month) %>%
  summarize(tot_capts = sum(captured)) %>% 
  mutate(tot_capts = ifelse(tot_capts > 0, 1, 0)) %>%
  ungroup(Month)

cc <- test_capt %>% group_by(Mark) %>% summarize(tot_capts = sum(tot_capts)) %>%
 filter(tot_capts > 0) %>% droplevels()

test_capt %<>% filter(Mark %in% cc$Mark)

uni_ind <- unique(test_capt$Mark)
ind_type <- numeric(length(unique(test_capt$Mark)))

for (i in 1:length(ind_type)) {
 t_ind <- test_capt[test_capt$Mark == uni_ind[i], ]
 
 if (t_ind$tot_capts[1] == 1) {
   ind_type[i] <- ifelse(sum(t_ind$tot_capts[2:3]) == 0, 1, 0)
 } else if (t_ind$tot_capts[2] == 1 & t_ind$tot_capts[1] == 0) {
   ind_type[i] <- ifelse(sum(t_ind$tot_capts[3]) == 0, 1, 0)
 } else if (t_ind$tot_capts[3] == 1 & t_ind$tot_capts[1] == 0 & t_ind$tot_capts[2] == 0) {
   ind_type[i] <- 2
 } else {
   ind_type[i] <- 2
 }
 
}

## Percent of individuals in 2021 that were never recaptured
mean(ind_type[ind_type != 2])
