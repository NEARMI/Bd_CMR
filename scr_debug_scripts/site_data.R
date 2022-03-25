capt_history %>% group_by(pop_spec) %>% slice(1) %>% 
  as.data.frame() %>% dplyr::select(Site, Species) %>% arrange(Species) %>% 
  left_join(.,
    capt_history %>% 
  group_by(Site, Species) %>% summarize(
  num_ind   = n_distinct(Mark)
, num_cap   = sum(captured)
, num_swab  = sum(swabbed)
    )
    ) %>% left_join(.
    , 
    capt_history %>% 
  group_by(Site, Species, Mark) %>% summarize(
  tot_swab = sum(swabbed)
    ) %>% ungroup(Mark) %>%
        summarize(prop_swabbed = round(length(which(tot_swab > 0))/n(), 2))
      ) %>% left_join(.,
    capt_history %>% 
  group_by(Site, Species, Mark) %>% summarize(
  tot_merc = ifelse(length(which(!is.na(merc))) > 0, 1, 0)
    ) %>% ungroup(Mark) %>%
        summarize(prop_merc = round(length(which(tot_merc == 1))/n(), 2))    
        ) %>% left_join(.,
          capt_history %>% group_by(Site, Species, Mark) %>% 
  mutate(
    capt_num = cumsum(captured)
  , swab_num = cumsum(swabbed)) %>%
  summarize(
    tot_cap        = max(capt_num)
  , swab_after_cap = ifelse(max(which(captured == 1)) > min(which(swabbed == 1)), 1, 0)
  ) %>% ungroup(Mark) %>% summarize(
    prop_cap_after_swab = (length(which(swab_after_cap == 1)) / n()) %>% round(2)
  , tot_cap        = (length(which(tot_cap > 1)) / n()) %>% round(2)
  ))
