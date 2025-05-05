ups <- unique(capt_history$pop_spec) 

ups_val <- 22

{
print(ups[ups_val])
  
capt_history %>% 
  filter(pop_spec == ups[ups_val]) %>%
  dplyr::select(size, len, merc) %>%
  distinct() %>%
  mutate(entry = seq()) %>%
  summarize(
    uni_len = length(which(!is.na(len))) / n()
  , uni_siz = length(which(!is.na(size))) / n()
  ) %>% print()
  
capt_history %>% 
  filter(pop_spec == ups[ups_val]) %>%
  dplyr::select(size, len, merc) %>%
  distinct() %>%
  mutate(entry = seq()) %>%
  pivot_longer(-entry, values_to = "meas", names_to = "var") %>% {
    ggplot(., aes(x = meas)) + geom_histogram(bins = 50) + facet_wrap(~var, scales = "free")
  }
}

capt_history %>% 
  filter(pop_spec == ups[ups_val]) %>%
  dplyr::select(size, len, merc) %>%
  distinct() %>% {
    ggplot(., aes(size, len)) + geom_point()
  }

ind.size <- capt_history %>% 
  group_by(Mark, pop_spec) %>% 
  summarize(size = mean(size, na.rm = T)) %>%
  ungroup() %>%
  group_by(pop_spec) %>%
  mutate(size_mean = mean(size, na.rm = T)) %>%
  mutate(size = ifelse(!is.na(size), size, size_mean)) %>%
  dplyr::select(-size_mean) %>%
  mutate(size = scale(size)[, 1])

ind.size <- ind.size$size

ind.size %>% {
  ggplot(., aes(x = size)) + geom_histogram(bins = 50) + facet_wrap(~pop_spec, scales = "free")
}

ind.len <- capt_history %>% 
  group_by(Mark, pop_spec) %>% 
  summarize(len = mean(len, na.rm = T)) %>%
  ungroup() %>%
  group_by(pop_spec) %>%
  mutate(len_mean = mean(len, na.rm = T)) %>%
  mutate(len = ifelse(!is.na(len), len, len_mean)) %>%
  dplyr::select(-len_mean) %>%
  mutate(len = scale(len)[, 1])

ind.len[ind.len$pop_spec == "LilyPond.BCF", ]$len     <- 0
ind.len[ind.len$pop_spec == "MatthewsPond.BCF", ]$len <- 0

ind.len %>% {
  ggplot(., aes(x = len)) + geom_histogram(bins = 50) + facet_wrap(~pop_spec, scales = "free")
}

ind.merc <- capt_history %>% 
  group_by(Mark, pop_spec) %>% 
  summarize(merc = mean(merc, na.rm = T)) %>%
  ungroup() %>%
  group_by(pop_spec) %>%
  mutate(merc_mean = mean(merc, na.rm = T)) %>%
  mutate(merc = ifelse(!is.na(merc), merc, merc_mean)) %>%
  dplyr::select(-merc_mean) %>%
  mutate(merc = scale(merc)[, 1])

ind.merc %>% {
  ggplot(., aes(x = merc)) + geom_histogram(bins = 50) + facet_wrap(~pop_spec, scales = "free")
}


capt_history.p %>% group_by(Site, Species, capture_date) %>% slice(1) %>% 
  filter(Site == "Blackrock-C") %>% 
  as.data.frame() %>% dplyr::select(capture_date)

sampling.effort %>% filter(Site == "Blackrock-C") %>% as.data.frame()



sampling.effort %>% group_by(Site, Species) %>% summarize(n_in_effort = n()) %>% as.data.frame() %>% 
  left_join(. 
    , 
    capt_history.p %>% group_by(Site, Species, capture_date) %>% slice(1) %>% ungroup( ) %>% 
  group_by(Site, Species) %>% summarize(n_in_data = n()) %>% as.data.frame()
    )
