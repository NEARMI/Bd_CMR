mehg_range <- data.all %>% 
  group_by(Site) %>%
  summarize(
    lwr   = quantile(MeHg_conc_ppb, 0.025, na.rm = T)
  , lwr_n = quantile(MeHg_conc_ppb, 0.200, na.rm = T)
  , mid   = quantile(MeHg_conc_ppb, 0.500, na.rm = T)
  , upr_n = quantile(MeHg_conc_ppb, 0.800, na.rm = T)
  , upr   = quantile(MeHg_conc_ppb, 0.975, na.rm = T)
  ) %>% arrange(desc(mid)) %>%
  mutate(Site = factor(Site, levels = Site))

data.all %>% {
  ggplot(., aes(x = MeHg_conc_ppb)) + geom_histogram(aes(fill = SubSite), alpha = 0.7) + 
    facet_wrap(~Site, scales = "free")
}

mehg_range %>% {
  ggplot(., aes(y = Site, x = mid)) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0, size = 1.5) +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3, size = 0.5) +
    geom_point(size = 3) +
    xlab("MeHg")
}

capt_history %>% filter(!is.na(merc), !is.na(log_bd_load), !is.na(len), swabbed == 1) %>% droplevels() %>% {
  ggplot(., aes(log_bd_load, merc)) + 
    geom_point(aes(colour = len)) + 
    facet_wrap(~pop_spec, scales = "free") 
}

capt_history %>% 
  filter(!is.na(merc), !is.na(log_bd_load), !is.na(len), swabbed == 1) %>% 
  droplevels() %>% 
  dplyr::select(merc, len, log_bd_load) %>%
  as.matrix() %>% pairs()

