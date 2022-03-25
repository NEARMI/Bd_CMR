capt_history.temp  <- capt_history %>% filter(pop_spec == this_pop)
capt_history.slice <- capt_history.temp %>% 
  group_by(capture_date) %>% 
  filter(captured == 1) %>% summarize(n_capts = n()) %>%
  arrange(n_capts) %>%
  mutate(rough_real_order = seq(n())) %>% 
  arrange(capture_date)

stan.p_pred_baseline <- stan.fit.samples$beta_p %>% reshape2::melt(.)

stan.p_pred_var <- stan.fit.samples$p_day_dev %>%
  reshape2::melt(.) %>% 
  rename(day = Var2, eps = value) %>% 
  left_join(., stan.p_pred_baseline) %>%
  mutate(pred_p = plogis(eps + value)) %>%
  group_by(day) %>%
  summarize(
    mid = quantile(pred_p, 0.50)
  , lwr = quantile(pred_p, 0.025)
  , upr = quantile(pred_p, 0.975)
  ) %>% arrange(day) %>%
  mutate(day = plyr::mapvalues(day, from = unique(day), to = as.character(unique(capt_history.temp$capture_date))))

stan.p_pred_var %<>% arrange(mid) %>% mutate(
  capture_date = as.Date(day)
, pred_order = seq(n())
  ) %>% dplyr::select(-day) %>% arrange(capture_date)

stan.p_pred_var %<>% left_join(., capt_history.slice)

stan.p_pred_var %<>% mutate(
    population = this_pop
  , location   = this_loc
  , species    = this_spec 
)
