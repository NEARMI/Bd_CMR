data.all %>% filter(pop_spec == "RANA.SanFrancisquito") %>%
  mutate(year = as.factor(Year)) %>% {
  ggplot(., aes(year, bd_load)) + geom_violin() + scale_y_log10()
}

data.all %>% filter(pop_spec == "RANA.SanFrancisquito") %>% 
  mutate(year = as.factor(Year)) %>% {
  ggplot(., aes(julian, bd_load)) + geom_line(aes(colour = year
    , group = interaction(Mark, year))) +
      scale_y_log10()
  }

data.all %>% filter(pop_spec == "RANA.SanFrancisquito") %>% 
  mutate(year = as.factor(Year)) %>% {
  ggplot(., aes(CaptureDate, bd_load)) + geom_line(aes(group = Mark)) +
      scale_y_log10()
  }

check_mark <- (data.all %>% filter(pop_spec == "RANA.SanFrancisquito"
  , BdResult %in% c("pos", "neg")) %>%
  group_by(Mark) %>% summarize(n_swab = n()) %>% 
  arrange(desc(n_swab)) %>% filter(n_swab > 4))$Mark

data.all %>% filter(Mark %in% check_mark) %>% {
  ggplot(., aes(CaptureDate, bd_load)) + 
    geom_line(aes(group = Mark)) +
    geom_point() +
      scale_y_log10() + facet_wrap(~Mark)
  }


nswab <- data.all %>% group_by(pop_spec, Year, Mark) %>% 
  filter(BdResult %in% c("pos", "neg")) %>%
  droplevels() %>%
  summarize(n_swab = n()) 

data.all %>% left_join(., nswab) %>% filter(n_swab > 1) %>% 
  ungroup() %>% group_by(pop_spec, Year) %>%
  droplevels() %>% summarize(n_ind_re = n_distinct(Mark)) %>% as.data.frame()

unique(data.all$pop_spec)
loc <- "NOVI.ScotiaBarrens"
{

san <- data.all %>% 
  filter(pop_spec == loc) %>% 
  filter(BdResult %in% c("pos", "neg")) %>%
  group_by(Mark, Year) %>% 
  mutate(lbd = log(bd_load + 1)) %>%
  summarize(
    mean_bd = exp(mean(lbd))
  ) 

san.f <- expand.grid(
  Mark = unique(san$Mark)
, Year = unique(san$Year)
)

san.f %<>% left_join(., san) %>% 
  mutate(captured = ifelse(is.na(mean_bd), 0, 1)) %>%
  arrange(Mark, Year) %>%
  group_by(Mark) %>%
  mutate(lag_capt = lead(captured, 1))

## Proportion of animals new

cap <- data.all %>% filter(pop_spec == loc) %>%
  group_by(Mark, Year) %>% 
  summarize(
    capt = 1
  )

cap.f <- expand.grid(
  Mark = unique(cap$Mark)
, Year = unique(cap$Year)
)

cap.f %<>% left_join(., cap) %>% 
  group_by(Mark, Year) %>% 
  mutate(capt = ifelse(is.na(capt), 0, 1)) %>% 
  arrange(Mark, Year) %>%
  ungroup() %>%
  group_by(Mark) %>%
  mutate(ccapt = cumsum(capt)) %>%
  mutate(first_capt = Year[min(which(ccapt == 1))]) %>% 
  mutate(cap1 = ifelse(Year == first_capt, 1, 0)) %>%
  ungroup() %>%
  group_by(Year) %>% 
  summarize(
    n_ind = sum(capt)
  , n_new = sum(cap1)
  ) %>% 
  mutate(prop_new = n_new / n_ind) %>%
  left_join(.,
    san %>%
  group_by(Year) %>%
  summarize(bbd = mean(log(mean_bd)))
    ) %>% 
  mutate(pnl = lead(prop_new))

cap.f %>% {
    ggplot(., aes(bbd, pnl)) + geom_point()
  }
}

san.f %>% {
    ggplot(., aes(mean_bd, lag_capt)) + geom_jitter(height = 0.2) +
      facet_wrap(~Year) + geom_smooth() + scale_x_log10()
  }
  
