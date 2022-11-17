capt_history

real_data %>% {
  ggplot(., aes(JD, TargetCopies.swab)) + 
    geom_point(aes(colour = IndividualID)) +
    geom_line(aes(colour = IndividualID)) +
    scale_y_log10() +
    xlab('Julian Day') + ylab('Bd Load') +
    facet_wrap(~IndividualID)
}

capt_history.bd_load %<>% mutate(JD = as.POSIXlt(capture_date)$yday)
unique(capt_history.bd_load$pop_spec)

capt_history.gg <- capt_history.bd_load %>% 
  group_by(Mark, pop_spec, Year) %>% 
  summarize(n_meas = n()) %>% 
  arrange(desc(n_meas)) %>% 
  ungroup() %>% 
  group_by(pop_spec) %>%
 # filter(n_meas >= (max(n_meas, na.rm = T) - 1)) %>%
  slice(c(1:20)) %>%
  arrange(desc(n_meas)) %>%
  ungroup() %>%
  left_join(., capt_history.bd_load) %>% 
  ungroup()

capt_history.gg %>% 
 filter(pop_spec == c(
   "ANBO.TwoMedicine"
   )) %>% {
  ggplot(., aes(JD, bd_load)) + 
   geom_point() +
   geom_line() +
   scale_y_log10() +
   xlab('Julian Day') + 
   ylab('Bd Load') +
   facet_wrap(~Mark)
 }

capt_history.gg %>% 
  group_by(Mark, pop_spec) %>%
  summarize(bd_load = mean(bd_load)) %>%
  ungroup() %>%
  group_by(pop_spec) %>%
  arrange(desc(bd_load)) %>%
  mutate(ind_rank = seq(n())) %>% 
  dplyr::select(-c(bd_load)) %>%
  left_join(., capt_history.gg) %>%
  filter(pop_spec == "NOVI.ScotiaBarrens") %>% arrange(desc(ind_rank)) %>%
  mutate(Mark = factor(Mark, levels = unique(Mark))) %>% {
  ggplot(., aes(Mark, bd_load)) + 
   geom_line(aes(group = Mark)) +
   geom_point(aes(colour = JD), size = 3) +  
   scale_colour_continuous(type = "viridis") +
   scale_y_log10() +
   xlab('Individual') + 
   ylab('Bd Load') +
      theme(axis.text.x = element_text(angle = 300, hjust = 0))
  }

capt_history.gg %>%
  filter(pop_spec %in% c("BCF.MatthewsPond")) %>% {
  ggplot(., aes(JD, bd_load)) + 
   geom_point(aes(group = Mark, colour = pop_spec)) +
   geom_line(aes(group = Mark, colour = pop_spec)) +
   scale_colour_brewer(palette = "Dark2", name = "Population") +
   scale_y_log10() +
   xlab('Julian Day') + 
   ylab('Bd Load') 
  }


real_data %>% left_join(.
  , real_data %>% filter(!is.na(TargetCopies.swab)) %>% 
  group_by(IndividualID, year) %>% 
  summarize(n_meas = n()) %>% arrange(desc(n_meas))) %>%
  ungroup() %>% 
  ## Select just those individuals that were swabbed for Bd many times within a given year
  filter(
  #   n_meas == max(n_meas, na.rm = T)
      n_meas >= (max(n_meas, na.rm = T) - 2)
    ) %>%
  mutate(log_bd = log(TargetCopies.swab)) %>%
  left_join(., real_data %>% group_by(year, JD) %>% 
      summarize(temp = mean(Temp, na.rm = T)) %>% mutate(temp = round(temp, 1)))




####
## Fit example model
####

glm(
  n_year ~ bd_load
  , family = "binomial"
  , data = {
    capt_history.bd_load %>% 
  filter(pop_spec == "NOVI.SPR") %>% 
  mutate(log_bd_load = ifelse(Month == 5, log_bd_load, NA)) %>%
  group_by(Mark) %>%
  summarize(
    n_year = n_distinct(Year)
  , bd_load = mean(log_bd_load)) %>%
  mutate(n_year = ifelse(n_year > 1, 1, 0))
  }
) %>% confint() %>% as.data.frame() %>% 
  mutate(
    param = c("int", "bd_load")
  ) %>% filter(param == "bd_load") %>% {
    ggplot(., aes(param, ymin = `2.5 %`, ymax = `97.5 %`)) + 
      geom_errorbar(width = 0.2) +
      xlab("NOVI.SPR") + theme(axis.text.x = element_blank())
  }

