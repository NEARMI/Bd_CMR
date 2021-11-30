expdat %>% 
  group_by(ind) %>%
  summarize(num_swab = sum(bd_swabbed)) %>% 
  filter(num_swab == 5)

expdat %>% filter(ind == 27) %>% {
  ggplot(., aes(times, log_bd_load)) + 
    geom_point(lwd = 1, alpha = 0.3) +
    geom_point(
      data = expdat %>% filter(ind == 27, bd_swabbed == 1)
    , lwd = 3, alpha = 1, colour = "firebrick3") +
    facet_wrap(~periods)
}

expdat %>% {
  ggplot(., aes(times, log_bd_load)) + 
    geom_line(aes(group = ind), lwd = 0.5, alpha = 0.3) +
    facet_wrap(~periods)
}
