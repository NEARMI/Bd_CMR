yeti_vals <- expand.grid(
  pop_size = 250
, times    = c(2, 5, 8)
, bd_eff   = c(-0.05, -0.15, -0.25, -0.35)
, ind_var  = c(0, 1, 2)
) %>% mutate(param_set = seq(n()))

yeti_out <- readRDS("CMR_power.rds")

yeti.param_est <- yeti_out[[1]]
yeti.ind_est   <- yeti_out[[2]]

yeti.param_est %>% left_join(., yeti_vals) %>% {
  ggplot(., aes(mid, param_set)) + 
    geom_point() + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr)) + 
    facet_wrap(~param, scales = "free") 
}


yeti.param_est %>% left_join(., yeti_vals) %>% 
  filter(param == "beta_offseason[2]") %>%
  mutate(
    ind_var = as.factor(ind_var)
  , times   = as.factor(times)) %>% {
  ggplot(., aes(mid, param_set)) + 
    geom_point(aes(colour = times, shape = ind_var)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr, colour = times, linetype = ind_var)) + 
    facet_wrap(~bd_eff, scales = "free") +
      geom_vline(aes(xintercept = bd_eff))
  }

yeti.param_est %>% left_join(., yeti_vals) %>% 
  filter(param == "beta_offseason[1]") %>%
  mutate(
    ind_var = as.factor(ind_var)
  , times   = as.factor(times)) %>% {
  ggplot(., aes(mid, param_set)) + 
    geom_point(aes(colour = times, shape = ind_var)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr, colour = times, linetype = ind_var)) + 
    facet_wrap(~bd_eff, scales = "free") +
      geom_vline(xintercept = 2)
  }
