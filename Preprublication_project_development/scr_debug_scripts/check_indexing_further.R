###
## Another check that the estimates for bd, phi, p, and chi line up with the 
## appropriate days and that they are estimated sensibly
###

capt_history.phi %>% group_by(Mark, Year) %>% slice(1) %>% ungroup() %>% mutate(pred_bd  = stan.fit.samples$X %>% colMeans()) %>%
  filter(Mark == 91) %>%
  dplyr::select(-Month, -Year, -month_year, -Region, -State, -Site, -Species, -merc, -len, -swabbed, -captured) %>% as.data.frame()
  
  capt_history.phi %>% 
    mutate(
  pred_phi = stan.fit.samples$phi %>% colMeans()
    ) %>% filter(Mark == 91) %>%
  dplyr::select(-Month, -Year, -month_year, -Region, -State, -Site, -Species, -merc, -len, -swabbed, -captured) %>% as.data.frame()
  
  capt_history.p %>% 
    mutate(
  pred_p   = stan.fit.samples$p   %>% colMeans()
, pred_chi = stan.fit.samples$chi %>% colMeans()
    ) %>% filter(Mark == 91) %>%
  dplyr::select(-Month, -Year, -month_year, -Region, -State, -Site, -Species, -merc, -len, -swabbed, -captured) %>% as.data.frame()
  