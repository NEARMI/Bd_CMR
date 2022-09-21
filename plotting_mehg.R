#########################################################
## Plot diagnostics for MeHg/Bd only fit (no survival) ##
#########################################################

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

beta_est <- stan.fit.summary[grep("beta_mehg", dimnames(stan.fit.summary)[[1]]), ] %>%
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%', param = "Var1")

beta_est$param <- c("ANBO:Male", "RANA:Male", "Female", "Unknown", "Length", "Drawdown")

beta_est %<>% mutate(param = factor(param, levels = c(
  rev(c("ANBO:Male", "RANA:Male", "Female", "Unknown", "Drawdown", "Length"))
)))

beta_est %>% {
    ggplot(., aes(mid, param)) + geom_point() +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
      ylab("Parameter") + 
      geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
      theme(
        axis.text.y  = element_text(size = 12)
      , strip.text.x   = element_text(size = 12)) +
      xlab("Coefficient Estimate")
}

beta_est <- stan.fit.summary[grep("beta_bd", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%', param = "Var1")

beta_est$param <- c("ANBO:Male", "RANA:Male", "Female", "Unknown", "Temperature", "Length", "MeHg")

beta_est %<>% mutate(param = factor(param, levels = c(
  rev(c("ANBO:Male", "RANA:Male", "Female", "Unknown", "Temperature", "Length", "MeHg"))
)))

beta_est %>% {
    ggplot(., aes(mid, param)) + geom_point() +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
      ylab("Parameter") + 
      geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
      theme(
        axis.text.y  = element_text(size = 12)
      , strip.text.x   = element_text(size = 12)) +
      xlab("Coefficient Estimate")
}

capt_history %>% {
  ggplot(., aes(len, merc)) + 
      geom_point(aes(colour = Species, shape = Sex)) + 
      facet_wrap(~Site, scales = "free") +
      ylab("MeHg Concentration") +
      xlab("Length") +
      scale_colour_brewer(palette = "Dark2")
}

bd.mehg <- capt_history.bd_load %>% 
  dplyr::select(Mark, pop_spec, Species, bd_load, Year, Sex) %>%
  full_join(.
  , capt_history %>% 
    group_by(Mark, pop_spec, Species, Year, Sex) %>% 
    summarize(merc = mean(merc, na.rm = T))
  ) %>% distinct()

bd.mehg %<>% left_join(
  .
, data.frame(
   Mark      = seq(n_distinct(capt_history$Mark))
 , bd_dev    = stan.fit.summary[grep("bd_ind", dimnames(stan.fit.summary)[[1]]), 1]
 , bd_dev_sd = stan.fit.summary[grep("bd_ind", dimnames(stan.fit.summary)[[1]]), 3]
)
)

bd.mehg$location <- apply(bd.mehg$pop_spec %>% matrix(), 1
  , FUN = function (x) strsplit(x, "[.]")[[1]][2])

bd.mehg %>% {
    ggplot(., aes(merc, bd_load)) + geom_point(aes(colour = Species, shape = Sex)) + 
      facet_wrap(~location, scales = "free") +
      scale_y_continuous(trans = "pseudo_log"
        , breaks = c(0, 10, 100, 1000, 1E4, 1E5, 1E6)) +
      xlab("MeHg Concentration") +
      ylab("Bd Load") +
      scale_colour_brewer(palette = "Dark2")
  }

bd.mehg %>% filter(merc > 30) %>% group_by(Mark) %>% 
  mutate(num_entry = n()) %>% 
filter(
  num_entry > 1
, pop_spec == "RANA.JonesPond"
) %>% {
    ggplot(., aes(merc, bd_load)) + 
    geom_line(aes(group = Mark), alpha = 0.4, size = 0.4) +
    geom_point(aes(colour = as.factor(Year))) + 
      facet_wrap(~pop_spec, scales = "free") +
      scale_y_continuous(trans = "pseudo_log"
        , breaks = c(0, 10, 100, 1000, 1E4, 1E5, 1E6)) +
     scale_x_log10() +
      xlab("MeHg Concentration") +
      ylab("Bd Load") +
    scale_colour_brewer(palette = "Dark2", name = "Year")
}

bd.mehg %>% group_by(Mark) %>% 
  mutate(num_entry = n()) %>% 
filter(
  num_entry > 1
, pop_spec == "RANA.JonesPond"
) %>% {
    ggplot(., aes(bd_load, bd_dev)) + 
    geom_line(aes(group = Mark), alpha = 0.4, size = 0.4) +
    geom_point(aes(colour = as.factor(Year))) + 
      facet_wrap(~pop_spec, scales = "free") +
      scale_x_continuous(trans = "pseudo_log"
        , breaks = c(0, 10, 100, 1000, 1E4, 1E5, 1E6)) +
      ylab("Bd Individual-Level Random Effect Deviate") +
      xlab("Bd Load") +
    scale_colour_brewer(palette = "Dark2", name = "Year")
}


pred_frame <- data.frame(stan.fit.summary[grep("ind_len", dimnames(stan.fit.summary)[[1]]), c(4, 6, 8)]) %>%
  rename(lwr_len = `X2.5.`, mid_len = `X50.`, upr_len = `X97.5.`) %>% 
  cbind(.
, data.frame(stan.fit.summary[grep("ind_mehg", dimnames(stan.fit.summary)[[1]]), c(4, 6, 8)]) %>% 
  rename(lwr_mehg = `X2.5.`, mid_mehg = `X50.`, upr_mehg = `X97.5.`)) %>% 
  cbind(.
, data.frame(stan.fit.summary[grep("bd_ind", dimnames(stan.fit.summary)[[1]]), c(4, 6, 8)]) %>% 
  rename(lwr_bd = `X2.5.`, mid_bd = `X50.`, upr_bd = `X97.5.`)) %>%
  mutate(Mark = seq(n())) %>% 
  left_join(
  ., bd.mehg %>% dplyr::select(Mark, merc) %>% group_by(Mark) %>% slice(1)
) %>% mutate(
  merc = ifelse(is.na(merc), 0, 1)
) %>% left_join(
  ., capt_history %>% dplyr::select(Mark, Sex) %>% group_by(Mark) %>% slice(1)
)

pred_frame  %>% {
  ggplot(., aes(mid_len, mid_mehg)) + 
    geom_errorbar(aes(ymin = lwr_mehg, ymax = upr_mehg), size = 0.4, alpha = 0.4) +
    geom_errorbarh(aes(xmin = lwr_len, xmax = upr_len), size = 0.4, alpha = 0.4) +
    geom_point(aes(colour = as.factor(merc), shape = Sex), alpha = 0.7) +
    scale_colour_manual(values = c("firebrick3", "dodgerblue2"), name = "Missing?", labels = c("Yes", "No")) +
    xlab("Length") + ylab("MeHg")
}

