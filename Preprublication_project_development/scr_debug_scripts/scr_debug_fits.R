###
## Playing with a few different random effect options (using JONES POND RALU as an example set)
###

stan.fit.summary <- stan.fit.summary_dr_ir_pr

temp_betas4 <- stan.fit.summary[grep("beta", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>%
  mutate(model = "p_random_day_ind_phi")
  
temp_betas <- rbind(temp_betas1, temp_betas2, temp_betas3, temp_betas4)

param_names <- apply(
  matrix((temp_betas %>% filter(model == unique(model)[1]))$Var1 %>% as.character())
, 1
, FUN = function(x) strsplit(x, "[[]")[[1]][1]
)

temp_betas %<>% 
  group_by(model) %>%
  mutate(Var1 = as.character(Var1)) %>% 
  mutate(Var1 = plyr::mapvalues(Var1, from = unique(temp_betas$Var1), to = param_names)) %>%
  rename(params = Var1) %>%
  group_by(params, model) %>%
  mutate(param_lev = seq(n())) %>% 
  relocate(param_lev, .after = params)

temp_betas %>% 
    filter(
      (params == "beta_phi" & param_lev == 1) | 
        (params == "beta_offseason" & param_lev == 2)) %>% {
    ggplot(., aes(model, mid)) +
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3) +
      geom_point() +
      scale_color_brewer(palette = "Dark2", name = "Species") +
      xlab("Model Type") +
      ylab("Estimate") +
      facet_wrap(~params, scales = "free") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme(
        axis.text.x = element_text(angle = 300, hjust = 0)
      )
}

ind_p_est <- stan.fit.samples$p_ind_dev %>% reshape2::melt() %>%
  group_by(Var2) %>%
  summarize(
    lwr_p = quantile(value, 0.025)
  , mid_p = quantile(value, 0.500)
  , upr_p = quantile(value, 0.975)
  ) %>% 
  rename(ind = Var2) %>%
  mutate(ind = factor(ind, levels = ind)) %>% 
  arrange(mid_p) %>% mutate(
   ord_p  = seq(n())
  )

ind_bd_est <- stan.fit.samples$bd_delta_eps %>% reshape2::melt() %>%
  group_by(Var2) %>%
  summarize(
    lwr_bd = quantile(value, 0.025)
  , mid_bd = quantile(value, 0.500)
  , upr_bd = quantile(value, 0.975)
  ) %>% 
  rename(ind = Var2) %>%
  mutate(ind = factor(ind, levels = ind)) %>% 
  arrange(mid_bd) %>% mutate(
   ord_bd  = seq(n())
  )

ind_phi_est <- stan.fit.samples$phi_ind_dev %>% reshape2::melt() %>%
  group_by(Var2) %>%
  summarize(
    lwr_phi = quantile(value, 0.025)
  , mid_phi = quantile(value, 0.500)
  , upr_phi = quantile(value, 0.975)
  ) %>% 
  rename(ind = Var2) %>%
  mutate(ind = factor(ind, levels = ind)) %>% 
  arrange(mid_phi) %>% mutate(
   ord_phi  = seq(n())
  )

ind_est <- left_join(ind_bd_est, ind_p_est) %>% left_join(., ind_phi_est)

ind_p_est %>% mutate(ind = factor(ind, levels = ind)) %>%
    filter(ind %in% sample(seq(1, n_distinct(ind)), 300)) %>% {
  ggplot(., aes(mid_p, ind)) + 
    geom_errorbarh(aes(xmin = lwr_p, xmax = upr_p), height = 0.3) + geom_point()
}

ind_bd_est %>% mutate(ind = factor(ind, levels = ind)) %>%
  filter(ind %in% sample(seq(1, n_distinct(ind)), 300)) %>% {
  ggplot(., aes(mid_bd, ind)) + 
    geom_errorbarh(aes(xmin = lwr_bd, xmax = upr_bd), height = 0.3) + geom_point()
}

ind_est %>% {
  ggplot(., aes(mid_bd, mid_p)) + geom_point()
  }

day_p_est <- stan.fit.samples$p_day_dev %>% reshape2::melt() %>%
  group_by(Var2) %>%
  summarize(
    lwr_p = quantile(value, 0.025)
  , mid_p = quantile(value, 0.500)
  , upr_p = quantile(value, 0.975)
  ) %>% 
  rename(day = Var2) %>%
  mutate(day = factor(day, levels = day)) %>% 
  arrange(mid_p) %>% mutate(
   ord_p  = seq(n())
  )

day_p_est %>% mutate(day = factor(day, levels = day)) %>% {
  ggplot(., aes(mid_p, day)) + geom_errorbarh(aes(xmin = lwr_p, xmax = upr_p), height = 0.3) + geom_point()  
}

