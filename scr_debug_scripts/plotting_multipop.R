###################################################
## Plot diagnostics for a joint population model ##
###################################################

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

####
## Plotting Setup
####

stan.fit.summary[grep("beta", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>% as.data.frame()

stan.fit.summary[grep("inseason_pop", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%')

stan.fit.summary[grep("offseason_pop", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%')

stan.fit.summary[grep("bd_pop_year", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%')

stan.fit.summary[grep("p_pop", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%')

spec_in_pop <- (capt_history.p %>% group_by(pop_spec) %>% slice(1))$Species %>% as.numeric()

pred.vals <- expand.grid(
  bd   = seq(0, 14, by = 1)
, len  = seq(-3, 3, by = 0.5)
, mehg = seq(-3, 3, by = 0.5)
, pop  = seq(1, dim(stan.fit.samples$offseason_pop_eps)[2], by = 1)
) %>% mutate(
  spec = spec_in_pop[pop]
)

pred.est <- matrix(data = 0, nrow = nrow(pred.vals), ncol = dim(stan.fit.samples[[1]])[1])

for (j in 1:nrow(pred.est)) {
   pred.est[j, ] <- plogis(
    stan.fit.samples$beta_offseason_int[, pred.vals$spec[j]] +
    stan.fit.samples$beta_offseason_sex[, 1] +
    stan.fit.samples$offseason_pop[, pred.vals$pop[j]] + 
    (stan.fit.samples$beta_offseason_bd[, pred.vals$spec[j]] + stan.fit.samples$offseason_pop_bd[, pred.vals$pop[j]]) * pred.vals$bd[j] +
    (stan.fit.samples$beta_offseason_len[, pred.vals$spec[j]] + stan.fit.samples$offseason_pop_len[, pred.vals$pop[j]]) * pred.vals$len[j] +
    stan.fit.samples$beta_offseason_mehg[, pred.vals$spec[j]] * pred.vals$mehg[j]
    )
}

pred.vals <- cbind(pred.vals, pred.est)

pred.vals %<>% pivot_longer(., c(-bd, -len, -mehg, -pop, -spec), names_to = "iter", values_to = "est")

pred.vals.gg <- pred.vals %>%
  group_by(bd, len, mehg, spec, pop) %>%
  summarize(
    lwr   = quantile(est, 0.025)
  , lwr_n = quantile(est, 0.200)
  , mid   = quantile(est, 0.500)
  , upr_n = quantile(est, 0.800)
  , upr   = quantile(est, 0.975)
  )

pred.vals.gg %<>% mutate(pop =  as.factor(pop), spec = as.factor(spec))

pred.vals.gg %>% filter(len == 0, mehg == 0) %>% {
  ggplot(., aes(bd, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = pop), alpha = 0.3) +
    geom_line(aes(colour = pop), size = 1) + 
    scale_colour_discrete() +
    scale_fill_discrete() +
    facet_grid(~spec*pop)
}

pred.vals.gg %>% filter(mehg == 0, bd == 7) %>% {
  ggplot(., aes(len, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = pop), alpha = 0.3) +
    geom_line(aes(colour = pop), size = 1) + 
    scale_colour_discrete() +
    scale_fill_discrete() +
    facet_grid(~spec*pop)
}

beta_est <- stan.fit.summary[grep("beta", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%')

param_names <- apply(
  matrix(beta_est$Var1 %>% as.character())
, 1
, FUN = function(x) strsplit(x, "[[]")[[1]][1]
)

beta_est %<>% 
  mutate(Var1 = as.character(Var1)) %>% 
  mutate(Var1 = plyr::mapvalues(Var1, from = unique(beta_est$Var1), to = param_names)) %>%
  rename(params = Var1) %>%
  group_by(params) %>%
  mutate(param_lev = seq(n())) %>% 
  relocate(param_lev, .after = params) %>%
  mutate(param_lev = as.character(param_lev))


beta_est %>% filter(params %in% c("beta_bd", "beta_inseason", "beta_p_spec", "ind_len_beta", "beta_mehg_spec")) %>% 
  mutate(param_lev = plyr::mapvalues(param_lev, from = c(1,2,3,4), to = c("AMCI", "ANBO", "BCF", "RANA"))) %>% {
    ggplot(., aes(mid, param_lev)) + geom_point() +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3) +
      facet_wrap(~params, scales = "free") +
      ylab("Species")
  }

beta_est %>% filter(params %notin% c("beta_bd", "beta_inseason", "beta_p_spec", "ind_len_beta", "beta_mehg_spec")) %>% 
  mutate(param_lev = plyr::mapvalues(param_lev, from = c(1,2,3,4), to = c("AMCI", "ANBO", "BCF", "RANA"))) %>% {
    ggplot(., aes(mid, param_lev)) + geom_point() +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3) +
      geom_vline(xintercept = 0) +
      facet_wrap(~params, scales = "free") +
      ylab("Species")
  }

stan.ind_pred_var <- stan.fit.samples$X %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value) %>%
  group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% arrange(mid) %>%
  mutate(ind = factor(ind, levels = ind))

capt_history.temp <- capt_history %>% 
  filter(pop_spec == this_pop)
capt_history.slice <- capt_history.temp %>% group_by(X_stat_index) %>% slice(1)

stan.ind_pred_var <- cbind(
  Year = capt_history.slice$Year
, ind  = capt_history.slice$Mark
, Rep  = capt_history.slice$SecNumConsec
, t(stan.fit.samples$X)
) %>% as.data.frame() %>% 
  reshape2::melt(c("Year", "ind", "Rep")) %>% 
  dplyr::select(-Rep) %>%
  rename(iter = variable, eps = value) %>% 
  distinct() %>%
  group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% 
  arrange(mid) %>%
  mutate(ind = factor(ind, levels = ind))

stan.ind_pred_var %<>% mutate(
    population = this_pop
  , location   = this_loc
  , species    = this_spec  
)

ind_order.r <- capt_history.temp %>% 
  group_by(Mark) %>% 
  summarize(
   n_swabs = sum(swabbed)
) %>% left_join(
  ., capt_history.temp %>%
  filter(swabbed == 1) %>%
  group_by(Mark) %>% 
  summarize(
   tot_bd = mean(log_bd_load)
)) %>% mutate(
  given_mean = ifelse(
    is.na(tot_bd)
  , 1
  , 0
  ))

no_swabs <- (ind_order.r %>% filter(given_mean == 0))$Mark

ind_order.r %<>% 
  filter(given_mean == 0) %>%
  arrange(desc(tot_bd)) %>% 
  mutate(order_real = seq(n()), Mark = as.character(Mark))

ind_order.p <- stan.ind_pred_var %>% 
  filter(ind %in% no_swabs) %>%
  arrange(desc(mid)) %>% 
  mutate(order_pred = seq(n()), ind = as.character(ind)) %>% 
  rename(Mark = ind)

ind_order <- left_join(ind_order.r, ind_order.p)

ind_order %<>% mutate(
    population = this_pop
  , location   = this_loc
  , species    = this_spec 
)

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

pop_size_est <- stan.fit.samples$pop_size %>% 
  reshape2::melt() %>% 
  mutate(value = ifelse(value == 0, NA, value)) %>%
  group_by(Var2) %>% 
  summarize(
    lwr   = quantile(value, 0.025, na.rm = T)
  , lwr_n = quantile(value, 0.200, na.rm = T)
  , mid   = quantile(value, 0.500, na.rm = T)
  , upr_n = quantile(value, 0.800, na.rm = T)
  , upr   = quantile(value, 0.975, na.rm = T)
  ) %>% 
  rename(Sample_Date = Var2) %>% 
  mutate(Sample_Date = as.character(Sample_Date)) %>%
  mutate(
    Sample_Date = plyr::mapvalues(
      Sample_Date
    , from = Sample_Date
    , to   = as.character(unique(capt_history$capture_date)))) %>%
      mutate(sdate = seq(1, n())) %>% filter(sdate > 1)

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


