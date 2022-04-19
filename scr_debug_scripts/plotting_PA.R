####
## Alt plotting script for PA -- continuous Bd curves throughout the season
####

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

## check correct indexing of informing X and correspondence of max calculation with X indices

Xfirst      <- (capt_history.p %>% ungroup() %>% mutate(row_index = seq(n())) %>% group_by(X_stat_index) %>% slice(1))$row_index
which_index <- 77

stan.fit.samples$X[50, Xfirst[which_index]:(Xfirst[which_index + 1] - 1)] %>% max()
stan.fit.samples$X_max[50, which_index]

## Add predicted bd and plot that 
capt_history.p %<>% mutate(
  pred_bd = stan.fit.samples$X %>% colMeans()
)

capt_history.p %>% filter(swabbed == 1) %>% {
  ggplot(., aes(yday, log_bd_load)) + geom_point() +
    geom_line(aes(group = Mark), alpha = 0.3) +
    geom_point(aes(yday, pred_bd), alpha = 0.5, colour = "firebrick3") +
    facet_wrap(~Year, ncol = 1)
}

capt_history.p %>% filter(swabbed == 1) %>% {
  ggplot(., aes(log_bd_load, pred_bd)) + geom_point() +
    geom_line(aes(group = Mark), alpha = 0.3) +
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~Year, ncol = 1)
}

## Individual Bd trajectories
stan.ind_pred_var <- stan.fit.samples$X %>% reshape2::melt(.) %>%
  group_by(Var2) %>%
  summarize(
    mid = quantile(value, 0.50)
  , lwr = quantile(value, 0.025)
  , upr = quantile(value, 0.975)
  ) %>% cbind(., (capt_history.p %>% dplyr::select(capture_date, Mark, Year, swabbed, yday_s, log_bd_load)))

stan.ind_pred_var %>% filter(Mark < 20) %>% {
  ggplot(., aes(yday_s, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line(aes(group = Mark)) +
    geom_point(
      data = stan.ind_pred_var %>% filter(Mark < 20, swabbed == 1)
      , aes(yday_s, log_bd_load), colour = "firebrick2") +
    facet_wrap(~Mark)
}




## Betas
beta_est <- stan.fit.summary[grep("beta", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%')

## Between-season survival
pred.vals <- expand.grid(
  bd   = seq(0, 14, by = 1)
, len  = seq(-3, 3, by = 0.5)
)

pred.est <- matrix(data = 0, nrow = nrow(pred.vals), ncol = dim(stan.fit.samples[[1]])[1])

for (j in 1:nrow(pred.est)) {
   pred.est[j, ] <- plogis(
    stan.fit.samples$beta_offseason_sex[, 1] +
    stan.fit.samples$beta_offseason[, 1] * pred.vals[j, ]$bd + 
    stan.fit.samples$beta_offseason[, 2] * pred.vals[j, ]$len
   )
}

pred.vals <- cbind(pred.vals, pred.est)

pred.vals %<>% pivot_longer(., c(-bd, -len), names_to = "iter", values_to = "est")

pred.vals.gg <- pred.vals %>%
  group_by(bd, len) %>%
  summarize(
    lwr   = quantile(est, 0.025)
  , lwr_n = quantile(est, 0.200)
  , mid   = quantile(est, 0.500)
  , upr_n = quantile(est, 0.800)
  , upr   = quantile(est, 0.975)
  )

stan.ind_pred_var <- cbind(
data.frame(
  per  = capt_history.p$X_stat_index
, ind  = capt_history.p$Mark
) %>% distinct()
, t(stan.fit.samples$X_max)
) %>% as.data.frame() %>% 
  reshape2::melt(c("per", "ind")) %>% 
  rename(iter = variable, eps = value) %>% 
  distinct() %>%
  group_by(ind, per) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% group_by(ind) %>% slice(1) %>%
  dplyr::select(-per)

ind_order.r <- capt_history %>% 
  group_by(Mark) %>% 
  summarize(
   n_swabs = sum(swabbed)
) %>% left_join(
  ., capt_history %>%
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
  ungroup( ) %>%
  filter(ind %in% no_swabs) %>%
  arrange(desc(mid)) %>% 
  mutate(order_pred = seq(n()), ind = as.character(ind)) %>% 
  rename(Mark = ind)

ind_order <- left_join(ind_order.r, ind_order.p)

ind_order %>% {
  ggplot(., aes(order_real, order_pred)) + geom_point()
}

