####
## Recover parameters for cluster fits
####

pred_coef        <- as.data.frame(stan.fit.summary[grep("beta"
  , dimnames(stan.fit.summary)[[1]]), c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(
  param     = rownames(.)
, param_set = yeti_set)

stan.ind_pred_eps <- stan.fit.samples$bd_ind_eps %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value)
stan.ind_pred_var <- stan.fit.samples$bd_ind_sigma %>%
  reshape2::melt(.) %>% rename(sd = value) %>%
  left_join(., stan.ind_pred_eps) %>%
  mutate(eps = eps * sd) %>% 
  group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% ungroup()

ind_order.r <- expdat.all %>% 
  group_by(ind) %>% 
  summarize(
   n_swabs = sum(bd_swabbed)
) %>% left_join(
  ., expdat.all %>%
  filter(bd_swabbed == 1) %>%
  group_by(ind) %>% 
  summarize(
   mean_bd = mean(log_bd_load)
)) %>% arrange(desc(mean_bd)) %>% 
  mutate(order_real = seq(n()), ind = as.character(ind)) %>% 
  filter(!is.na(mean_bd))

ind_order.p <- stan.ind_pred_var %>% 
  arrange(desc(mid)) %>% 
  mutate(order_pred = seq(n()), ind = as.character(ind))

ind_order <- left_join(ind_order.r, ind_order.p) %>% mutate(param_set = yeti_set)

if (yeti_set == 1) {
  pred_coef.f <- pred_coef
  ind_order.f <- ind_order
} else {
  pred_coef.f <- rbind(pred_coef.f, pred_coef)
  ind_order.f <- rbind(ind_order.f, ind_order)
}

print(paste("Through", yeti_set, "Parameter Sets", sep = " "))

