####
## CMR Diagnostics
####

## --- Recovery of simulated coefficients? --- ##

## Primary Bd effects
pred_coef        <- as.data.frame(stan.fit.summary[grep("beta"
  , dimnames(stan.fit.summary)[[1]]), c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(param = rownames(.))

pred_coef %>% mutate(param = factor(param, levels = param)) %>% {
  ggplot(., aes(param, mid)) + 
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    geom_point(data = data.frame(
      param = pred_coef$param
    , mid   = c(
      rep(NA, 4)
    , rev(bd_mort[1, ])
    , NA
    , c(background_mort, p_mort[1])
    , rev(bd_detect[1, ])
      )
    ), colour = "firebrick3") +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
    xlab("Parameter") + 
    ylab("Estimate") +
    scale_x_discrete(labels = c(
      "Bd intercept"
    , "Bd over time"
    , "Bd over time quadratic"
    , "Bd over temp"
    , "Survival intercept"
    , "Survival slope"
    , "Survival timegaps"
    , "Survival offseason background"
    , "Survival offseason bd"
    , "Detection intercept"
    , "Detection slope"
    )) +
    theme(axis.text.x = element_text(size = 11, angle = 300, hjust = 0))
}

## -- survival over Bd load within season -- ##

stan.pred        <- apply(stan.fit.samples$beta_phi, 1
  , FUN = function(x) plogis(x[1] + x[2] * one_pop$bd_probs$log_bd_load)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "mortality")
stan.pred        %<>% mutate(log_bd_load = plyr::mapvalues(log_bd_load
  , from = unique(log_bd_load), to = one_pop$bd_probs$log_bd_load))
stan.pred        %<>% group_by(log_bd_load) %>% 
  summarize(
    lwr = quantile(mortality, c(0.025))
  , mid = quantile(mortality, c(0.5))
  , upr = quantile(mortality, c(0.975))
  )

stan.pred %>% {
  ggplot(., aes(log_bd_load, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line() +
    geom_line(data = (one_pop$bd_probs %>% mutate(mort = mort ))
      , aes(log_bd_load, mort)
      , colour = "dodgerblue4", lwd = 2) +
    xlab("Log of Bd Load") + 
    ylab("Predicted weekly bd mortality")
}

## -- survival over Bd between season -- ##

stan.pred        <- apply(stan.fit.samples$beta_offseason, 1
  , FUN = function(x) plogis(x[1] + x[2] * one_pop$bd_probs$log_bd_load)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "mortality")
stan.pred        %<>% mutate(log_bd_load = plyr::mapvalues(log_bd_load
  , from = unique(log_bd_load), to = one_pop$bd_probs$log_bd_load))
stan.pred        %<>% group_by(log_bd_load) %>% 
  summarize(
    lwr = quantile(mortality, c(0.025))
  , mid = quantile(mortality, c(0.5))
  , upr = quantile(mortality, c(0.975))
  )

stan.pred %>% {
  ggplot(., aes(log_bd_load, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line() +
    geom_line(data = data.frame(
      log_bd_load = stan.pred$log_bd_load
    , mort = plogis(background_mort[1] + p_mort[1] * stan.pred$log_bd_load))
      , aes(log_bd_load, mort)
      , colour = "dodgerblue4", lwd = 2) +
    xlab("Log of Bd Load") + 
    ylab("Predicted mortality between seasons")
}

## -- detection over Bd load -- ##
stan.pred        <- apply(stan.fit.samples$beta_p, 1
  , FUN = function(x) plogis(x[1] + x[2] * one_pop$bd_probs$log_bd_load)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "detect")
stan.pred        %<>% mutate(log_bd_load = plyr::mapvalues(log_bd_load
  , from = unique(log_bd_load), to = one_pop$bd_probs$log_bd_load))
stan.pred        %<>% group_by(log_bd_load) %>% 
  summarize(
    lwr = quantile(detect, c(0.10))
  , mid = quantile(detect, c(0.5))
  , upr = quantile(detect, c(0.90))
  )

stan.pred %>% {
  ggplot(., aes(log_bd_load, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line() +
    geom_line(data = one_pop$bd_probs, aes(log_bd_load, detect)
      , colour = "dodgerblue4", lwd = 2) +
    xlab("Log of Bd Load") + ylab("Predicted detection probability")
}

## --- Individual random effect estimates --- ##

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

stan.ind_pred_var %>% 
  arrange(mid) %>% 
  mutate(ind = factor(ind, levels = ind)) %>% {
  ggplot(., aes(ind, mid)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr)) +
    geom_point() +
    xlab("Individual") + 
    ylab("Random Effect Deviate") +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = 1, colour = "firebrick3") +
    scale_colour_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(size = 8))
}

ind_order.r <- expdat.all %>% 
  group_by(ind) %>% 
  summarize(
   n_swabs = sum(bd_swabbed)
) %>% left_join(
  ., expdat.all %>%
  filter(bd_swabbed == 1) %>%
  group_by(ind) %>% 
  summarize(
   max_bd = max(log_bd_load)
)) %>% arrange(desc(max_bd)) %>% 
  mutate(order_real = seq(n()), ind = as.character(ind)) %>% 
  filter(!is.na(max_bd))

ind_order.p <- stan.ind_pred_var %>% 
  arrange(desc(mid)) %>% 
  mutate(order_pred = seq(n()), ind = as.character(ind))

ind_order <- left_join(ind_order.r, ind_order.p)

ggplot(ind_order, aes(order_real, order_pred)) + 
  geom_point(size = 2) +
  xlab("Individual Ordered by Average of Measured Bd") +
  ylab("Predicted Order")

if (n_pop > 1) {

## --- Population random effect estimates --- ##

stan.ind_pred_eps <- stan.fit.samples$phi_delta_pop_eps %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value)
stan.ind_pred_var <- stan.fit.samples$phi_delta_pop_sigma %>%
  reshape2::melt(.) %>% rename(sd = value) %>%
  left_join(., stan.ind_pred_eps) %>%
  mutate(eps = eps * sd) %>% 
  group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% ungroup()

### Need to make this functioning given the new multi-population structure

stan.ind_pred_var %>% 
  arrange(mid) %>% 
  mutate(ind = factor(ind, levels = ind)) %>% {
  ggplot(., aes(ind, mid)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr)) +
    xlab("Population") + 
    ylab("Random Effect Deviate") +
    geom_hline(yintercept = 0
      , linetype = "dashed", lwd = 1, colour = "firebrick3") +
    scale_colour_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(size = 8))
  }

}

## mean survival per period
ind_occ_phi.all %<>% mutate(pred_phi = colMeans(stan.fit.samples$phi))
ind_occ_phi.all %<>% left_join(.
  , expdat.all %>% filter(sampling_days == 1) %>% rename(sampling_events_phi = all_times) %>%
    dplyr::select(ind, sampling_events_phi, cum_surv)
  )

ind_occ_p.all %<>% 
     mutate(ind_per = interaction(ind, periods_occ)) %>% 
     mutate(ind_per = factor(ind_per, levels = unique(ind_per))) %>% 
     mutate(ind_per = as.numeric(ind_per))

ind_occ_phi.all %>% filter(phi_zeros == 0) %>% 
  dplyr::select(offseason, time_gaps, pred_phi) %>%
  mutate(offseason = as.factor(offseason)) %>% {
  ggplot(., aes(x = pred_phi)) + 
      geom_histogram(bins = 50) +
      facet_wrap(~offseason, scales = "free")
  }

ind_occ_phi.all %>% filter(phi_zeros == 0) %>% 
  dplyr::select(offseason, time_gaps, pred_phi) %>%
  mutate(time_gaps = as.factor(time_gaps)) %>% {
  ggplot(., aes(x = pred_phi)) + 
      geom_histogram(bins = 100) +
      facet_wrap(~time_gaps, scales = "free")
  }
