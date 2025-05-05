## bd_beta real vs predicted
reshape2::melt(stan.fit.samples$bd_pop_eps) %>% 
  left_join(., reshape2::melt(stan.fit.samples$bd_pop_sigma) %>% rename(pop_sd = value)) %>%
  left_join(., reshape2::melt(stan.fit.samples$beta_bd) %>% filter(Var2 == 1) %>% 
              dplyr::select(-Var2) %>% rename(mean_beta = value)) %>%
  mutate(value.out = value * pop_sd + mean_beta) %>%
  group_by(Var2) %>%
  summarize(
    lwr = quantile(value.out, 0.025, na.rm = T)
  , mid = quantile(value.out, 0.500, na.rm = T)
  , upr = quantile(value.out, 0.975, na.rm = T)
  ) %>% 
  mutate(real = bd_beta[, 1]) %>%
  mutate(
    samp = lapply(samp, sum) %>% unlist()
  , ind  = all_ind[, 1]) %>% {
    ggplot(., aes(mid, real)) + geom_point(aes(size = samp, colour = ind)) +
      geom_errorbar(aes(ymin = lwr, ymax = upr)) +
      theme_classic() +
      xlab("Estimate") +
      ylab("Simulated")
  }

## beta offseason intercept
hist(stan.fit.samples$beta_offseason[, 1], breaks = 50)

## bd_mort real vs predicted
reshape2::melt(stan.fit.samples$offseason_pop_eps) %>% 
  left_join(., reshape2::melt(stan.fit.samples$offseason_pop_sigma) %>% rename(pop_sd = value)) %>%
  left_join(., reshape2::melt(stan.fit.samples$beta_offseason) %>% filter(Var2 == 2) %>% 
              dplyr::select(-Var2) %>% rename(mean_beta = value)) %>%
  mutate(value.out = value * pop_sd + mean_beta) %>%
  group_by(Var2) %>%
  summarize(
      lwr = quantile(value.out, 0.025, na.rm = T)
    , mid = quantile(value.out, 0.500, na.rm = T)
    , upr = quantile(value.out, 0.975, na.rm = T)
  ) %>% 
  mutate(realv = p_mort[, 1]) %>%
  mutate(
    samp = lapply(samp, sum) %>% unlist()
    , ind  = all_ind[, 1]) %>% {
      ggplot(., aes(Var2, mid)) + 
        geom_point(aes(size = samp, colour = ind)) +
        geom_point(aes(Var2, realv), shape = 2, lwd = 4) +
        geom_errorbar(aes(ymin = lwr, ymax = upr)) +
        theme_classic() +
        xlab("Estimate") +
        ylab("Simulated")
    }

## Other way to look at offseason survival
ind_occ_phi.all %>% 
  mutate(pred = colMeans(stan.fit.samples$phi)) %>% 
  filter(offseason == 1) %>% filter(pred > 0) %>%
  left_join(.
  , expdat.all %>% dplyr::select(ind, log_bd_load, all_times, periods) %>% rename(sampling_events_phi = all_times)) %>% {
    ggplot(.,aes(log_bd_load,pred)) + geom_point() +
      geom_line(
        data = data.frame(
          log_bd_load = seq(21)
        , pred        = plogis(background_mort[1, 1] + p_mort[6, 1] * seq(21))
        )
      , aes(log_bd_load, pred)
      )
  }

## bd_detect real vs predicted
reshape2::melt(stan.fit.samples$p_pop_eps) %>% 
  left_join(., reshape2::melt(stan.fit.samples$p_pop_sigma) %>% rename(pop_sd = value)) %>%
  left_join(., reshape2::melt(stan.fit.samples$beta_p) %>% filter(Var2 == 1) %>% 
              dplyr::select(-Var2) %>% rename(mean_beta = value)) %>%
  mutate(value.out = value * pop_sd + mean_beta) %>%
  group_by(Var2) %>%
  summarize(
    lwr = quantile(value.out, 0.025, na.rm = T)
    , mid = quantile(value.out, 0.500, na.rm = T)
    , upr = quantile(value.out, 0.975, na.rm = T)
  ) %>% 
  mutate(real = bd_detect[, 2]) %>%
  mutate(
    samp = lapply(samp, sum) %>% unlist()
    , ind  = all_ind[, 1]) %>% {
      ggplot(., aes(mid, real)) + geom_point(aes(size = samp, colour = ind)) +
        geom_errorbar(aes(ymin = lwr, ymax = upr)) +
        theme_classic() +
        xlab("Estimate") +
        ylab("Simulated")
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
      xlab("Individual") + 
      ylab("Random Effect Deviate") +
      geom_hline(yintercept = 0, linetype = "dashed", lwd = 1, colour = "firebrick3") +
      scale_colour_brewer(palette = "Dark2") +
      theme(axis.text.x = element_text(size = 8))
  }

ind_order <- expdat.all %>% group_by(ind) %>% summarize(
  tot_bd = max(log_bd_load)
) %>% arrange(desc(tot_bd)) %>% mutate(order_real = seq(n()), ind = as.character(ind)) %>% 
  left_join(.
            , {
              stan.ind_pred_var %>% arrange(desc(mid)) %>% 
                mutate(order_pred = seq(n()), ind = as.character(ind))
            })

ggplot(ind_order, aes(order_real, order_pred)) + geom_point()
