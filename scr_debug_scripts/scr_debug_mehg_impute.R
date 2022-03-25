data.frame(
  ind.len
, ind.hg
) %>% plot()

stan.fit.summary[grep("beta", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%')

stan.fit.summary_np[grep("beta", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%')

rgamma(1000
  , shape = stan.fit.samples$inverse_phi_mehg[1]
  , rate  = stan.fit.samples$rate_mehg_have[1, 2]
    ) %>% hist()

rgamma(1000
  , shape = stan.fit.samples_np$inverse_phi_mehg[1]
  , rate  = stan.fit.samples_np$rate_mehg_have[1, 2]
    ) %>% hist()

stan.fit.samples$mu_mehg_have[, 3] %>% hist()
stan.fit.samples$rate_mehg_have[, 3] %>% hist()

ind.len[hg.mis[3]]

stan.fit.samples$mu_mehg_mis[, 3] %>% hist()
stan.fit.samples$rate_mehg_mis[, 3] %>% hist()

rgamma(1000
  , shape = stan.fit.samples$inverse_phi_mehg[2]
  , rate  = stan.fit.samples$rate_mehg_mis[1, 2]
    ) %>% hist()

n_cap <- capt_history %>% group_by(capture_date) %>% summarize(num_capt = sum(captured))
n_cap <- n_cap[-1, ]

stan.fit.samples$pop_size %>% 
  reshape2::melt() %>% 
  group_by(Var2) %>% 
  summarize(
    lwr   = quantile(value, 0.025)
  , lwr_n = quantile(value, 0.200)
  , mid   = quantile(value, 0.500)
  , upr_n = quantile(value, 0.800)
  , upr   = quantile(value, 0.975)
  ) %>% rename(Sample_Date = Var2) %>% 
  mutate(Sample_Date = as.character(Sample_Date)) %>%
  mutate(
    Sample_Date = plyr::mapvalues(
      Sample_Date
    , from = Sample_Date
    , to   = as.character(unique(capt_history$capture_date)))) %>%
      mutate(sdate = seq(1, n())) %>% filter(sdate > 1) %>% {
    ggplot(., aes(sdate, mid)) + 
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.5) +
      geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.2) +
      geom_line() +
      scale_x_continuous(
        breaks = seq(1, 12)
      , labels = as.character(unique(capt_history$capture_date))[-1]
        ) +
      geom_point(data = n_cap, aes(capture_date, num_capt), colour = "firebrick3", size = 3) +
      xlab("Date") +
      ylab("Population Estimate") +
      theme(axis.text.x = element_text(angle = 300, hjust = 0))
  }

