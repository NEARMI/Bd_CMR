####
## CMR model Diagnostics
####

stan.fit.summary[grep("beta", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>% {
    ggplot(., aes(Var1, mid)) + 
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3) +
      geom_point() +
      geom_hline(yintercept = 0) +
      xlab("Parameter") +
      ylab("Estimate") +
      theme(axis.text.x = element_text(angle = 300, hjust = 0))
      
  }

## Check that phi is estimated when it should be
capt_history.phi %<>% mutate(
  pred_phi = colMeans(stan.fit.samples$phi)
)

capt_history.phi %>%
  group_by(offseason, time_gaps, phi_ones) %>% 
  summarize(
    mean_surv = mean(pred_phi)
  )

capt_history.phi %>% filter(pred_phi != 0) %>% {
  ggplot(., aes(phi_ones, pred_phi)) + 
    geom_jitter(width = 0.1) +
    ggtitle("Should see a mix of values at phi_ones = 0 and only 1s at phi_ones = 1")
}

## Also check that detection is predicted when it should be, that the gamma index causes p to be low prior to the individual
 ## being detected for the first time, and that chi is predicted when it should be
capt_history.p %<>% mutate(
  pred_p   = colMeans(stan.fit.samples$p)
, pred_chi = colMeans(stan.fit.samples$chi)
)

capt_history.p %>% {
  ggplot(., aes(p_zeros, pred_p)) + 
    geom_jitter(width = 0.1) +
    ggtitle("Should see pretty low values at p_zeros = 0 and overall higher values at p_zeros = 1")
}

## Also check that estimated bd reasonably corresponds to raw values (as a loose check of correct bd indexing)
capt_history.bd_load %<>% mutate(bd_pred = colMeans(stan.fit.samples$X)[capt_history.bd_load$X_stat_index])

capt_history.bd_load %>% {
  ggplot(., aes(log_bd_load, bd_pred)) + geom_point(size = 2) + 
    xlab("Measured bd") + 
    ylab("Predicted bd") +
    geom_abline(slope = 1, intercept = 0)
}

## -- Apparent survival between primary periods within a year as a function of bd and size -- ##

outval   <- seq(1, 10, by = 1)
out.pred <- matrix(nrow = stan.length, ncol = length(outval))

for (i in 1:ncol(out.pred)) {
  out.pred[, i] <- plogis(
    stan.fit.samples$beta_phi[, 1] + 
    stan.fit.samples$beta_phi[, 2] * outval[i] +
    stan.fit.samples$beta_phi[, 3] * 0 # + stan.fit.samples$beta_phi[, 4] * 0
    )
}

out.pred.bd <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "bd")

outval   <- seq(-2, 2, by = .1)
out.pred <- matrix(nrow = stan.length, ncol = length(outval))

for (i in 1:ncol(out.pred)) {
  out.pred[, i] <- plogis(
    stan.fit.samples$beta_phi[, 1] + 
    stan.fit.samples$beta_phi[, 2] * 5 +
    stan.fit.samples$beta_phi[, 3] * outval[i] # + stan.fit.samples$beta_phi[, 4] * 0
    )
}

out.pred.s <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "size")

out.pred <- rbind(out.pred.bd, out.pred.s)

out.pred %>%
  group_by(gap, variable) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% {
  ggplot(., aes(gap, mid)) + 
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
      geom_line(lwd = 2) +
      facet_wrap(~variable, scales = "free") +
    xlab("bd load") +
    ylab("Probability animal remains in population")
  }

## -- survival between years as a function of bd, size, and mercury -- ##

outval   <- matrix(seq(1, 10, by = 1))
out.pred <- matrix(nrow = stan.length, ncol = length(outval))

for (i in 1:ncol(out.pred)) {
   out.pred[, i] <- plogis(
    stan.fit.samples$beta_offseason_year[, 2] + 
    stan.fit.samples$beta_offseason[, 1] * outval[i] +
    stan.fit.samples$beta_offseason[, 2] * 0 + # + stan.fit.samples$beta_offseason[, 4] * 0
    stan.fit.samples$offseason_pop[, 3]
    )
}

out.pred.bd <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "bd")

outval   <- matrix(seq(-2, 2, by = .1))
out.pred <- matrix(nrow = stan.length, ncol = length(outval))

for (i in 1:ncol(out.pred)) {
   out.pred[, i] <- plogis(
    stan.fit.samples$beta_offseason_year[, 2] + 
    stan.fit.samples$beta_offseason[, 1] * 5 +
    stan.fit.samples$beta_offseason[, 2] * outval[i]#  + stan.fit.samples$beta_offseason[, 4] * 0
    )
}

out.pred.s <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "size")

outval   <- matrix(seq(-2, 2, by = .1))
out.pred <- matrix(nrow = stan.length, ncol = length(outval))

for (i in 1:ncol(out.pred)) {
   out.pred[, i] <- plogis(
    stan.fit.samples$beta_offseason[, 1] + 
    stan.fit.samples$beta_offseason[, 2] * 5 +
    stan.fit.samples$beta_offseason[, 3] * 0 +
    stan.fit.samples$beta_offseason[, 4] * outval[i] 
    )
}

out.pred.m <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "Methylmercury")

out.pred <- rbind(out.pred.bd, out.pred.s, out.pred.m)

out.pred %>%
  group_by(gap, variable) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% {
  ggplot(., aes(gap, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line(lwd = 2) +
    facet_wrap(~variable, scales = "free") +
    theme(panel.spacing = unit(2, "lines")) +
    xlab("") +
    ylab("Probability animal remains in population")
  }

## -- detection as a function of bd, mercury, size -- ##

outval   <- matrix(seq(1, 10, by = 1))
out.pred <- matrix(nrow = stan.length, ncol = length(outval))

for (i in 1:ncol(out.pred)) {
   out.pred[, i] <- plogis(
    stan.fit.samples$beta_p[, 1] + 
    stan.fit.samples$beta_p[, 2] * outval[i] +
    stan.fit.samples$beta_p[, 3] * 0 # + stan.fit.samples$beta_p[, 4] * 0
    )
}

out.pred.bd <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "bd")

outval   <- matrix(seq(-2, 2, by = .1))
out.pred <- matrix(nrow = stan.length, ncol = length(outval))

for (i in 1:ncol(out.pred)) {
   out.pred[, i] <- plogis(
    stan.fit.samples$beta_p[, 1] + 
    stan.fit.samples$beta_p[, 2] * 5 +
    stan.fit.samples$beta_p[, 3] * outval[i] # + stan.fit.samples$beta_p[, 4] * 0
    )
}

out.pred.s <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "size")

outval   <- matrix(seq(-2, 2, by = .1))
out.pred <- matrix(nrow = stan.length, ncol = length(outval))

for (i in 1:ncol(out.pred)) {
   out.pred[, i] <- plogis(
    stan.fit.samples$beta_p[, 1] + 
    stan.fit.samples$beta_p[, 2] * 5 +
    stan.fit.samples$beta_p[, 3] * 0 +
    stan.fit.samples$beta_p[, 4] * outval[i] 
    )
}

out.pred.m <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "Methylmercury")

out.pred <- rbind(out.pred.bd, out.pred.s, out.pred.m)

out.pred %>%
  group_by(gap, variable) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% {
  ggplot(., aes(gap, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line(lwd = 2) +
    facet_wrap(~variable, scales = "free") +
    theme(panel.spacing = unit(2, "lines")) +
    xlab("") +
    ylab("Detection Probability")
  }

stan.fit.samples$beta_p_y %>% reshape2::melt() %>% 
  mutate(Year = as.factor(Var2)) %>% 
  mutate(Year = plyr::mapvalues(Year, from = unique(Year), to = c(2018, 2019, 2020, 2021))) %>% {
  ggplot(., aes(x = value, y = (..count..)/sum(..count..))) + 
    geom_histogram(aes(colour = Year, fill = Year), bins = 75, alpha = 0.3) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
      facet_wrap(~Year)
}

stan.fit.samples$beta_p_m %>% reshape2::melt() %>% 
  mutate(Month = as.factor(Var2)) %>% 
  mutate(Month = plyr::mapvalues(Month, from = unique(Month), to = c("June", "July", "Aug"))) %>% {
  ggplot(., aes(x = value, y = (..count..)/sum(..count..))) + 
    geom_histogram(aes(colour = Month, fill = Month), bins = 75, alpha = 0.3) +
    scale_colour_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
      facet_wrap(~Month)
}

## -- estimated individual variation in bd -- ##

stan.ind_pred_var <- stan.fit.samples$X %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value) %>%
  group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% arrange(mid) %>%
  mutate(ind = factor(ind, levels = ind))

stan.ind_pred_var <- cbind(
  Year = capt_history$Year
, ind  = capt_history$Mark
, Rep  = capt_history$SecNumConsec
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

stan.ind_pred_var %>% 
  filter(ind == 45) %>% {
  ggplot(., aes(x = eps)) +
    geom_histogram(bins = 50) +
    facet_wrap(~Year)
  }

capt_history %>% filter(Mark == 87) %>% as.data.frame()

ind_order.r <- capt_history %>% 
  group_by(Mark) %>% 
  summarize(
   n_swabs = sum(swabbed)
) %>% left_join(
  .,capt_history %>%
  filter(swabbed == 1) %>%
  group_by(Mark) %>% 
  summarize(
   tot_bd = mean(log_bd_load)
)) %>% 
  filter(!is.na(tot_bd)) %>%
  arrange(desc(tot_bd)) %>% 
  mutate(order_real = seq(n()), Mark = as.character(Mark))

ind_order.p <- stan.ind_pred_var %>% 
  arrange(desc(mid)) %>% 
  mutate(order_pred = seq(n()), ind = as.character(ind)) %>% 
  rename(Mark = ind)

ind_order <- left_join(ind_order.r, ind_order.p)

ind_order %>% {
  ggplot(., aes(order_real, mid)) + 
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr)) +
    geom_point(aes(order_real, tot_bd), colour = "firebrick3") +
    xlab("Individual") + 
    ylab("Random Effect Deviate") +
    theme(axis.text.x = element_text(size = 8)) +
    geom_hline(yintercept = mean(ind_order$mid)) +
    geom_hline(yintercept = mean(ind_order$tot_bd, na.rm = T), colour = "firebrick3")
}

low_ind  <- unique((capt_history %>% filter(Mark %in% (ind_order %>% arrange(order_real) %>% tail(10))$Mark))$Mark)
high_ind <- unique((capt_history %>% filter(Mark %in% (ind_order %>% arrange(order_real) %>% head(10))$Mark))$Mark)

ind_order$ind_type <- 0
ind_order[ind_order$Mark %in% low_ind, ]$ind_type  <- 1
ind_order[ind_order$Mark %in% high_ind, ]$ind_type <- 2

ggplot(ind_order, aes(order_real, order_pred)) + 
  geom_point(aes(colour = as.factor(ind_type)), size = 2) +
  scale_color_brewer(palette = "Dark2", name = "Individual type", labels = c("Middle", "Lowest Bd", "Highest Bd")) +
  xlab("Individual Ordered by Average of Measured Bd") +
  ylab("Predicted Order")

capt_history %<>% left_join(., )

## -- estimated population variation in all parameters -- ##

stan.fit.summary[grep("offseason_pop", dimnames(stan.fit.summary)[[1]]), ][15:27, ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>% 
  mutate(Var1 = plyr::mapvalues(Var1, from = Var1, to = as.character(u_sites))) %>% {
    ggplot(., aes(Var1, mid)) + 
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3) +
      geom_point() +
      geom_hline(yintercept = 0) +
      xlab("Population") +
      ylab("Estimate") + 
      theme(axis.text.x = element_text(angle = 300, hjust = 0))
  }

stan.fit.summary[grep("bd_ind_sigma", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>% 
  mutate(Var1 = plyr::mapvalues(Var1, from = Var1, to = as.character(u_sites))) %>% {
    ggplot(., aes(Var1, mid)) + 
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3) +
      geom_point() +
      geom_hline(yintercept = 0) +
      xlab("Population") +
      ylab("Estimate") + 
      theme(axis.text.x = element_text(angle = 300, hjust = 0))
  }

stan.fit.summary[grep("inseason_pop", dimnames(stan.fit.summary)[[1]]), ][15:27, ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>% 
  mutate(Var1 = plyr::mapvalues(Var1, from = Var1, to = as.character(u_sites))) %>% {
    ggplot(., aes(Var1, mid)) + 
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3) +
      geom_point() +
      geom_hline(yintercept = 0) +
      xlab("Population") +
      ylab("Estimate") + 
      theme(axis.text.x = element_text(angle = 300, hjust = 0))
  }

stan.fit.summary[grep("bd_pop_year", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>% 
  mutate(Var1 = plyr::mapvalues(Var1, from = Var1, to = seq(n()))) %>% {
    ggplot(., aes(Var1, mid)) + 
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3) +
      geom_point() +
      geom_hline(yintercept = 0) +
      xlab("Population") +
      ylab("Estimate") + 
      theme(axis.text.x = element_text(angle = 300, hjust = 0))
  }
