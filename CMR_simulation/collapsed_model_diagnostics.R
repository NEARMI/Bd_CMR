
## -- survival over time -- ##

out.pred <- matrix(nrow = 1000, ncol = 10)
ttg <- seq(1, 10, by = 1)

for (i in 1:ncol(out.pred)) {
  out.pred[, i] <- plogis(stan.fit.samples$beta_phi + stan.fit.samples$beta_timegaps * ttg[i])
}

reshape2::melt(out.pred) %>% rename(iter = Var1, gap = Var2) %>% {
  ggplot(., aes(gap, value)) + geom_line(aes(group = iter), alpha = 0.3) +
    xlab("Weeks between samples") +
    ylab("Probability animal remains in population")
}

## -- survival over season -- ##
## -- Calculate offseason survival probability as a function of individual bd deviate -- ##

ind.pred <- reshape2::melt(
  plogis(sweep(stan.fit.samples$bd_ind, 1, stan.fit.samples$beta_offseason, "*"))
)

ind.pred %>% filter(Var2 < 10) %>% {
  ggplot(., aes(x = value)) + geom_histogram(bins = 50) + 
    facet_wrap(~Var2)
}

hist(ind.pred$value, breaks = 1000)

stan.ind_pred_eps <- stan.fit.samples$bd_delta_eps %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value)
stan.ind_pred_var <- stan.fit.samples$bd_delta_sigma %>%
  reshape2::melt(.) %>% left_join(., stan.ind_pred_eps) %>%
  mutate(eps = eps * value) %>% group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% arrange(mid) %>%
  mutate(ind = factor(ind, levels = ind))

## Not doing a very good job of making these look different than 0
stan.ind_pred_var %>% {
  ggplot(., aes(as.factor(ind), mid)) + geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr)) +
    xlab("Individual") + 
    ylab("Random Effect Deviate") +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = 1, colour = "firebrick3") +
    theme(
      axis.text.x = element_text(size = 8)
    )
}

ind_order <- expdat.all %>% group_by(ind) %>% summarize(
  tot_bd = sum(log_bd_load)
) %>% arrange(desc(tot_bd)) %>% mutate(order_real = seq(n()), ind = as.character(ind)) %>% 
  left_join(.
    , {
      stan.ind_pred_var %>% arrange(desc(mid)) %>% mutate(order_pred = seq(n()), ind = as.character(ind))
    })

ind.pred %>% filter(Var2 %in% c(ind_order$ind[c(1:10)])) %>% {
  ggplot(., aes(x = value)) + 
    geom_histogram(bins = 50, fill = "dodgerblue3", alpha = 0.5) + 
    facet_wrap(~Var2)
}

ind.pred %>% filter(Var2 %in% c(ind_order$ind[c(87:97)])) %>% {
  ggplot(., aes(x = value)) + 
    geom_histogram(bins = 50, fill = "firebrick3", alpha = 0.5) + 
    facet_wrap(~Var2)
}


### Empirical debugging

ind_order <- capt_history %>% group_by(Mark) %>% summarize(
  tot_bd = sum(log_bd_load)
) %>% arrange(desc(tot_bd)) %>% mutate(order_real = seq(n()), Mark = as.character(Mark)) %>% 
  left_join(.
    , {
      stan.ind_pred_var %>% arrange(desc(mid)) %>% 
        mutate(order_pred = seq(n()), ind = as.character(ind)) %>% rename(Mark = ind)
    })

low_ind  <- unique((capt_history %>% filter(Mark %in% (ind_order %>% arrange(order_real) %>% tail(10))$Mark))$Mark)
high_ind <- unique((capt_history %>% filter(Mark %in% (ind_order %>% arrange(order_real) %>% head(10))$Mark))$Mark)

capt_history %>% filter(Mark %in% (ind_order %>% arrange(order_real) %>% tail(10))$Mark) %>% 
  filter(swabbed == 1) %>% {
  ggplot(., aes(week, bd_load)) +
    geom_line(aes(colour = as.factor(Mark))) + 
    geom_point(aes(colour = as.factor(Mark))) + 
      scale_y_log10() +
    facet_wrap(~year) 
  }

ind_order$ind_type <- 0
ind_order[ind_order$Mark %in% low_ind, ]$ind_type  <- 1
ind_order[ind_order$Mark %in% high_ind, ]$ind_type <- 2

ggplot(ind_order, aes(order_real, order_pred)) + 
  geom_point(aes(colour = as.factor(ind_type)), size = 2) +
  scale_color_brewer(palette = "Dark2", name = "Individual type", labels = c("Middle", "Lowest Bd", "Highest Bd")) +
  xlab("Individual Ordered by Average of Measured Bd") +
  ylab("Predicted Order")

ind.pred %>% filter(Var2 %in% low_ind) %>% {
  ggplot(., aes(x = value)) + 
    geom_histogram(bins = 50, fill = "dodgerblue3", alpha = 0.5) + 
    facet_wrap(~Var2)
}

ind.pred %>% filter(Var2 %in% high_ind) %>% {
  ggplot(., aes(x = value)) + 
    geom_histogram(bins = 50, fill = "firebrick3", alpha = 0.5) + 
    facet_wrap(~Var2)
}

###### ###### ###### ###### ###### ###### ###### 

stan.ind_pred_eps <- stan.fit.samples$phi_delta_eps %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value)
stan.ind_pred_var <- stan.fit.samples$phi_delta_sigma %>%
  reshape2::melt(.) %>% left_join(., stan.ind_pred_eps) %>%
  mutate(eps = eps * value) %>% group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% arrange(mid) %>%
  mutate(ind = factor(ind, levels = ind))

## Not doing a very good job of making these look different than 0
stan.ind_pred_var %>% {
  ggplot(., aes(as.factor(ind), mid)) + geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr)) +
    xlab("Individual") + 
    ylab("Random Effect Deviate") +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = 1, colour = "firebrick3") +
    theme(
      axis.text.x = element_text(size = 8)
    )
}

