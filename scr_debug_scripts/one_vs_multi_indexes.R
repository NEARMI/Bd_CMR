
stan.fit <- readRDS(paste(paste("fits/stan_fit", which.dataset, sep = "_"), "Rds", sep = "."))
stan.fit <- readRDS(paste(paste("fits/stan_fit", "ANBO", sep = "_"), "Rds", sep = "."))

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

stan.fit.summary.s <- stan.fit.summary
stan.fit.samples.s <- stan.fit.samples

stan.fit.summary.m <- stan.fit.summary
stan.fit.samples.m <- stan.fit.samples

##########

stan.fit.summary.s[grep("beta", dimnames(stan.fit.summary.s)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%')

stan.fit.summary[grep("beta", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') 

#########

compare_pops <- data.frame(
  pop = c(rep("single", 1500), rep("multi", 1500))
, est = c(
  stan.fit.samples.s$beta_offseason[, 2]
, (stan.fit.samples.m$beta_offseason[, 1] + stan.fit.samples.m$offseason_pop_bd[, 5])))

compare_pops %>% {
  ggplot(., aes(x = est)) + 
    geom_histogram(aes(colour = pop, fill = pop), alpha = 0.4, bins = 50) + 
    facet_wrap(~pop, ncol = 1)
}

compare_pops <- data.frame(
  pop = c(rep("single", 1500), rep("multi", 1500))
, est = c(
  stan.fit.samples.s$beta_phi
, (stan.fit.samples.m$beta_inseason[, 1] + stan.fit.samples.m$inseason_pop[, 5])))

compare_pops %>% {
  ggplot(., aes(x = est)) + 
    geom_histogram(aes(colour = pop, fill = pop), alpha = 0.4, bins = 50) + 
    facet_wrap(~pop, ncol = 1)
}

compare_pops <- data.frame(
  pop = c(rep("single", 1500), rep("multi", 1500))
, est = c(
  stan.fit.samples.s$beta_offseason[, 3]
, stan.fit.samples.m$beta_offseason[, 2]
  ))

compare_pops %>% {
  ggplot(., aes(x = est)) + 
    geom_histogram(aes(colour = pop, fill = pop), alpha = 0.4, bins = 50) + 
    facet_wrap(~pop, ncol = 1)
}

###########

which_pop_days <- (capt_history.p %>% group_by(date_fac) %>% slice(1))$pop_spec %>% as.numeric()

stan.fit.samples.s$p_day_dev %>% 
  reshape2::melt() %>% 
  group_by(Var2) %>% summarize(
    mid_s = quantile(value, 0.50)
  ) %>% rename(day = Var2) %>% left_join(
    .
    , 
  stan.fit.samples.m$p_day_dev[, which(which_pop_days == 5)] %>% 
  reshape2::melt() %>% 
  group_by(Var2) %>% summarize(
    mid_m = quantile(value, 0.50)
  ) %>% rename(day = Var2)
  ) %>% {
    ggplot(., aes(mid_s, mid_m)) + geom_point() + geom_abline(intercept = 0, slope = 1)
  }

###########

which_pop_ind <- (capt_history.p %>% group_by(Mark) %>% slice(1))$pop_spec %>% as.numeric()

stan.fit.samples.s$bd_ind %>% 
  reshape2::melt() %>% 
  group_by(Var2) %>% summarize(
    mid_s = quantile(value, 0.50)
  ) %>% rename(day = Var2) %>% left_join(
    .
    , 
  stan.fit.samples.m$bd_ind[, which(which_pop_ind == 5)] %>% 
  reshape2::melt() %>% 
  group_by(Var2) %>% summarize(
    mid_m = quantile(value, 0.50)
  ) %>% rename(day = Var2)
  ) %>% {
    ggplot(., aes(mid_s, mid_m)) + geom_point() + geom_abline(intercept = 0, slope = 1)
  }

##########

capt_history.r    <- capt_history %>% filter(Species == "ANBO") %>% droplevels()
mark_check.all    <- capt_history %>% group_by(pop_spec) %>% summarize(n_mark = n_distinct(Mark))
which_mark_single <- 7
mark_check        <- mark_check.all$n_mark[1] + mark_check.all$n_mark[2] + 
  mark_check.all$n_mark[3] + mark_check.all$n_mark[4] + which_mark_single

data.frame(
   mult_pop = stan.fit.samples.m$phi[, which(capt_history.phi$Mark == mark_check)] %>% colMeans()
 , sing_pop = stan.fit.samples.s$phi[, ((which_mark_single - 1)*19 + 1):(19 * which_mark_single)] %>% colMeans()
) %>% {
  ggplot(., aes(sing_pop, mult_pop)) + geom_jitter()
}

