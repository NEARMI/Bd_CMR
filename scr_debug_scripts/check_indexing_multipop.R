## Checking phi, p, and chi indexing for individuals only captured once

 ## !!!! This will be an important debugging script to show in the final repo so dont delete

stan.fit <- readRDS(paste("fits/stan_fit_multipop_2022-03-26", "Rds", sep = "."))

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

shinystan::launch_shinystan(stan.fit)

shinystan::launch_shinystan(shinystan::as.shinystan(stan.fit, pars = names(stan.fit.samples)[-c(6, 26, 28, 36, 37, 40, 41, 45, 46, 53, 54, 55)]))

capt_history.r <- capt_history %>% filter(pop_spec == which.dataset)

mark_check.all <- capt_history %>% group_by(pop_spec) %>% summarize(n_mark = n_distinct(Mark))

## Mark
mark_check <- mark_check.all$n_mark[1] + 
  mark_check.all$n_mark[2] + 
  mark_check.all$n_mark[3] +
  mark_check.all$n_mark[4] +
  mark_check.all$n_mark[5] +
  mark_check.all$n_mark[6] + 

mark_check <- sample(sum(mark_check.all$n_mark), 1)
capt_history.phi %>% 
  filter(Mark == mark_check) %>% 
  dplyr::select(-Month, -Year, -month_year, -pop_spec, -Region, -State, -Site, -PrimNum) %>%
  mutate(pred_phi = stan.fit.samples$phi[, which(capt_history.phi$Mark == mark_check)] %>% colMeans()) %>%
  as.data.frame()

mark_check <- sample(sum(mark_check.all$n_mark), 1)
capt_history %>% 
  filter(Mark == mark_check) %>% 
  dplyr::select(-Month, -Year, -month_year, -pop_spec, -Region, -State, -Site) %>%
  mutate(
    pred_phi = c(stan.fit.samples$phi[, which(capt_history.phi$Mark == mark_check)] %>% colMeans(), 1)
  , pred_p   = stan.fit.samples$p[, which(capt_history.p$Mark == mark_check)] %>% colMeans()
  , pred_chi = stan.fit.samples$chi[, which(capt_history.p$Mark == mark_check)] %>% colMeans()) %>% 
  as.data.frame() 

## Loop from first to last observation
first_p_one <- (capture_range$first[mark_check] + 1)
last_e      <- capture_range$final[mark_check]

## entries that are 1 ~ bernoulli( )
entry_one <- phi_first_index[mark_check] - 1 + first_p_one - 1
entry_two <- phi_first_index[mark_check] - 1 + last_e - 1

stan.fit.samples$phi[1:5, which(capt_history.phi$Mark == mark_check)]
stan.fit.samples$phi[1:5, entry_one]
stan.fit.samples$phi[1:5, entry_two]

## entries for y ~ bernoulli
entry_one <- p_first_index[mark_check] - 1 + first_p_one
entry_two <- p_first_index[mark_check] - 1 + last_e

stan.fit.samples$p[1:5, which(capt_history.p$Mark == mark_check)]
stan.fit.samples$p[1:5, entry_one]
stan.fit.samples$p[1:5, entry_two]

## entries of chi
stan.fit.samples$chi[, which(capt_history.p$Mark == mark_check)] %>% colMeans()

p_first_index[mark_check]:(p_first_index[mark_check] + 25 - 1)


p_first_index[mark_check] - 1 + last_e

stan.fit.samples$chi[1:5, which(capt_history.p$Mark == mark_check)]
stan.fit.samples$chi[1:5, p_first_index[mark_check] - 1 + last_e]


####
## Checking estimated individual sizes and population MeHg
####

capt_history %>% group_by(pop_spec) %>% summarize(
  mean(merc, na.rm = T)
)

stan.fit.samples$beta_mehg_spec[, 3] %>% hist()

stan.fit.samples$mehg_pop[, 1] %>% hist()

stan.fit.samples$mehg_pop_est[, 1] %>% hist()

len.mis %>% as.array()

stan.fit.samples$ind_len[, 689] %>% hist()
stan.fit.samples$ind_len_scaled[, 688] %>% hist()


####
## Checking each day detection deviate
####

stan.fit.samples$p_per_day %>% reshape2::melt() %>% 
  rename(day = Var2) %>%
  group_by(day) %>% 
  summarize(
    lwr   = quantile(value, 0.025)
  , lwr_n = quantile(value, 0.200)
  , mid   = quantile(value, 0.500)
  , upr_n = quantile(value, 0.800)
  , upr   = quantile(value, 0.975)
  ) %>% mutate(
    pop = (capt_history.p %>% group_by(date_fac) %>% slice(1))$pop_spec %>% as.numeric()
  ) %>% {
    ggplot(., aes(mid, day)) +
      geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0, size = 1) +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.5) +
      geom_point() +
      facet_wrap(~pop, scales = "free")
  }


####
## And check back with all of the beta estimates as well
####

stan.fit.summary[grep("beta", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>% as.data.frame()


####
## What is going on with Bd in each population?
####

unique(capt_history$pop_spec)
spectoplot <- 2
poptoplot  <- 6

outval   <- matrix(seq(1, 14, by = 1))
out.pred <- matrix(nrow = dim(stan.fit.samples[[1]])[1], ncol = length(outval))

for (j in 1:ncol(out.pred)) {
   out.pred[, j] <- plogis(
     stan.fit.samples$beta_offseason_int[, spectoplot] + stan.fit.samples$offseason_pop[, poptoplot] +
       (stan.fit.samples$beta_offseason_bd[, spectoplot] + stan.fit.samples$offseason_pop_bd[, poptoplot]) * outval[j] +
       (stan.fit.samples$beta_offseason_len[, spectoplot] + stan.fit.samples$offseason_pop_len[, poptoplot]) * 0 +
       stan.fit.samples$beta_offseason_mehg[, spectoplot] * stan.fit.samples$mehg_pop_est_scaled[, poptoplot]
    )
}

out.pred <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "Bd")

out.pred.off <- out.pred %>%
  group_by(gap, variable) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% mutate(
    when       = "offseason"
  )

out.pred.off %>% {ggplot(., aes(gap, mid)) + geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) + geom_line()}
