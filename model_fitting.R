################################################
## Try and fit a CMR model with the newt data ##
################################################

####
## About this script
####

## 1) see and start with 'data_exploration.R' for data info
## 2) see 'CMR_simulations.R' for a CMR model with a bd submodel that fits reasonably well with simulated data
## 3) see data_structure.Rmd / html for the data structures needed to fit the model

####
## Notes as of OCT 20:
####

## After some reasonable success with fitting the model to simulated data, time to start on real data. 
 ## First step is to get the data into the correct structure...

####
## Packages and functions
####
needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan")
lapply(needed_packages, require, character.only = TRUE)
source("../ggplot_theme.R")
set.seed(10002)
'%notin%' <- Negate('%in%')

Bd_Newts_AllSites   <- read.csv("Bd_Newts_AllSites_10.1.21.csv")

## Stupid dates
if (length(grep("/", Bd_Newts_AllSites$Date[1])) > 0) {

date_convert <- apply(matrix(Bd_Newts_AllSites$Date), 1, FUN = function (x) {
  a <- strsplit(x, split = "/")[[1]]
  b <- a[3]
  b <- strsplit(b, "")[[1]][c(3, 4)] %>% paste(collapse = "")
  paste(c(a[c(1, 2)], b), collapse = "/")
})

Bd_Newts_AllSites$Date <- date_convert
Bd_Newts_AllSites      %<>% mutate(Date = as.Date(Date, "%m/%d/%y"))

} else {
  
Bd_Newts_AllSites      %<>% mutate(Date = as.Date(Date))
  
}

## Some parameters
red_ind        <- FALSE    ## TRUE   = reduce the number of individuals for debugging purposes

if (red_ind) {
red_total_capt <- 4        ## minimum number of recaptures to keep an individual in an analysis
}

## Just select one site for now for a trial fit
A11 <- Bd_Newts_AllSites %>% 
  filter(Site == "A11") %>% 
  group_by(year) %>% 
  mutate(week = ceiling(julian / 7))

## Find the first and last week ever sampled in this population
week_range <- A11 %>% 
  summarize(
    min_week = min(week)
  , max_week = max(week)
  ) %>% summarize(
    min_week = min(min_week)
  , max_week = max(max_week)
  ) %>% unlist()

## Find the unique weeks sampled in each year
sampled_weeks <- A11 %>% 
  summarize(week = unique(week)) %>%
  mutate(sampled = 1)

## Construct an "all possible combinations" data frame and parse it down
capt_history <- expand.grid(
  week = seq(from = week_range["min_week"], to = week_range["max_week"], by = 1)
, year = unique(A11$year)
, Mark = unique(A11$Mark)
) %>% left_join(., sampled_weeks) %>% 
  mutate(sampled = ifelse(is.na(sampled), 0, 1)) %>%
  left_join(.
    , {A11 %>% dplyr::select(week, year,  Mark, month, copies.swab)}
    ) %>% rename(
      captured = month
    , bd_load  = copies.swab) %>% 
    mutate(
     captured = ifelse(is.na(captured), 0, 1)
   , swabbed  = ifelse(is.na(bd_load), 0, 1)) %>%
    group_by(Mark, week, year) %>%
    summarize(
      sampled  = sum(sampled)
    , captured = sum(captured, na.rm = T)
    , swabbed  = sum(swabbed, na.rm = T)
    , bd_load  = sum(bd_load, na.rm = T)
    ) %>% 
   mutate(
      sampled  = ifelse(sampled > 1, 1, sampled)
    , captured = ifelse(captured > 1, 1, captured)
    , swabbed  = ifelse(swabbed > 1, 1, swabbed)
    , log_bd_load = log(bd_load + 1)                           ### eeek!
    , log_bd_load = ifelse(is.na(log_bd_load), 0, log_bd_load) ### eeek X2!!
    ) %>% 
   mutate(Mark = as.factor(Mark)) %>% 
   mutate(Mark = as.numeric(Mark))

## Add a column for each unique sampling week and all weeks
capt_history %<>% ungroup() %>% arrange(year, week, Mark) %>% mutate(
   week_year  = interaction(week, year)
 , week_year  = as.factor(week_year)
 , week_year  = as.numeric(week_year)
 , year_f     = as.numeric(as.factor(year)) - 1
 , cont_weeks = (52 * year_f) + week
)

## Quick double check to see how collapsing to week held up
capt_history %>% filter(swabbed == 1) %>%  {
  ggplot(., aes(week, log_bd_load)) + 
    geom_line(aes(group = Mark)) +
    facet_wrap(~year) +
    xlab("Week") +
    ylab("Bd Load") 
}

## Try a small subset to see if the model will run (regardless of whether it will fit or not)
if (red_ind) {
  capt_history.well_meas <- capt_history %>% 
  group_by(Mark) %>%
  summarize(total_capt = sum(captured)) %>% 
  filter(total_capt > red_total_capt) %>% 
  droplevels()
  
  capt_history %<>% filter(Mark %in% capt_history.well_meas$Mark)
} 

## individuals' measured bd 
capt_history.bd_load <- capt_history %>% 
  ungroup() %>%
  arrange(Mark, week_year) %>%
  filter(swabbed == 1)

## first and last _OF THE CAPTURE EVENTS_ in which each individual was captured
capture_range  <- capt_history %>% 
  group_by(Mark) %>% 
  filter(sampled == 1) %>%  
  summarize(
    first = min(which(captured == 1))
  , final = max(which(captured == 1))) %>% 
  dplyr::select(first, final) %>% 
  mutate(
    first = ifelse(is.infinite(first), 0, first)
  , final = ifelse(is.infinite(final), 0, final)
    )

capture_total <- capt_history %>% 
  filter(sampled == 1) %>% 
  group_by(week_year) %>% 
  summarize(total_capt = sum(captured))

capt_history %>% mutate(event = week) %>% {
  ggplot(., aes(week, Mark, fill = as.factor(captured))) + 
    geom_tile(aes(alpha = sampled)) +
    geom_point(data = capt_history %>% mutate(event = week) %>% 
        filter(swabbed == 1), aes(x = week, y = Mark, z = NULL), lwd = 0.7) +
    facet_wrap(~year) +
    xlab("Week of the year") +
    ylab("Individual") +
    scale_fill_manual(
        values = c("dodgerblue4", "firebrick4")
      , name   = "Detected?"
      , labels = c("No", "Yes")) +
    guides(alpha = FALSE) +
    theme(
      axis.text.y = element_text(size = 6)
  #  , axis.text.x = element_text(size = 8, angle = 300)
    , legend.text = element_text(size = 12)
    , legend.key.size = unit(.55, "cm")
    ) 
}

####
## Data in the needed structure for the stan model
####

## Numbers and lengths of things
n_periods <- length(unique(capt_history$year))
n_ind     <- length(unique(capt_history$Mark))  
n_times.w <- length(seq(week_range[1], week_range[2]))
n_times.a <- length(seq(week_range[1], week_range[2])) * n_periods
n_occ     <- sampled_weeks %>% group_by(year) %>%
  summarize(n_occ = length(unique(week)))
n_occ     <- n_occ$n_occ

## Vectors for detection
capt_history.p   <- capt_history %>% filter(sampled == 1) %>% 
  arrange(Mark) %>% ungroup()

p_first_index <- (capt_history.p %>% mutate(index = seq(n())) %>% 
  group_by(Mark) %>% 
  summarize(first_index = min(index)))$first_index

## check to make sure things are aligning properly
if ((p_first_index[2] - p_first_index[1]) != sum(n_occ)) {
  print("Gaps between p indexes doesn't match number of sampling events")
}

## determine the first period in which each individual was present
first_capt <- capt_history.p %>% group_by(Mark, year) %>% 
  summarize(capt = sum(captured)) %>% 
  mutate(capt = cumsum(capt)) %>% 
  mutate(capt = ifelse(capt > 0, 1, 0)) 

## time periods in which we do not know if an individual was present or not
p_zeros <- matrix(data = 0, nrow = n_ind, ncol = sum(n_occ))
for (i in 1:n_ind) {
  p_zeros[i, ] <- rep(first_capt[first_capt$Mark == unique(first_capt$Mark)[i], ]$capt, n_occ)
}
p_zeros   <- (p_zeros %>% reshape2::melt() %>% arrange(Var1))$value

capt_history.p$p_zeros <- p_zeros

## Vectors for survival
last_week        <- max(capt_history[capt_history$sampled == 1, ]$week_year)
capt_history.phi <- capt_history %>% filter(week_year != last_week, sampled == 1) %>% 
  arrange(Mark) %>% ungroup()

## Determine the number of time periods that elapse between back to back samples
the_weeks <- capt_history %>% 
  filter(sampled == 1) %>% 
  summarize(cont_weeks = unique(cont_weeks))
the_weeks <- the_weeks$cont_weeks
time_gaps <- (the_weeks - lag(the_weeks, 1))[-1]

## Weeks between sampling events
capt_history.phi %<>% mutate(time_gaps = rep(time_gaps, n_ind))

## Offseason vector (not the best strategy, but ok)
capt_history.phi %<>% mutate(offseason = ifelse(time_gaps > 20, 1, 0))

phi_first_index <- (capt_history.phi %>% mutate(index = seq(n())) %>% group_by(Mark) %>% 
  summarize(first_index = min(index)))$first_index

## check to make sure things are aligning properly
if ((phi_first_index[2] - phi_first_index[1]) != (sum(n_occ) - 1)) {
  print("Gaps between phi indexes doesn't match number of sampling events minus 1")
}

## Indices for which entries of phi must be 0
phi_zeros <- matrix(data = 0, nrow = n_ind, ncol = sum(n_occ) - 1)
for (i in 1:n_ind) {
  phi_zeros[i, ] <- c(
    rep(1, capture_range$first[i] - 1)
  , rep(0, ncol(phi_zeros) - (capture_range$first[i] - 1))
    )
}
phi_zeros   <- (phi_zeros %>% reshape2::melt() %>% arrange(Var1))$value

capt_history.phi$phi_zeros <- phi_zeros

####
## Other needed covaraites 
####

## OCT 20: Hmm, not good. Temp will need to be imputed?
temp <- expand.grid(
  year = unique(A11$year)
, week = seq(week_range[1], week_range[2])
)

temp_have <- A11 %>% 
  group_by(year, week) %>% 
  summarize(temp = mean(Site_temp, na.rm = T))

## For the first run just get 2018 temp to be 2020 temp
temp_have[1:8, 3] <- temp_have[c(21, 9, 10, 11, 12, 15, 18, 21), 3]

temp <- left_join(temp, temp_have)

## For now do a really ugly imputation to just get the model running
temp %<>% group_by(week) %>% mutate(temp = ifelse(is.na(temp), mean(temp, na.rm = T), temp))
temp[c(7:9), ]$temp   <- (temp[c(10:12), ]$temp - temp[c(4:6), ]$temp) + temp[c(4:6), ]$temp
temp[c(22:24), ]$temp <- (temp[c(25:27), ]$temp - temp[c(19:21), ]$temp) + temp[c(19:21), ]$temp

temp %<>% arrange(year, week) 
temp <- as.matrix(temp[, 3])

####
## Run the stan model
####

stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop           = 1
 , n_periods       = n_periods
 , n_ind           = n_ind
 , n_times         = n_times.a
 , times_within    = n_times.w
 , ind_occ         = sum(n_occ) * n_ind
 , ind_occ_min1    = (sum(n_occ) - 1) * n_ind
  
  ## short vector indexes 
 , time              = rep(seq(n_times.w), n_periods)
 , time_per_period   = matrix(data = seq(n_times.w * n_periods), nrow = n_times.w, ncol = n_periods)
 , periods           = rep(seq(n_periods), each = n_times.w)
 , ind_occ_size      = rep(sum(n_occ), n_ind)
 , ind_occ_min1_size = rep(sum(n_occ) - 1, n_ind)
 , ind_in_pop        = rep(1, n_ind)

 , phi_first_index   = phi_first_index
 , p_first_index     = p_first_index
  
  ## long vector indexes
 , ind_occ_min1_rep    = capt_history.phi$Mark
 , sampling_events_phi = capt_history.phi$week_year
 , offseason           = capt_history.phi$offseason
 , pop_phi             = rep(1, nrow(capt_history.phi))
 , phi_zeros           = capt_history.phi$phi_zeros

 , ind_occ_rep       = capt_history.p$Mark
 , sampling_events_p = capt_history.p$week
 , periods_occ       = as.numeric(as.factor(capt_history.p$year))
 , pop_p             = rep(1, nrow(capt_history.p))
 , p_zeros           = capt_history.p$p_zeros

  ## covariates
 , N_bd            = nrow(capt_history.bd_load)
 , X_bd            = capt_history.bd_load$log_bd_load  
 , ii_bd           = capt_history.bd_load$Mark
 , tt_bd           = capt_history.bd_load$week_year
 , temp            = temp
 , time_gaps       = capt_history.phi$time_gaps
  
  ## Capture data
 , N_y             = nrow(capt_history.p)
 , y               = capt_history.p$captured
  
 , first           = capture_range$first
 , last            = capture_range$final

  )

stan.fit  <- stan(
  file    = "CMR_simulation/CMR_empirical.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )
  
shinystan::launch_shinystan(stan.fit)

saveRDS(stan.fit, "stan.fit.empirical.Rds")

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

bd_levels <- log(c(seq(1, 10, by = 0.5) %o% 10^(0:5)))

####
## Recovery of simulated coefficients?
##   *NOTE*: High success with keeping all individuals, only moderate success with dropping individuals
####

## Primary Bd effects
pred_coef        <- as.data.frame(stan.fit.summary[1:4, c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(param = rownames(.))

pred_coef %>% {
  ggplot(., aes(param, mid)) + geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    xlab("Parameter") + ylab("Estimate") +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.5) +
    scale_x_discrete(labels = c(
      "Detection intercept"
    , "Detection slope"
    , "Survival intercept"
    , "Survival slope"
      )) +
    theme(axis.text.x = element_text(size = 11))
}

## survival over Bd load
stan.pred        <- apply(stan.fit.samples$beta_phi, 1
  , FUN = function(x) plogis(x[1] + x[2] * bd_levels)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "mortality")
stan.pred        %<>% mutate(log_bd_load = plyr::mapvalues(log_bd_load
  , from = unique(log_bd_load), to = bd_levels))
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
    xlab("Log of Bd Load") + ylab("Predicted mortality probability")
}

## detection over Bd load
stan.pred        <- apply(stan.fit.samples$beta_p, 1
  , FUN = function(x) plogis(x[1] + x[2] * bd_levels)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "detect")
stan.pred        %<>% mutate(log_bd_load = plyr::mapvalues(log_bd_load
  , from = unique(log_bd_load), to = bd_levels))
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
    xlab("Log of Bd Load") + ylab("Predicted detection probability")
}

## Change in Bd
pred_coef <- as.data.frame(stan.fit.summary[
  c(
    grep("bd_delta_mu", dimnames(stan.fit.summary)[[1]])
  , grep("bd_delta_sigma", dimnames(stan.fit.summary)[[1]])
  , grep("bd_add", dimnames(stan.fit.summary)[[1]])
  , grep("bd_obs", dimnames(stan.fit.summary)[[1]])
  , grep("start_mean", dimnames(stan.fit.summary)[[1]])
  , grep("start_var", dimnames(stan.fit.summary)[[1]])
  )
  , c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(param = rownames(.))

## Hard to directly compare coefficient estimates given the change in sccale, need to work on this
pred_coef %>% filter(param != "start_mean") %>% {
  ggplot(., aes(param, mid)) + geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    xlab("Parameter") + ylab("Estimate") +
    theme(axis.text.x = element_text(size = 11))
}

####
## What does predicted overall bd look like?
####
stan.pred     <- matrix(nrow = length(unique(capt_history$Date)), ncol = 1000, data = 0)
stan.pred.ind <- array(dim = c(length(unique(capt_history$Mark)), length(unique(capt_history$Date)), 1000), data = 0)
samp_occ      <- seq(length(unique(capt_history$Date)))

for (i in 1:nrow(stan.pred)) {
stan.pred[i, ] <- stan.fit.samples$beta_bd[, 1] + 
  stan.fit.samples$beta_bd[, 2] * samp_occ[i] +
  stan.fit.samples$beta_bd[, 3] * samp_occ[i]^2

for (j in 1:dim(stan.pred.ind)[1]) {
stan.pred.ind[j,i,] <- stan.fit.samples$bd_ind[, j] +
  stan.fit.samples$beta_bd[, 2] * samp_occ[i] +
  stan.fit.samples$beta_bd[, 3] * samp_occ[i]^2
}

}

stan.pred <- reshape2::melt(stan.pred)
names(stan.pred) <- c("occ", "iter", "value")

stan.pred %<>% 
  group_by(occ) %>%
  summarize(
    lwr = quantile(value, 0.025)
  , mid = quantile(value, 0.50)
  , upr = quantile(value, 0.975)
  )

stan.pred.ind <- reshape2::melt(stan.pred.ind)
names(stan.pred.ind) <- c("ind", "occ", "iter", "value")

stan.pred.ind %<>% 
  group_by(occ, ind) %>%
  summarize(
    lwr = quantile(value, 0.025)
  , mid = quantile(value, 0.50)
  , upr = quantile(value, 0.975)
  )

ggplot(stan.pred, aes(occ, mid)) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  geom_line() +
  geom_line(data = stan.pred.ind
    , aes(occ, mid, group = ind), lwd = 0.5, alpha = 0.5) +
  geom_line(
    data = (capt_history %>% 
        filter(swabbed == 1) %>% 
        mutate(occ = week))
  , aes(occ, log_bd_load, group = Mark)
  , colour = "red"
  ) + xlab("Sampling Occasion") +
  ylab("Bd Load")

####
## Add individual random effect estimates
####

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

## Random effects in phi and p
stan.ind_pred_eps <- stan.fit.samples$eps_alpha_phi %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value)
stan.ind_pred_var <- stan.fit.samples$sigma_alpha_phi %>%
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
