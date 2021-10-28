################################################
## Try and fit a CMR model with the newt data ##
################################################

####
## About this script
####

## Extension of model_fitting.R for multiple populations

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

Bd_Newts_AllSites %>% group_by(Site) %>% summarize(n_y = length(unique(year)))
unique((Bd_Newts_AllSites %>% filter(Site == "A04"))$year)

## Just select one site for now for a trial fit
A11 <- Bd_Newts_AllSites %>% 
  filter(SA == "PA") %>% 
  filter(Site == "A04" | Site == "A11" | Site == "A05" | Site == "A10") %>%
  group_by(year) %>%  
  mutate(week = ceiling(julian / 7)) %>% 
  filter(!is.na(Site)) %>% 
  arrange(Site)

n_sites <- length(unique(A11$Site))
u_sites <- unique(A11$Site)

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
  group_by(Site, year) %>%
  summarize(week = unique(week)) %>%
  mutate(sampled = 1) %>%
  arrange(Site, week) 

## Construct an "all possible combinations" data frame and parse it down
all_ind <- 0
for (i in 1:n_sites) {
  
capt_history.t <- expand.grid(
  week = seq(from = week_range["min_week"], to = week_range["max_week"], by = 1)
, year = unique(A11$year)
, Site = u_sites[i]
, Mark = unique((A11 %>% filter(Site == u_sites[i]))$Mark)) %>% 
  left_join(., (sampled_weeks %>% filter(Site == u_sites[i]))) %>% 
  left_join(.
  , (A11 %>% filter(Site == u_sites[i]) %>%
      dplyr::select(week, year,  Mark, month, copies.swab)
     )
  ) %>% mutate(sampled = ifelse(is.na(sampled), 0, 1)) %>%
   rename(
      captured = month
    , bd_load  = copies.swab) %>% 
    mutate(
     captured = ifelse(is.na(captured), 0, 1)
   , swabbed  = ifelse(is.na(bd_load), 0, 1)) %>%
    group_by(Mark, week, year, Site) %>%
    summarize(
      sampled  = sum(sampled)
    , captured = sum(captured, na.rm = T)
    , swabbed  = sum(swabbed, na.rm = T)
    , bd_load  = sum(bd_load, na.rm = T)
    ) %>% 
   mutate(
      sampled     = ifelse(sampled > 1, 1, sampled)
    , captured    = ifelse(captured > 1, 1, captured)
    , swabbed     = ifelse(swabbed > 1, 1, swabbed)
    , log_bd_load = log(bd_load + 1)                           ### eeek!
    , log_bd_load = ifelse(is.na(log_bd_load), 0, log_bd_load) ### eeek X2!!
    )

capt_history.t %<>% mutate(Mark = as.factor(Mark)) %>%
  mutate(Mark = as.numeric(Mark))
n_inds <- max(capt_history.t$Mark)
capt_history.t %<>% mutate(Mark = Mark + all_ind)
all_ind <- all_ind + n_inds

capt_history.t %<>% ungroup() %>% 
  arrange(year, week, Mark, Site) %>% 
  mutate(
   week_year  = interaction(week, year)
 , week_year  = as.factor(week_year)
 , week_year  = as.numeric(week_year)
 , year_f     = as.numeric(as.factor(year)) - 1
 , cont_weeks = (52 * year_f) + week
)

if (i == 1) {
capt_history <- capt_history.t
} else {
capt_history <- rbind(capt_history, capt_history.t)
}
   
print(paste("Through", i, "of", n_sites, "sites", sep = " "))

}
 
capt_history %<>% arrange(Mark, year, Site, week)

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
    facet_grid(Site~year) +
    xlab("Week of the year") +
    ylab("Individual") +
    scale_fill_manual(
        values = c("dodgerblue4", "firebrick4")
      , name   = "Detected?"
      , labels = c("No", "Yes")) +
    guides(alpha = FALSE) +
    theme(
      axis.text.y = element_text(size = 6)
    , legend.text = element_text(size = 12)
    , legend.key.size = unit(.55, "cm")
    ) 
}

test_ind <- capt_history %>% group_by(Site) %>%
  summarize(n_ind = length(unique(Mark))) %>% 
  dplyr::select(-Site) %>% as.matrix()
test_ind <- ifelse(test_ind > 100, 100, test_ind)

rand_inds <- matrix(
  nrow = 4,
  ncol = 15
)
for (i in 1:n_sites) {
  temp_dat       <- capt_history %>% filter(Site == u_sites[i])
  rand_inds[i, ] <- sort(sample(unique(temp_dat$Mark), test_ind[i, 1]))
  temp_dat  <- temp_dat %>% filter(Mark %in% rand_inds[i, ])
  if (i == 1) {
    temp_dat.a <- temp_dat
  } else {
   temp_dat.a <- rbind(temp_dat.a, temp_dat) 
  }
}

rand_inds     <- unique(temp_dat.a$Mark)
capt_history  <- temp_dat.a %>% mutate(Mark = as.factor(Mark)) %>% mutate(Mark = as.numeric(Mark))
capture_range <- capture_range[rand_inds, ]
capt_history.bd_load <- capt_history %>% 
  ungroup() %>%
  arrange(Mark, week_year) %>%
  filter(swabbed == 1)

####
## Data in the needed structure for the stan model
####

## Numbers and lengths of things
n_periods <- length(unique(capt_history$year))
n_ind     <- length(unique(capt_history$Mark)) 
n_ind.per <- capt_history %>% group_by(Site) %>%
  summarize(n_ind = length(unique(Mark))) %>% 
  dplyr::select(-Site) %>% as.matrix()
n_times.w <- length(seq(week_range[1], week_range[2]))
n_times.a <- length(seq(week_range[1], week_range[2])) * n_periods
n_occ     <- sampled_weeks %>% group_by(year, Site) %>%
  summarize(n_occ = length(unique(week))) %>% 
  mutate(Site = factor(Site, levels = u_sites)) %>%
  arrange(Site) %>%
  pivot_wider(values_from = n_occ, names_from = Site) %>% 
  arrange(year) %>%
  ungroup() %>%
  dplyr::select(-year) %>% as.matrix()
n_occ[is.na(n_occ)] <- 0

## Vectors for detection
capt_history.p   <- capt_history %>% 
  filter(sampled == 1) %>% 
  ungroup()

p_first_index <- (capt_history.p %>% mutate(index = seq(n())) %>% 
  group_by(Mark) %>% 
  summarize(first_index = min(index)))$first_index

## determine the first period in which each individual was present
first_capt <- capt_history.p %>% 
  group_by(Mark, year, Site) %>% 
  summarize(capt = sum(captured)) %>% 
  mutate(capt = cumsum(capt)) %>% 
  mutate(capt = ifelse(capt > 0, 1, 0)) 

## time periods in which we do not know if an individual was present or not
for (k in 1:n_sites) {
p_zeros <- matrix(data = 0, nrow = n_ind.per[k, 1], ncol = sum(n_occ[, u_sites[k]]))
for (i in 1:n_ind.per[k, 1]) {
  tdat <- first_capt %>% filter(Site == u_sites[k])
  tdat %<>% filter(Mark == unique(tdat$Mark)[i])
  rep.t <- n_occ[, u_sites[k]]; rep.t <- rep.t[which(rep.t != 0)]
  p_zeros[i, ] <- rep(tdat$capt, rep.t)
  p_zeros[i, ] <- ifelse(cumsum(p_zeros[i, ]) > 0, 1, 0)
}
p_zeros.t   <- (p_zeros %>% reshape2::melt() %>% arrange(Var1))$value

if (k == 1) {
p_zeros.a <- p_zeros.t
} else {
p_zeros.a <- c(p_zeros.a, p_zeros.t)
}
}

capt_history.p$p_zeros <- p_zeros.a

## Vectors for survival
last_week        <- capt_history %>% group_by(Site) %>% 
  filter(sampled == 1) %>% summarize(last_week = max(week_year))
capt_history.phi <- capt_history %>% 
  left_join(., last_week) %>%
  filter(week_year != last_week, sampled == 1) %>% ungroup()

## Determine the number of time periods that elapse between back to back samples
the_weeks.a <- capt_history %>% 
  filter(sampled == 1) %>% 
  group_by(Site) %>%
  summarize(cont_weeks = unique(cont_weeks))

for (k in 1:n_sites) {
  the_weeks.t <- (the_weeks.a %>% filter(Site == u_sites[k]))$cont_weeks
  time_gaps.t <- (the_weeks.t - lag(the_weeks.t, 1))[-1]
  
  time_gaps.ta <- rep(time_gaps.t, n_ind.per[k, 1])
  
  if (k == 1) {
   time_gaps   <- time_gaps.t
   time_gaps.a <- time_gaps.ta
  } else {
   time_gaps   <- c(time_gaps, time_gaps.t)
   time_gaps.a <- c(time_gaps.a, time_gaps.ta)
  }
}

## Weeks between sampling events
capt_history.phi %<>% mutate(time_gaps = time_gaps.a)

## Offseason vector (not the best strategy, but ok)
capt_history.phi %<>% mutate(offseason = ifelse(time_gaps > 20, 1, 0))

phi_first_index <- (capt_history.phi %>% mutate(index = seq(n())) %>% 
    group_by(Mark) %>% 
    summarize(first_index = min(index)))$first_index

## Indices for which entries of phi must be 0
for (k in 1:n_sites) {
phi_zeros <- matrix(data = 0, nrow = n_ind.per[k, 1], ncol = sum(n_occ[, u_sites[k]]) - 1)

for (i in 1:n_ind.per[k, 1]) {
  tdat <- first_capt %>% filter(Site == u_sites[k])
  tdat %<>% filter(Mark == unique(tdat$Mark)[i])
  this_ind <- tdat$Mark[1]

  phi_zeros[i, ] <- c(
    rep(1, capture_range$first[this_ind] - 1)
  , rep(0, ncol(phi_zeros) - (capture_range$first[this_ind] - 1))
  )
  
}

phi_zeros.t   <- (phi_zeros %>% reshape2::melt() %>% arrange(Var1))$value

if (k == 1) {
phi_zeros.a <- phi_zeros.t
} else {
phi_zeros.a <- c(phi_zeros.a, phi_zeros.t)
}
}

capt_history.phi$phi_zeros <- phi_zeros.a

####
## Other needed covaraites 
####

## OCT 20: Hmm, not good. Temp will need to be imputed?
temp <- expand.grid(
  year = unique(A11$year)
, week = seq(week_range[1], week_range[2])
, Site = unique(A11$Site)
)

temp_have <- A11 %>% 
  group_by(year, week, Site) %>% 
  summarize(temp = mean(Site_temp, na.rm = T))

temp %<>% left_join(., temp_have)

temp.lm <- lm(
  temp ~ week
, data = temp
)

predvals <- data.frame(
  week     = sort(unique(temp$week))
, predvals = predict(temp.lm, newdata = data.frame(week = sort(unique(temp$week))))
  )

temp %<>% left_join(., predvals) %>% mutate(temp = ifelse(is.na(temp), predvals, temp))

temp <- temp %>% arrange(year, week, Site) %>% 
  dplyr::select(-predvals) %>% pivot_wider(names_from = Site, values_from = temp) %>%
  dplyr::select(-year, -week) 
temp <- as.matrix(temp)

####
## Run the stan model
####

stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop           = n_sites
 , n_periods       = n_periods
 , n_ind           = n_ind
 , n_times         = n_times.a
 , times_within    = n_times.w
 , ind_occ         = sum(colSums(n_occ) * c(n_ind.per))
 , ind_occ_min1    = sum((colSums(n_occ) - 1) * c(n_ind.per))
  
  ## short vector indexes 
 , time              = rep(seq(n_times.w), n_periods)
 , time_per_period   = matrix(data = seq(n_times.w * n_periods), nrow = n_times.w, ncol = n_periods)
 , periods           = rep(seq(n_periods), each = n_times.w)
 , ind_occ_size      = rep(colSums(n_occ), n_ind.per)
 , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)
 , ind_in_pop        = rep(1, n_ind)

 , phi_first_index   = phi_first_index
 , p_first_index     = p_first_index
  
  ## long vector indexes
 , ind_occ_min1_rep    = capt_history.phi$Mark
 , sampling_events_phi = capt_history.phi$week_year
 , offseason           = capt_history.phi$offseason
 , pop_phi             = as.numeric(as.factor(capt_history.phi$Site))
 , phi_zeros           = capt_history.phi$phi_zeros

 , ind_occ_rep       = capt_history.p$Mark
 , sampling_events_p = capt_history.p$week
 , periods_occ       = as.numeric(as.factor(capt_history.p$year))
 , pop_p             = as.numeric(as.factor(capt_history.p$Site))
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
  file    = "CMR_simulation/CMR_empirical_pr.stan"
, data    = stan_data
, chains  = 1
, refresh = 20
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )
  
# saveRDS(stan.fit, "stan.fit.empirical.Rds")
# stan.fit <- readRDS("stan.fit.empirical.Rds")

shinystan::launch_shinystan(stan.fit)

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

bd_levels <- log(c(seq(1, 8, by = 0.5) %o% 10^(0:5)))

####
## Recovery of simulated coefficients?
##   *NOTE*: High success with keeping all individuals, only moderate success with dropping individuals
####

## Primary Bd effects
pred_coef        <- as.data.frame(stan.fit.summary[
    c(
      grep("beta_bd", dimnames(stan.fit.summary)[[1]])
    , grep("beta_p", dimnames(stan.fit.summary)[[1]])
    )
  , c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(param = rownames(.))

pred_coef %>% {
  ggplot(., aes(param, mid)) + geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    xlab("Parameter") + ylab("Estimate") +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.5) +
    theme(axis.text.x = element_text(size = 11))
}

## Other Bd effects
pred_coef        <- as.data.frame(stan.fit.summary[
    c(
      grep("beta_offseason", dimnames(stan.fit.summary)[[1]])
    , grep("beta_timegaps", dimnames(stan.fit.summary)[[1]])
    )
  , c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(param = rownames(.))

pred_coef %>% {
  ggplot(., aes(param, mid)) + geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    xlab("Parameter") + ylab("Estimate") +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.5) +
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
    xlab("Log of Bd Load") + 
    ylab("Predicted detection probability")
}

between_seas <- stan.fit.samples$phi[, (phi_first_index[-c(1:19)] + 8 + 12)]
between_seas <- reshape2::melt(between_seas)
between_seas %<>% filter(value != 0)

between_seas %>% {
  ggplot(., aes(x = value)) + geom_histogram(aes(y = stat(count) / sum(count)), bins = 500) +
    xlab("Between Year Survival") +
    ylab("Frequency")
}

####
## What does predicted overall bd look like?
####
stan.pred     <- matrix(nrow = length(unique(capt_history$week)), ncol = 1000, data = 0)
stan.pred.ind <- array(dim = c(length(unique(capt_history$Mark)), length(unique(capt_history$week)), 1000), data = 0)
samp_occ      <- seq(length(unique(capt_history$week)))

for (i in 1:nrow(stan.pred)) {
stan.pred[i, ] <- stan.fit.samples$beta_bd[, 1] + 
  stan.fit.samples$beta_bd[, 2] * samp_occ[i] +
  stan.fit.samples$beta_bd[, 3] * samp_occ[i]^2 
  stan.fit.samples$beta_bd[, 4] * temp[i, 1]

for (j in 1:dim(stan.pred.ind)[1]) {
stan.pred.ind[j,i,] <- stan.fit.samples$bd_ind[, j] + 
  stan.fit.samples$beta_bd[, 1] +
  stan.fit.samples$beta_bd[, 2] * samp_occ[i] +
  stan.fit.samples$beta_bd[, 3] * samp_occ[i]^2 
  stan.fit.samples$beta_bd[, 4] * temp[i, 1]
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

stan.pred %<>% mutate(occ = plyr::mapvalues(occ, from = occ
  , to = unique(capt_history$week)
    ))

stan.pred.ind <- reshape2::melt(stan.pred.ind)
names(stan.pred.ind) <- c("ind", "occ", "iter", "value")

stan.pred.ind %<>% 
  group_by(occ, ind) %>%
  summarize(
    lwr = quantile(value, 0.025)
  , mid = quantile(value, 0.50)
  , upr = quantile(value, 0.975)
  )

stan.pred.ind %<>% mutate(occ = plyr::mapvalues(occ
  , from = unique(stan.pred.ind$occ)
  , to = unique(capt_history$week)
    ))

ggplot(stan.pred, aes(occ, mid)) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line() +
  geom_line(data = stan.pred.ind
    , aes(occ, mid, group = ind), lwd = 0.5, alpha = 0.2) +
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
