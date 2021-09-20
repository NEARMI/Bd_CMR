################################################
## Try and fit a CMR model with the newt data ##
################################################

## First, see and start with 'data_exploration.R' for data info
## Second, see 'individual_CRM.R' for a CMR model with a bd submodel that fits reasonablly well with simulated data

## There is basically no chance the model is actually going to work, but just for kicks (and to 
 ## set the data up in the right structure to run the model)

capt_history %<>% mutate(log_bd_load = ifelse(is.na(log_bd_load), 0, log_bd_load))

## Try a small subset to see if the model will run (regardless of whether it will fit or not)
capt_history.well_meas <- capt_history %>% 
  group_by(Mark) %>%
  summarize(total_capt = sum(captured)) %>% 
  filter(total_capt > 4) %>% 
  droplevels()

capt_history %<>% filter(Mark %in% capt_history.well_meas$Mark)

capt_history.matrix <- matrix(
  nrow = length(unique(capt_history$Mark))
, ncol = length(unique(capt_history$Date))
, data = capt_history$captured
, byrow = T)

capt_history.bd_load <- matrix(
  nrow = length(unique(capt_history$Mark))
, ncol = length(unique(capt_history$Date))
, data = capt_history$log_bd_load
, byrow = T)

capt_history.bd_measured <- matrix(
  nrow = length(unique(capt_history$Mark))
, ncol = length(unique(capt_history$Date))
, data = capt_history$swabbed
, byrow = T)

capture_range  <- capt_history %>% group_by(Mark) %>% 
  summarize(
    first = min(which(captured == 1))
  , final = max(which(captured == 1))) %>% 
  dplyr::select(first, final) %>% 
  mutate(
    first = ifelse(is.infinite(first), 0, first)
  , final = ifelse(is.infinite(final), 0, final)
    )

capture_total <- capt_history %>% group_by(Date) %>% summarize(total_capt = sum(captured))

capt_history %>% mutate(event = as.factor(Date)) %>% {
  ggplot(., aes(event, Mark, fill = as.factor(captured))) + 
    geom_tile(alpha = 0.8) +
    xlab("Sampling Event") + ylab("Individual") +
    scale_fill_manual(
        values = c("dodgerblue4", "firebrick4")
      , name   = "Detected?"
      , labels = c("No", "Yes")) +
    theme(
      axis.text.y = element_text(size = 6)
    , axis.text.x = element_text(size = 8, angle = 300)
    , legend.text = element_text(size = 12)
    , legend.key.size = unit(.55, "cm")
    ) 
}

stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin

stan_data     <- list(
  ## bookkeeping params
   n_ind         = length(unique(capt_history$Mark))
 , n_occasions   = length(unique(capt_history$Date))
  ## Capture data
 , y             = capt_history.matrix
 , first         = capture_range$first
 , last          = capture_range$final
 , n_captured    = capture_total$total_capt
  ## Covariate associated parameters
 , X_bd          = capt_history.bd_load
 , X_measured    = capt_history.bd_measured
  )

stan.fit  <- stan(
  file    = "CMR_mine/CMR_ind_some_bd_empirical.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

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
    scale_x_discrete(labels = c("Detection intercept", "Detection slope", "Survival intercept", "Survival slope")) +
    theme(axis.text.x = element_text(size = 11))
}

## survival over Bd load
stan.pred        <- apply(stan.fit.samples$beta_phi, 1
  , FUN = function(x) plogis(x[1] + x[2] * bd_levels)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "mortality")
stan.pred        %<>% mutate(log_bd_load = plyr::mapvalues(log_bd_load
  , from = unique(log_bd_load), to = log(bd_levels)))
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
  , from = unique(log_bd_load), to = log(bd_levels)))
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
## Individual random effect estimates
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
