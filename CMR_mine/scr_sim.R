####
## Explore simulation parameters in a convenient way
####

bd_beta   <- c(10, 30, -50) ## Intercept and slope for mean response
                            ##  ** SEP 23: This is pretty un-dynamic and probably needs to be updated.
                             ## for now, based on the patterns in the real data estimate bd based on some quadratic function of time
bd_sigma  <- 2              ## observation noise
bd_theta  <- c(2)           ## random effect variance covariance

## Response of individuals to Bd load
bd_mort   <- c(decay = -0.45, offset = 6)     ## logistic response coefficients for mortality across log(bd_load)
bd_detect <- c(decay = 0.1, offset = -0.5)     ## logistic response coefficients for detection across log(bd_load)

temp_beta <- c(5)

expdat <- expand.grid(
  times = seq(times)
, ind   = factor(seq(ind))
  )
expdat %<>% mutate(
  bd_load   = simulate(~poly(times, 2) + (1 | ind)
  , nsim    = nsim
  , family  = Gamma(link = "log")
  , newdata = expdat
  , newparams = list(
      beta  = bd_beta
    , sigma = bd_sigma
    , theta = bd_theta
    )
  )$sim_1
) %>% mutate(
  bd_load     = round(bd_load, digits = 0)
, log_bd_load = round(log(bd_load), digits = 0)
) %>% mutate(
  log_bd_load = ifelse(is.infinite(log_bd_load), 0, log_bd_load)
)

expdat %>% {
  ggplot(., aes(times, bd_load)) + 
    geom_line(aes(group = ind)) +
    scale_y_log10()
}

## bd effects simulation
bd_mort   <- c(decay = -0.3, offset = 6)     ## logistic response coefficients for mortality across log(bd_load)
bd_detect <- c(decay = 0.1, offset = -0.5)     ## logistic response coefficients for detection across log(bd_load)
  
## establish "true" relationship between load and mortality probability and detection probability
 ## Note: here "mort" is actually survival probability from time t to t+1
bd_probs <- data.frame(
  log_bd_load  = bd_range
, mort         = exp(bd_mort["offset"] + bd_mort["decay"] * bd_range) /
    (1 + exp(bd_mort["offset"] + bd_mort["decay"] * bd_range))
) %>% left_join(.
  , data.frame(
  log_bd_load  = bd_range
, detect       = exp(bd_detect["offset"] + bd_detect["decay"] * bd_range) /
    (1 + exp(bd_detect["offset"] + bd_detect["decay"] * bd_range))
)
  )

bd_probs %>% pivot_longer(cols = c(2, 3)) %>% {
  ggplot(., aes(log_bd_load, value)) +
  geom_line(aes(colour = name)) + 
  scale_colour_manual(
      name   = "Relationship"
    , values = c("firebrick3", "dodgerblue3")
    , labels = c("Detection", "Survival")) + 
  ylab("Prediction") + 
  xlab("log of Bd load") +
  theme(legend.key.size = unit(0.75, "cm"))
}

expdat <- expand.grid(
  periods = seq(periods)
, times   = seq(times[1])
, ind     = factor(seq(ind))
  )

expdat %<>% mutate(
  bd_load   = simulate(~periods + poly(times, 2) + (1 | ind)
  , nsim    = nsim
  , family  = Gamma(link = "log")
  , newdata = expdat
  , newparams = list(
      beta  = c(8, 2, 30, -50)
    , sigma = bd_sigma
    , theta = bd_theta
    )
  )$sim_1
) %>% mutate(
  bd_load     = round(bd_load, digits = 0)
, log_bd_load = round(log(bd_load), digits = 0)
) %>% mutate(
  log_bd_load = ifelse(is.infinite(log_bd_load), 0, log_bd_load)
)

expdat %>% {
  ggplot(., aes(times, bd_load)) + 
    geom_line(aes(group = ind)) +
    scale_y_log10() +
    facet_wrap(~periods)
}
