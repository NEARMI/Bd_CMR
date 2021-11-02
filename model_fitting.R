########################################
## Fit a CMR model with the newt data ##
########################################

####
## Notes as of OCT 28:
####

## 1) Multi-population fit went ok, but with some divergent transitions.
 ## There are definitely some weird populations, for example, A05, in which every individual was only captured a single time each

### Next steps:
 ## 1) Add back the summary of bd informing between season survival for the long-form model
 ## 2) Step back to a single population and start adding other covariates
   ## -- Other individual-level covaraites
   ## -- more complicated detection model
 ## 3) Play with decreasing the number of sampling periods and try just estimating between season survival as 
  ## a function of average estimated bd load

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
red_ind        <- TRUE    ## TRUE   = reduce the number of individuals for debugging purposes
if (red_ind) {
  red_ind.count <- 50
}

## Quick look at some of the data
Bd_Newts_AllSites %>% group_by(Site) %>% summarize(n_y = length(unique(Mark))) %>% {
  ggplot(., aes(x = n_y)) + 
    geom_histogram(bins = nrow(.)) +
    xlab("Number of Individuals") +
    ylab("Number of Populations")
}

## Select a subset of sites
A11 <- Bd_Newts_AllSites %>% 
  filter(SA == "PA") %>% 
#  filter(Site == "A11") %>%
# filter(Site == "A04" | Site == "A11" | Site == "A05" | Site == "A10") %>%
  group_by(year) %>%  
  mutate(week = ceiling(julian / 7)) %>% 
  filter(!is.na(Site)) %>% 
  arrange(Site)

n_sites <- length(unique(A11$Site))
u_sites <- unique(A11$Site)

## Find the first and last week and year ever sampled in each population
week_range <- A11 %>% 
  group_by(Site) %>% 
  summarize(
    min_week = min(week)
  , max_week = max(week)
  )

year_range <- A11 %>% 
  group_by(Site) %>% 
  summarize(
    n_years = length(unique(year))
  )

## Find the unique weeks sampled in each year
sampled_weeks <- A11 %>% 
  group_by(Site, year) %>%
  summarize(week = unique(week)) %>%
  mutate(sampled = 1) %>%
  arrange(Site, week) 

## Construct an "all possible combinations" data frame (all of the possible windows in time in which an individual
 ## could have been captured) and then parse it down
all_ind <- 0
for (i in 1:n_sites) {
  
## Extract a given site and sampling characteristics of that site
week_range.i    <- week_range %>% filter(Site == u_sites[i])
A11.i           <- A11 %>% filter(Site == u_sites[i])
sampled_weeks.i <- sampled_weeks %>% filter(Site == u_sites[i])
  
capt_history.t <- 
  ## First create that "all possible combinations" data frame
  expand.grid(
  week = seq(from = week_range.i$min_week, to = week_range.i$max_week, by = 1)
, year = unique(A11.i$year)
, Site = u_sites[i]
, Mark = unique(A11.i$Mark)) %>% 
  ## Add in which weeks were sampled and which individuals were sampled
  left_join(., sampled_weeks.i) %>% 
  left_join(., (A11.i %>% dplyr::select(week, year,  Mark, month, copies.swab))) %>% 
  ## left_join figures out which days were sampled and which were not
  mutate(sampled = ifelse(is.na(sampled), 0, 1)) %>%
  rename(
  ## just using a random non-na column to find captures (convenient given how left-join works)
    captured = month
  , bd_load  = copies.swab) %>% 
  ## convert other nas to 0s
  mutate(
    captured = ifelse(is.na(captured), 0, 1)
  , swabbed  = ifelse(is.na(bd_load), 0, 1)) %>%
  ## collapsing to week (there are a _very_ few number of individuals captured multiple times in the same week,
    ## so this will have an extremely negligable effect)
  group_by(Mark, week, year, Site) %>%
  summarize(
    sampled  = sum(sampled)
  , captured = sum(captured, na.rm = T)
  , swabbed  = sum(swabbed , na.rm = T)
  , bd_load  = sum(bd_load , na.rm = T)
  ) %>% 
  mutate(
    sampled     = ifelse(sampled  > 1, 1, sampled)
  , captured    = ifelse(captured > 1, 1, captured)
  , swabbed     = ifelse(swabbed  > 1, 1, swabbed)
  ## Need to do better here
  , log_bd_load = log(bd_load + 1)                           ### eeek!
  , log_bd_load = ifelse(is.na(log_bd_load), 0, log_bd_load) ### eeek X2!!
  )

## Jump through a few hoops to name unique individuals. 
 ## NOTE: this is an issue if individuals move populations
capt_history.t %<>% mutate(Mark = as.factor(Mark)) %>%
  mutate(Mark = as.numeric(Mark))
n_inds <- max(capt_history.t$Mark)
capt_history.t %<>% mutate(Mark = Mark + all_ind)
all_ind <- all_ind + n_inds

## Add in other needed indexing columns
capt_history.t %<>% ungroup() %>% 
  arrange(year, week, Mark, Site) %>% 
  mutate(
## weeks counted from the first week in each year
   week_year  = interaction(week, year)
 , week_year  = as.factor(week_year)
 , week_year  = as.numeric(week_year)
## calculating continuous weeks from the first capture opportunity onward
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
 
## Sort data frame in the appropriate order (counting through consecutive weeks one individual at a time)
capt_history %<>% arrange(Mark, year, Site, week)

## individuals' measured bd 
capt_history.bd_load <- capt_history %>% 
  ungroup() %>%
  arrange(Mark, week_year) %>%
  filter(swabbed == 1)

## first and last _OF THE CAPTURE EVENTS_ in which each individual was captured
 ## (possible min and max will vary by which population individuals are in)
capture_range  <- capt_history %>% 
  group_by(Mark) %>% 
  filter(sampled == 1) %>%  
  summarize(
    first = min(which(captured == 1))
  , final = max(which(captured == 1))) %>% 
  dplyr::select(first, final) %>% 
  ## Remove all individuals in the data set that were never captured (in case there are any for w/e reason)
  filter(!is.infinite(first) | !is.infinite(final))

## Plot the data (as long as there is a reasonable number of individuals)
if (length(unique(capt_history$Mark)) < 300) {
capt_history %>% 
  mutate(event = week) %>% {
  ggplot(., aes(week, Mark, fill = as.factor(captured))) + 
    geom_tile(aes(alpha = sampled)) +
    geom_point(data = 
        capt_history %>% 
        mutate(event = week) %>% 
        filter(swabbed == 1)
      , aes(x = week, y = Mark, z = NULL), lwd = 0.7) +
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
}
  
## For debugging try a reduced number of individuals
if (red_ind) {

## Select a random subset of individuals 
test_ind <- capt_history %>% group_by(Site) %>%
  summarize(n_ind = length(unique(Mark))) %>% 
  dplyr::select(-Site) %>% as.matrix()
test_ind <- ifelse(test_ind > red_ind.count, red_ind.count, test_ind)

## pull out the appropriate individuals from each location
for (i in 1:n_sites) {
  temp_dat  <- capt_history %>% filter(Site == u_sites[i])
  # rand_inds <- (capt_history %>% group_by(Mark) %>% summarize(num_capt = sum(captured)) %>% arrange(desc(num_capt)) %>% filter(num_capt > 3))$Mark
  rand_inds <- sort(sample(unique(temp_dat$Mark), test_ind[i, 1]))
  temp_dat  <- temp_dat %>% filter(Mark %in% rand_inds)
  if (i == 1) {
    temp_dat.a <- temp_dat
  } else {
   temp_dat.a <- rbind(temp_dat.a, temp_dat) 
  }
}

## and adjust the names of the individuals
rand_inds     <- unique(temp_dat.a$Mark)
capt_history  <- temp_dat.a %>% mutate(Mark = as.factor(Mark)) %>% mutate(Mark = as.numeric(Mark))
capture_range <- capture_range[rand_inds, ]
capt_history.bd_load <- capt_history %>% 
  ungroup() %>%
  arrange(Mark, week_year) %>%
  filter(swabbed == 1)

}

####
## Data in the needed structure for the stan model
####

## total number of individuals
n_ind     <- length(unique(capt_history$Mark))             

## individuals per population
n_ind.per <- capt_history %>% group_by(Site) %>%
  summarize(n_ind = length(unique(Mark))) %>% 
  dplyr::select(-Site) %>% as.matrix()

## samoling occasions per population per year
n_occ     <- sampled_weeks %>% group_by(year, Site) %>%
  summarize(n_occ = length(unique(week))) %>% 
  mutate(Site = factor(Site, levels = u_sites)) %>%
  arrange(Site) %>%
  pivot_wider(values_from = n_occ, names_from = Site) %>% 
  arrange(year) %>%
  ungroup() %>%
  dplyr::select(-year) %>% as.matrix()
n_occ[is.na(n_occ)] <- 0

###
## Note: The way the "long-form / database-form" model works is to have long vectors of the
## data and outcomes and index vectors giving details/grouping associations about each data point
## Given the structure of the stan model, the easiest way to set up the correct structure is to 
## subest the complete data frame into three that are of the appropraite length and go from there
###

### --- Data for detection (.p for detection) --- ###
capt_history.p   <- capt_history %>% 
  filter(sampled == 1) %>% 
  ungroup()

## Which entries of p correspond to a new individual (the first entry for each individual)
p_first_index <- (capt_history.p %>% mutate(index = seq(n())) %>% 
  group_by(Mark) %>% 
  summarize(first_index = min(index)))$first_index

## determine the first period (for now year) in which each individual was _known_ to be present
first_capt <- capt_history.p %>% 
  group_by(Mark, year, Site) %>% 
  summarize(capt = sum(captured)) %>% 
  mutate(capt = cumsum(capt)) %>% 
## And then in all future times from the current time these individuals _could_ be here
  mutate(capt = ifelse(capt > 0, 1, 0)) 

## For each individual extract which time periods we do not know if an individual was present or not
for (k in 1:n_sites) {
p_zeros <- matrix(data = 0, nrow = n_ind.per[k, 1], ncol = sum(n_occ[, u_sites[k]]))
for (i in 1:n_ind.per[k, 1]) {
  tdat <- first_capt %>% filter(Site == u_sites[k])
  tdat %<>% filter(Mark == unique(tdat$Mark)[i])
  ## For each individual repeat 0s and 1s for each sampling occasions in all years that they were
   ## never captured (0s) and captured (1s) or after a first capture year (1s) 
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

## Add the p_zeros to the detection data frame for easier debugging
capt_history.p$p_zeros <- p_zeros.a

## These p_zeros are used to inform a scaling factor on detection probability (one scaling factor for each
 ## individual in each year). Need an index for these scaling factors
capt_history.p %<>% 
  mutate(gamma_index = paste(interaction(Mark, year),"a",sep="_")) %>% 
  mutate(gamma_index = factor(gamma_index, levels = unique(gamma_index))) %>%
  mutate(gamma_index = as.numeric(gamma_index))
  
### --- Data for survival (.phi for survival) --- ###

## phi not calculable on the last time step so drop it
last_week        <- capt_history %>% group_by(Site) %>% 
  filter(sampled == 1) %>% summarize(last_week = max(week_year))

capt_history.phi <- capt_history %>% 
  left_join(., last_week) %>%
  filter(week_year != last_week, sampled == 1) %>% ungroup()

## Determine the number of time periods that elapse between back to back samples.
 ## Do this with the data without the last date dropped (as phi is survival to the next)
  ## and then add to capt_history.phi
time_gaps <- (capt_history %>% 
  filter(sampled == 1) %>%
  group_by(Site, Mark) %>% 
  mutate(
   time_gaps =  cont_weeks - lag(cont_weeks, 1)
  ) %>% filter(!is.na(time_gaps)))$time_gaps
 
capt_history.phi %<>% mutate(time_gaps = time_gaps)

## Offseason vector 
capt_history.phi %<>% 
  group_by(Site) %>%
## Tricky/cool way to find the number of years -1 largest time gaps --- (should always work???)
  mutate(offseason = ifelse(time_gaps >= sort(unique(time_gaps), decreasing = TRUE)[length(unique(year)) - 1]
  , 1, 0)) %>% as.data.frame()

## Which entries of phi correspond to a new individual (the first entry for each individual)
phi_first_index <- (capt_history.phi %>% mutate(index = seq(n())) %>% 
    group_by(Mark) %>% 
    summarize(first_index = min(index)))$first_index

## Indices for which entries of phi must be 0. See above notes for p_zeros. Similar idea here, but 
 ## here the differentiation for 0s and 1s are for prior to and after an individual was captured for the first time
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

## Also add this one to the data frame for ease of debugging
capt_history.phi$phi_zeros <- phi_zeros.a

## Index for 
capt_history.phi %<>% 
  mutate(X_stat_index = paste(interaction(Mark, year),"a",sep="_")) %>% 
  mutate(X_stat_index = factor(X_stat_index, levels = unique(X_stat_index))) %>%
  mutate(X_stat_index = as.numeric(X_stat_index))

### --- Data for latent bd --- ###

## Latent bd is estimated over the whole time period and not just for the capture occasions,
 ## though bd on the capture occasions are used to determine detection and survival. Need to
  ## determine what entries of phi, and p correspond to the full time period bd. This is done here

## for calculating summaries of latent bd for between season survival
bd_first_index <- (capt_history %>% mutate(index = seq(n())) %>% 
  group_by(Mark, year, Site) %>% 
  summarize(first_index = min(index)))$first_index
bd_last_index  <- (capt_history %>% mutate(index = seq(n())) %>% 
  group_by(Mark, year, Site) %>% 
  summarize(last_index = max(index)))$last_index
  
## Index for every entry of bd (all time points)
capt_history %<>% mutate(index = seq(n()))

## Which of all of the time entries correspond to the correct phi entries 
phi_bd_index <- (left_join(
  capt_history.phi %>% dplyr::select(Mark, week, year, Site)
, capt_history     %>% dplyr::select(Mark, week, year, Site, index)
  ))$index

## to put it another way, the phi_bd_index of the latent bd vector is the bd associated with the nth row
 ## of capt_history.phi
capt_history.phi %<>% mutate(phi_bd_index = phi_bd_index)

## same thing for p
p_bd_index <- (left_join(
  capt_history.p %>% dplyr::select(Mark, week, year, Site)
, capt_history   %>% dplyr::select(Mark, week, year, Site, index)
  ))$index

capt_history.p %<>% mutate(p_bd_index = p_bd_index)

## And finally, what actual measured bd values inform the latent bd process?
x_bd_index <- (left_join(
  capt_history.bd_load %>% dplyr::select(Mark, week, year, Site)
, capt_history         %>% dplyr::select(Mark, week, year, Site, index)
  ))$index

capt_history.bd_load %<>% mutate(x_bd_index = x_bd_index)

####
## Finally, deal with any other needed covariates 
####

temp <- A11 %>% 
  group_by(year, week, Site) %>% 
  summarize(temp = mean(Site_temp, na.rm = T))

temp.lm <- lm(
  temp ~ week
, data = temp
)

temp.need <- capt_history %>% 
  group_by(week, year, Site) %>% 
  summarize(num_mark = length(unique(Mark))) %>% dplyr::select(-num_mark)

predvals <- data.frame(
  week     = sort(unique(temp.need$week))
, predvals = predict(temp.lm, newdata = data.frame(week = sort(unique(temp.need$week))))
  )

temp.need <- left_join(temp.need, temp) %>% left_join(., predvals) %>% 
  mutate(temp = ifelse(is.na(temp), predvals, temp)) %>% 
  dplyr::select(-predvals)

capt_history %<>% left_join(., temp.need)

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
 , n_ind           = n_ind
 , ind_per_period  = sum(year_range$n_years * n_ind.per)
  
 , ind_time        = sum((week_range[, 3] - week_range[, 2] + 1) * year_range$n_years * n_ind.per)
 , ind_occ         = sum(colSums(n_occ) * c(n_ind.per))
 , ind_occ_min1    = sum((colSums(n_occ) - 1) * c(n_ind.per))
  
  ## short vector indexes 
 , ind_occ_size      = rep(colSums(n_occ), n_ind.per)
 , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)

 , p_first_index     = p_first_index
 , phi_first_index   = phi_first_index
  
  ## long vector indexes
 , ind_occ_rep       = capt_history.p$Mark
 , sampling_events_p = capt_history.p$week
 , periods_occ       = as.numeric(as.factor(capt_history.p$year))
 , pop_p             = as.numeric(as.factor(capt_history.p$Site))
 , p_zeros           = capt_history.p$p_zeros
 , p_bd_index        = capt_history.p$p_bd_index
 , gamma_index       = capt_history.p$gamma_index
  
 , ind_occ_min1_rep    = capt_history.phi$Mark
 , sampling_events_phi = capt_history.phi$week_year
 , offseason           = capt_history.phi$offseason
 , pop_phi             = as.numeric(as.factor(capt_history.phi$Site))
 , phi_zeros           = capt_history.phi$phi_zeros
 , phi_bd_index        = capt_history.phi$phi_bd_index
 , X_stat_index        = capt_history.phi$X_stat_index

 , ind_bd_rep          = capt_history$Mark
 , sampling_events_bd  = capt_history$week
 , ind_in_pop          = as.numeric(as.factor(capt_history$Site))
 , temp                = capt_history$temp

  ## covariates
 , N_bd            = nrow(capt_history.bd_load)
 , X_bd            = capt_history.bd_load$log_bd_load  
 , x_bd_index      = capt_history.bd_load$x_bd_index
 , bd_first_index  = bd_first_index
 , bd_last_index   = bd_last_index
 , time_gaps       = capt_history.phi$time_gaps
  
  ## Capture data
 , N_y             = nrow(capt_history.p)
 , y               = capt_history.p$captured
  
 , first           = capture_range$first
 , last            = capture_range$final

  )

stan.fit  <- stan(
  file    = "CMR_empirical_long.stan"
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
# stan.fit.samples <- readRDS("CMR_simulation/phi.samples.Rds")

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

capt_history %<>% mutate(est_bd = colMeans(stan.fit.samples$X))

capt_history %>% filter(year == 2019) %>% {
  ggplot(., aes(cont_weeks, est_bd)) + geom_line(aes(group = Mark)) +
    geom_line(
      data = capt_history %>% filter(year == 2019) %>% filter(swabbed == 1)
     , aes(cont_weeks, log_bd_load, group = Mark), colour = "red"
      )
}

capt_history.phi %<>% mutate(phi.est = colMeans(stan.fit.samples$phi))

capt_history.phi %>% 
  filter(offseason == 1) %>% 
  filter(phi_zeros == 0) %>% {
  ggplot(., aes(x = phi.est)) +
      geom_histogram(bins = 50)
}


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
