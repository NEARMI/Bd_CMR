##########################################################
## Try the simple spotted frog model with the newt data ##
##########################################################

####
## Notes as of NOV 15:
####

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
# filter(SA == "PA") %>% 
  filter(Site == "A11") %>%
# filter(Site == "A04" | Site == "A11" | Site == "A05" | Site == "A10") %>%
  group_by(year) %>%  
  mutate(week = ceiling(julian / 7)) %>% 
  filter(!is.na(Site)) %>% 
  arrange(Site)

## Other covariates to include for populations
# Area
# Density

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

## Before converting Mark to a numeric, find the individual specific covariates to be used later
 ## Other covariates to include for individuals  
  # SVL
  # Sex
ind_cov <- A11.i %>% group_by(Mark, year, week) %>% 
  summarize(
    size = mean(SVL, na.rm = T)  ## returns the one value or a mean if caught multiple times in one week
  , sex  = Sex)

## -- Tidy up messy data -- ##
## Sex
uni_ind <- unique(ind_cov$Mark)
for (z in seq_along(uni_ind)) {
  
## Need to come back to this (NOV 3)
temp.ind_cor <- ind_cov %>% 
  filter(Mark == uni_ind[z]) %>% 
  group_by(Mark, year, week) %>%
  summarize(size = mean(size), sex = sex[1]) 

if (length(unique(temp.ind_cor$sex)) > 1) {
  tt        <- table(temp.ind_cor$sex, useNA = c("ifany"))
  real_sex  <- names(tt[which(tt == max(tt))])
  
  num_entry <- length(real_sex)
  an_na     <- any(is.na(real_sex))
  
  if ((num_entry == 2) & an_na) {
    temp.ind_cor$sex <- real_sex[which(!is.na(real_sex))]
  } else if ((num_entry == 2) & !an_na) {
    temp.ind_cor$sex <- "U"
  } else if ((num_entry == 1) & an_na) {
    temp.ind_cor$sex <- "U"
  } else if ((num_entry == 1) & !an_na) {
    temp.ind_cor$sex <- real_sex
  } else if (num_entry == 3) {
    temp.ind_cor$sex <- "U"
  }
  
} else {
  if (is.na(temp.ind_cor$sex)) {
    temp.ind_cor$sex <- "U"
  }
}

if (z == 1) {
  temp.ind_cor.a <- temp.ind_cor
} else {
  temp.ind_cor.a <- rbind(temp.ind_cor.a, temp.ind_cor)
}

}

## Size
temp.ind_cor.a %<>% group_by(Mark, year) %>% 
  mutate(size = mean(size, na.rm = T))
 ## NOV 3: for now (TO CHANGE to multiple imputation eventually)
need_size <- which(is.na(temp.ind_cor.a$size))
for (z in seq_along(need_size)) {
 temp.ind_cor.a[need_size[z], ]$size <- mean(temp.ind_cor.a[(temp.ind_cor.a$year == temp.ind_cor.a$year) & (temp.ind_cor.a$sex == temp.ind_cor.a$sex), ]$size, na.rm = T)   
}

capt_history.t %<>% left_join(., temp.ind_cor.a)

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

## Find clusters of primary periods
find_clusters <- capt_history %>% filter(Mark == 1)

for (i in 1:length(unique(find_clusters$year))) {
  fc.t <- find_clusters %>% filter(year == unique(find_clusters$year)[i])

  week.t <- fc.t[which(fc.t$sampled == 1), ]$week 
  rep.we <- which((week.t - lag(week.t, 1))[-1] == 1)
  
  clust.t <- numeric(length(week.t))
  
  pp <- 1
  for (j in 1:length(clust.t)){
    clust.t[j] <- pp
    if (j %notin% rep.we) {
      pp <- pp + 1
    }
  }
  
  samp.clust <- data.frame(
    week       = week.t
  , year       = unique(find_clusters$year)[i]
  , samp_clust = clust.t
  )
  
  if (i == 1) {
    samp.clust.a <- samp.clust
  } else {
    samp.clust.a <- rbind(samp.clust.a, samp.clust)
  }

}

sampled_weeks <- left_join(sampled_weeks, samp.clust.a)

n_occ.m <- sampled_weeks %>% 
  dplyr::select(-week) %>%
  group_by(year, samp_clust, Site) %>%
  summarize(n_occ = n()) %>%
  mutate(Site = factor(Site, levels = u_sites)) %>%
  arrange(Site) %>%
  pivot_wider(values_from = n_occ, names_from = Site) %>% 
  arrange(year) %>%
  ungroup() %>%
  dplyr::select(-year, -samp_clust) %>% 
  as.matrix()

capt_history %<>% left_join(., samp.clust.a)

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
  group_by(Mark, year, samp_clust, Site) %>% 
  summarize(capt = sum(captured)) %>% 
  ungroup(year, samp_clust) %>%
  ## And thenin all future times from the current time these individuals _could_ be here
  mutate(capt = cumsum(capt)) %>% 
  mutate(capt = ifelse(capt > 0, 1, 0)) 

## For each individual extract which time periods we do not know if an individual was present or not
for (k in 1:n_sites) {
  p_zeros <- matrix(data = 0, nrow = n_ind.per[k, 1], ncol = sum(n_occ[, u_sites[k]]))
  for (i in 1:n_ind.per[k, 1]) {
    tdat <- first_capt %>% filter(Site == u_sites[k])
    tdat %<>% filter(Mark == unique(tdat$Mark)[i])
    ## For each individual repeat 0s and 1s for each sampling occasions in all years that they were
    ## never captured (0s) and captured (1s) or after a first capture year (1s) 
    rep.t <- n_occ.m[, u_sites[k]]
    rep.t <- rep.t[which(rep.t != 0)]
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
## individual in each primary period). Need an index for these scaling factors
capt_history.p %<>% 
  mutate(gamma_index = paste(interaction(Mark, year, samp_clust),"a",sep="_")) %>% 
  mutate(gamma_index = factor(gamma_index, levels = unique(gamma_index))) %>%
  mutate(gamma_index = as.numeric(gamma_index))

### --- Data for survival (.phi for survival) --- ###

## phi not calculable on the last time step so drop it
last_week        <- capt_history %>%
  group_by(Site) %>% 
  filter(sampled == 1) %>% 
  summarize(last_week = max(week_year))

capt_history.phi <- capt_history %>% 
  left_join(., last_week) %>%
  filter(week_year != last_week, sampled == 1) %>% 
  ungroup()

## Determine the number of time periods that elapse between back to back samples.
## Do this with the data without the last date dropped (as phi is survival to the next)
## and then add to capt_history.phi
time_gaps <- (capt_history %>% 
    filter(sampled == 1) %>%
    group_by(Site, Mark, year) %>% 
    mutate(time_gaps =  samp_clust - lag(samp_clust, 1)) %>% 
    mutate(time_gaps = ifelse(is.na(time_gaps), 0, time_gaps)) %>%
    ungroup(year) %>%
    mutate(time_gaps = c(time_gaps[-1], NA)) %>% 
    filter(!is.na(time_gaps)))$time_gaps

capt_history.phi %<>% mutate(time_gaps = time_gaps)

## Offseason vector 
capt_history.phi %<>% 
  group_by(Site, Mark) %>%
  mutate(offseason = year - lag(year, 1)) %>% 
  mutate(offseason = ifelse(is.na(offseason), 0, offseason)) %>%
  mutate(offseason = c(offseason[-1], 0)) %>% 
  ungroup()

## Which entries of phi correspond to a new individual (the first entry for each individual)
phi_first_index <- (capt_history.phi %>% mutate(index = seq(n())) %>% 
    group_by(Mark) %>% 
    summarize(first_index = min(index)))$first_index

## Indices for which entries of phi must be 0. See above notes for p_zeros. Similar idea here, but 
## here the differentiation for 0s and 1s are for prior to and after an individual was captured for the first time
for (k in 1:n_sites) {
  phi_zeros <- matrix(data = 0, nrow = n_ind.per[k, 1], ncol = sum(n_occ[, u_sites[k]]) - 1)
  
  for (i in 1:n_ind.per[k, 1]) {
    tdat     <- first_capt %>% filter(Site == u_sites[k])
    tdat     %<>% filter(Mark == unique(tdat$Mark)[i])
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

## Index for the summary of individual bd in each time period
capt_history.phi %<>% 
  mutate(X_stat_index = paste(interaction(Mark, year),"a",sep="_")) %>% 
  mutate(X_stat_index = factor(X_stat_index, levels = unique(X_stat_index))) %>%
  mutate(X_stat_index = as.numeric(X_stat_index))

## periods where we assume survival is guaranteed
capt_history.phi %<>% mutate(phi_ones = ifelse(time_gaps == 1 | offseason == 1, 0, 1))

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
  , capt_history   %>% dplyr::select(Mark, week, year, Site, index)
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

## -- Individual specific covariates -- ##

## LOTS of problems with size. For now just 
ind.size <- (capt_history %>% group_by(Mark) %>%
  summarize(size = mean(size, na.rm = T)))$size
ind.size <- scale(ind.size)[, 1]

####
## Run the stan model
####

stan.iter     <- 2000
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin

stan_data     <- list(
  
  ## dimensional indexes 
  n_pop                = n_sites
  , n_ind              = n_ind
  , ind_per_period_bd  = max(capt_history.phi$X_stat_index)
  , ind_per_period_p   = max(capt_history.p$gamma_index) 
  
  , ind_time        = nrow(capt_history)
  , ind_occ         = nrow(capt_history.p)
  , ind_occ_min1    = nrow(capt_history.phi)
  
  ## short vector indexes 
  , ind_occ_size      = rep(colSums(n_occ), n_ind.per)
  , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)
  
  , p_first_index     = p_first_index
  , phi_first_index   = phi_first_index
  
  ## long vector indexes
  , ind_occ_rep       = capt_history.p$Mark
  , periods_occ       = as.numeric(as.factor(capt_history.p$year))
  , p_month           = scale(capt_history.p$week, center = T, scale = F)[, 1] # as.numeric(as.factor(capt_history.p$samp_clust))
  , pop_p             = as.numeric(as.factor(capt_history.p$Site))
  , p_zeros           = capt_history.p$p_zeros
  , p_bd_index        = capt_history.p$p_bd_index
  , gamma_index       = capt_history.p$gamma_index
  
  , ind_occ_min1_rep    = capt_history.phi$Mark
  , offseason           = capt_history.phi$offseason
  , phi_month           = as.numeric(as.factor(capt_history.phi$samp_clust))
  , phi_year            = as.numeric(as.factor(capt_history.phi$year))
  , pop_phi             = as.numeric(as.factor(capt_history.phi$Site))
  , phi_zeros           = capt_history.phi$phi_zeros
  , phi_ones            = capt_history.phi$phi_ones
  , phi_bd_index        = capt_history.phi$phi_bd_index
  , X_stat_index        = capt_history.phi$X_stat_index
  
  , ind_bd_rep          = capt_history$Mark
  , ind_in_pop          = as.numeric(as.factor(capt_history$Site))
  
  ## covariates
  , N_bd            = nrow(capt_history.bd_load)
  , X_bd            = capt_history.bd_load$log_bd_load  
  , X_ind           = capt_history.bd_load$Mark
  # , x_bd_index      = capt_history.bd_load$x_bd_index
  # , bd_first_index  = bd_first_index
  # , bd_last_index   = bd_last_index
  , time_gaps       = capt_history.phi$time_gaps
  
  , ind_size = ind.size
  # , ind_hg   = ind.hg
  
  ## Capture data
  , N_y             = nrow(capt_history.p)
  , y               = capt_history.p$captured
  
  , first           = capture_range$first
  , last            = capture_range$final
  
)

stan.fit  <- stan(
    file    = "CMR_simulation/CMR_collapsed_newts_exp.stan"
  , data    = stan_data
  , chains  = 1
  , cores   = 1
  , refresh = 20
  , iter    = stan.iter            
  , warmup  = stan.burn
  , thin    = stan.thin
  , control = list(adapt_delta = 0.92, max_treedepth = 12)
)

shinystan::launch_shinystan(stan.fit)

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

####
## Model Diagnostics
####

## First just check that phi is estimated when it should be
capt_history.phi %<>% mutate(
  pred_phi = colMeans(stan.fit.samples$phi)
)

capt_history.phi %>% 
  dplyr::select(Mark, year, samp_clust, time_gaps, offseason, phi_zeros, phi_ones, pred_phi) %>% 
  filter(Mark == 9) %>% 
  as.data.frame()

capt_history.p %<>% mutate(
  pred_p   = colMeans(stan.fit.samples$p)
, pred_chi = colMeans(stan.fit.samples$chi)
)

capt_history.p %>% 
  dplyr::select(Mark, year, samp_clust, time_gaps, p_zeros, pred_p, pred_chi) %>% 
  filter(Mark == 9) %>% as.data.frame()

## -- survival over time -- ##

ttg <- seq(1, 10, by = 1)
ttg <- seq(-2, 2, by = .1)
out.pred <- matrix(nrow = 1500, ncol = length(ttg))

for (i in 1:ncol(out.pred)) {
  out.pred[, i] <- plogis(
    stan.fit.samples$beta_phi[, 1] + 
    stan.fit.samples$beta_phi[, 2] * 7 +
    stan.fit.samples$beta_phi[, 3] * ttg[i]
    )
}

reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>%
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(ttg))) %>%
  group_by(gap) %>%
  summarize(
    lwr = quantile(value, 0.2)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.8)
  ) %>% {
  ggplot(., aes(gap, mid)) + 
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
      geom_line(lwd = 2) +
    xlab("Size") +
    ylab("Probability animal remains in population")
  }

## -- survival between years -- ##

ttg <- seq(1, 10, by = 1)
ttg <- seq(-2, 2, by = .1)
out.pred <- matrix(nrow = 1500, ncol = length(ttg))

for (i in 1:ncol(out.pred)) {
   out.pred[, i] <- plogis(
    stan.fit.samples$beta_offseason[, 1] + 
    stan.fit.samples$beta_offseason[, 2] * 7 +
    stan.fit.samples$beta_offseason[, 3] * ttg[i]
    )
}

reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>%
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(ttg))) %>%
  group_by(gap) %>%
  summarize(
    lwr = quantile(value, 0.2)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.8)
  ) %>% {
  ggplot(., aes(gap, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line(lwd = 2) +
    xlab("Size") +
    ylab("Probability animal remains in population")
  }

## -- detection -- ##

ttg <- seq(1, 10, by = 1)
ttg <- seq(-2, 2, by = .1)
out.pred <- matrix(nrow = 1500, ncol = length(ttg))

for (i in 1:ncol(out.pred)) {
   out.pred[, i] <- plogis(
    stan.fit.samples$beta_p[, 1] + 
    stan.fit.samples$beta_p[, 2] * ttg[i] +
    stan.fit.samples$beta_p[, 3] * 0 +
    stan.fit.samples$beta_p_y[, 3]
    )
}

reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>%
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(ttg))) %>%
  group_by(gap) %>%
  summarize(
    lwr = quantile(value, 0.2)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.8)
  ) %>% {
  ggplot(., aes(gap, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line(lwd = 2) +
    xlab("Bd") +
    ylab("Probability animal is detected")
  }

## Strong correlation here, check samples
data.frame(
  stan.fit.samples$beta_p[, 2]
, stan.fit.samples$beta_phi[, 2]
) %>% plot()

## -- individual variation in bd -- ##

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

stan.ind_pred_var <- stan.fit.samples$X %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value) %>%
  group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% arrange(mid) %>%
  mutate(ind = factor(ind, levels = ind))

stan.ind_pred_var %>% {
  ggplot(., aes(ind, mid)) + 
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr)) +
    xlab("Individual") + 
    ylab("Random Effect Deviate") +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = 1, colour = "firebrick3") +
    theme(axis.text.x = element_text(size = 8))
}

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
)) %>% arrange(desc(tot_bd)) %>% 
  mutate(order_real = seq(n()), Mark = as.character(Mark)) %>% 
  filter(!is.na(tot_bd))

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

ind.pred <- reshape2::melt(
  plogis(sweep(stan.fit.samples$bd_ind, 1, stan.fit.samples$beta_offseason[, 2], "*"))
)

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

ind.pred %>% filter(Var2 %in% (
  (ind_order %>% arrange(order_pred))$Mark[1:10]
)) %>% {
  ggplot(., aes(x = value)) + 
    geom_histogram(bins = 50, fill = "firebrick3", alpha = 0.5) + 
    facet_wrap(~Var2)
}


