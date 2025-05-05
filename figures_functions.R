####################################################################
## Functions to extract all needed estimates for manuscript plots ##
####################################################################

source("packages_functions.R")
source("ggplot_theme.R")

extract_estimates <- function(saved_model, sampling, spec_pop_plot_labels, spec_labs, spec_names, fit_mehg, fit_ind_mehg, multi_spec_red = F, which_fit, ...) {
  
stan.fit         <- readRDS(saved_model)
if (which_fit == "pop_mehg") {
  stan.fit.samples <- stan.fit$fitted_model %>% extract()
  capt_history.phi <- stan.fit$capt_history.phi
  capt_history.p   <- stan.fit$capt_history.p
  sampling         <- stan.fit$sampling
} else {
 # stan.fit.samples  <- stan.fit$fitted_model %>% extract()
  stan.fit.samples  <- stan.fit[[4]]
  capt_history.phi  <- stan.fit[[1]]
  capt_history.p    <- stan.fit[[2]]
  sampling          <- stan.fit[[3]]
}

these_specs <- unique(capt_history.phi$Species)
n_specs     <- length(these_specs)
these_pops  <- unique(capt_history.phi$pop_spec)
these_sexes <- c("F", "M")    ## skipping U for now

####
## Plotting Setup: Estimation of generated quantities
####

## Which species are associated with each population
spec_in_pop <- (capt_history.p %>% group_by(pop_spec) %>% slice(1))$Species %>% as.numeric()
## Full data frame of all the combinations of spec, sex, and population seen in the data
 ## (from which to create estimates from model samples)
spec_sex    <- data.frame(
  spec = rep(spec_in_pop, 2)
, sex  = rep(c("F", "M"), each = length(spec_in_pop))
, pop  = rep(seq(length(spec_in_pop)), 2)
)
## Create the model matrix of all of these combinations of species and sex
if (n_specs > 1) {
  if (multi_spec_red) {
spec_sex_mm <- model.matrix(~Sex, capt_history.phi)[, ] %>% as.data.frame() %>% distinct()
spec_sex_mm %<>% mutate(sex = ifelse(SexF == 1 & SexU != 1, "F", (ifelse(SexF != 1 & SexU == 1, "U", "M"))))
  } else {
spec_sex_mm <- model.matrix(~Species + Sex, capt_history.phi)[, ] %>% as.data.frame() %>% distinct()
spec_sex_mm %<>% mutate(sex = ifelse(SexF == 1 & SexU != 1, "F", (ifelse(SexF != 1 & SexU == 1, "U", "M"))))
spec_sex_mm %<>% mutate(spec = rep(seq(n_specs), each = 3))
  }
} else {
spec_sex_mm <- model.matrix(~Sex, capt_history.phi)[, ] %>% as.data.frame() %>% distinct()
spec_sex_mm %<>% mutate(sex = ifelse(SexF == 1 & SexU != 1, "F", (ifelse(SexF != 1 & SexU == 1, "U", "M"))))
}


list_of_CI   <- get_CI(stan.fit.samples, spec_sex, spec_sex_mm, spec_in_pop, fit_mehg, fit_ind_mehg
                      , n_specs, multi_spec_red = F, which_fit)
cleaned_CI   <- clean_CI(list_of_CI, capt_history.p, sampling, fit_mehg, fit_ind_mehg, spec_names
                        , spec_in_pop, which_fit)
pred.vals.gg <- get_envelopes(stan.fit.samples, spec_sex, spec_sex_mm, spec_in_pop, fit_mehg, fit_ind_mehg
                            , n_specs, multi_spec_red = F, which_fit
                            , these_pops, these_sexes, these_specs
                            , spec_pop_plot_labels, spec_names)


return(
  list(
    list_of_CI
  , cleaned_CI
  , pred.vals.gg
  )
)

}
get_CI            <- function(stan.fit.samples, spec_sex, spec_sex_mm, spec_in_pop, fit_mehg, fit_ind_mehg, n_specs, multi_spec_red, which_fit) {
  
int.est <- matrix(
  data = 0
, nrow = dim(stan.fit.samples$z_r)[1]
, ncol = dim(stan.fit.samples$z_r)[3]
)

bd.est <- matrix(
  data = 0
, nrow = dim(stan.fit.samples$z_r)[1]
, ncol = dim(stan.fit.samples$z_r)[3]
)

len.est <- matrix(
  data = 0
  , nrow = dim(stan.fit.samples$z_r)[1]
  , ncol = dim(stan.fit.samples$z_r)[3]
)

if (fit_mehg) {
  
mehg.est <- matrix(
  data = 0
, nrow = dim(stan.fit.samples$z_r)[1]
, ncol = dim(stan.fit.samples$z_r)[3]
)
  
}

if (fit_ind_mehg) {
  
bd.mehg.est <- matrix(
  data = 0
, nrow = dim(stan.fit.samples$z_r)[1]
, ncol = dim(stan.fit.samples$z_r)[3]
)
  
}
 
## Derived quantities here being average survival at the mean of all continuous predictors and
 ## the effect of Bd
for (j in 1:ncol(bd.est)) {
  if (n_specs > 1) {
  if (multi_spec_red) {
   int.est[, j] <- stan.fit.samples$beta_offseason_int[, 1] + stan.fit.samples$z_r[, 1, j]
   bd.est[, j]  <- stan.fit.samples$beta_offseason_bd + stan.fit.samples$z_r[, 2, j]
  } else {
  int.est[, j] <- (sweep(stan.fit.samples$beta_offseason_int[, 1:n_specs], 2
    , spec_sex_mm[spec_sex_mm$spec == spec_in_pop[j] & spec_sex_mm$sex == "M", 1:n_specs] %>% as.matrix() %>% c(), `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 1, j]
   bd.est[, j] <- (sweep(stan.fit.samples$beta_offseason_bd, 2
    , spec_sex_mm[spec_sex_mm$spec == spec_in_pop[j] & spec_sex_mm$sex == "M", 1:n_specs] %>% as.matrix() %>% c(), `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 2, j]
   len.est[, j] <- (sweep(stan.fit.samples$beta_offseason_len, 2
    , spec_sex_mm[spec_sex_mm$spec == spec_in_pop[j] & spec_sex_mm$sex == "M", 1:n_specs] %>% as.matrix() %>% c(), `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 3, j]
  }
  if (fit_mehg) {
    if (which_fit == "ind_mehg") {
    mehg.est[, j] <- (sweep(stan.fit.samples$beta_offseason_mehg[, 1:n_specs], 2
    , spec_sex_mm[spec_sex_mm$spec == spec_in_pop[j] & spec_sex_mm$sex == "M", 1:n_specs] %>% as.matrix() %>% c(), `*`) %>% rowSums()) +
      stan.fit.samples$z_r[, 4, j]
    } else {
    mehg.est[, j] <- (sweep(stan.fit.samples$beta_offseason_mehg[, 1:n_specs], 2
    , spec_sex_mm[spec_sex_mm$spec == spec_in_pop[j] & spec_sex_mm$sex == "M", 1:n_specs] %>% as.matrix() %>% c(), `*`) %>% rowSums())
    }
  }
  if (fit_ind_mehg) {
    bd.mehg.est[, j] <- (sweep(stan.fit.samples$beta_offseason_mehg_bd[, 1:n_specs], 2
    , spec_sex_mm[spec_sex_mm$spec == spec_in_pop[j] & spec_sex_mm$sex == "M", 1:n_specs] %>% as.matrix() %>% c(), `*`) %>% rowSums()) +
      stan.fit.samples$z_r[, 5, j]
  }
  } else {
   int.est[, j] <- stan.fit.samples$beta_offseason_int[, 1] + stan.fit.samples$z_r[, 1, j]
   bd.est[, j]  <- stan.fit.samples$beta_offseason_bd + stan.fit.samples$z_r[, 2, j] 
  if (fit_mehg) {
   mehg.est[, j] <- stan.fit.samples$beta_offseason_mehg
  }
  if (fit_ind_mehg) {
   bd.mehg.est[, j] <- stan.fit.samples$beta_offseason_mehg_bd + stan.fit.samples$z_r[, 5, j]
  }
    
  }
}

return(
list(
  int.est
, bd.est
, {if(fit_mehg){mehg.est}else{NULL}}
, {if(fit_ind_mehg){bd.mehg.est}else{NULL}}
, len.est
)
)
 
}
get_envelopes     <- function(stan.fit.samples, spec_sex, spec_sex_mm, spec_in_pop
                              , fit_mehg, fit_ind_mehg, n_specs, multi_spec_red, which_fit
                              , these_pops, these_sexes, these_specs
                              , spec_pop_plot_labels, spec_names) {
 
## loop over each species, sex, pop to estimate the survival of each of the individual types found in the data
for (k in 1:nrow(spec_sex)) {
  
## grid over which to predict
pred.vals <- expand.grid(
  bd   = scale(seq(0, 14, by = 1))[, 1]
, len  = seq(-2, 2, by = 0.4)
, mehg = seq(-2, 2, by = 0.4)
)

## matrix to house the full posterior 
pred.est <- matrix(data = 0, nrow = nrow(pred.vals), ncol = dim(stan.fit.samples$z_r)[1])

## Obtain the model matrix entry for this spec, sex, pop combination
if (n_specs > 1) {
  if (multi_spec_red) {
spec_sex_mm.t <- spec_sex_mm %>% filter(sex == spec_sex$sex[k]) %>% dplyr::select(-sex) %>% as.matrix()
  } else {
spec_sex_mm.t <- spec_sex_mm %>% filter(spec == spec_sex$spec[k], sex == spec_sex$sex[k]) %>% 
  dplyr::select(-spec, -sex) %>% as.matrix()
  }
} else {
spec_sex_mm.t <- spec_sex_mm %>% filter(sex == spec_sex$sex[k]) %>% dplyr::select(-sex) %>% as.matrix()
}

## loop over all posterior samples to calculate the outcome of interest (here survival)
for (j in 1:nrow(pred.est)) {
  
if (n_specs > 1) {
  if (multi_spec_red) {
 pred.est[j, ] <- plogis(
    (sweep(stan.fit.samples$beta_offseason_int, 2, spec_sex_mm.t, `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 1, spec_sex$pop[k]] + 
    (stan.fit.samples$beta_offseason_bd  + stan.fit.samples$z_r[, 2, spec_sex$pop[k]]) * pred.vals$bd[j] +
    (stan.fit.samples$beta_offseason_len + stan.fit.samples$z_r[, 3, spec_sex$pop[k]]) * pred.vals$len[j] +
    stan.fit.samples$beta_offseason_mehg * pred.vals$mehg[j]
 )
  } else {
    if (!fit_ind_mehg) {
 pred.est[j, ] <- plogis(
    (sweep(stan.fit.samples$beta_offseason_int, 2, spec_sex_mm.t, `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 1, spec_sex$pop[k]] + 
    (
      (sweep(stan.fit.samples$beta_offseason_bd, 2, spec_sex_mm.t[1:n_specs], `*`) %>% rowSums()) + 
        stan.fit.samples$z_r[, 2, spec_sex$pop[k]]
      ) * pred.vals$bd[j] +
    (
      (sweep(stan.fit.samples$beta_offseason_len, 2, spec_sex_mm.t[1:n_specs], `*`) %>% rowSums()) +
        stan.fit.samples$z_r[, 3, spec_sex$pop[k]]
      ) * pred.vals$len[j] +
    (
      sweep(stan.fit.samples$beta_offseason_mehg, 2, spec_sex_mm.t[1:n_specs], `*`) %>% rowSums()
      ) * pred.vals$mehg[j]
 )
    } else {
 pred.est[j, ] <- plogis(
    (sweep(stan.fit.samples$beta_offseason_int, 2, spec_sex_mm.t, `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 1, spec_sex$pop[k]] + 
    (
      (sweep(stan.fit.samples$beta_offseason_bd, 2, spec_sex_mm.t[1:n_specs], `*`) %>% rowSums()) + 
        stan.fit.samples$z_r[, 2, spec_sex$pop[k]]
      ) * pred.vals$bd[j] +
    (
      (sweep(stan.fit.samples$beta_offseason_len, 2, spec_sex_mm.t[1:n_specs], `*`) %>% rowSums()) +
        stan.fit.samples$z_r[, 3, spec_sex$pop[k]]
      ) * pred.vals$len[j] +
    (
      sweep(stan.fit.samples$beta_offseason_mehg, 2, spec_sex_mm.t[1:n_specs], `*`) %>% rowSums() +
        stan.fit.samples$z_r[, 4, spec_sex$pop[k]]
      ) * pred.vals$mehg[j] +
    (
      sweep(stan.fit.samples$beta_offseason_mehg_bd, 2, spec_sex_mm.t[1:n_specs], `*`) %>% rowSums() +
        stan.fit.samples$z_r[, 5, spec_sex$pop[k]]
      ) * pred.vals$mehg[j] * pred.vals$bd[j]
 )
    }
  }

} else {
    if (!fit_ind_mehg) {
 pred.est[j, ] <- plogis(
    (sweep(stan.fit.samples$beta_offseason_int, 2, spec_sex_mm.t, `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 1, spec_sex$pop[k]] + 
    (stan.fit.samples$beta_offseason_bd + stan.fit.samples$z_r[, 2, spec_sex$pop[k]]) * pred.vals$bd[j] +
    (stan.fit.samples$beta_offseason_len + stan.fit.samples$z_r[, 3, spec_sex$pop[k]]) * pred.vals$len[j]
 )
  } else {
 pred.est[j, ] <- plogis(
    (sweep(stan.fit.samples$beta_offseason_int, 2, spec_sex_mm.t, `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 1, spec_sex$pop[k]] + 
    (stan.fit.samples$beta_offseason_bd + stan.fit.samples$z_r[, 2, spec_sex$pop[k]]) * pred.vals$bd[j] +
    (stan.fit.samples$beta_offseason_len + stan.fit.samples$z_r[, 3, spec_sex$pop[k]]) * pred.vals$len[j] +
    (stan.fit.samples$beta_offseason_mehg + stan.fit.samples$z_r[, 4, spec_sex$pop[k]]) * pred.vals$mehg[j] +
    (stan.fit.samples$beta_offseason_mehg_bd + stan.fit.samples$z_r[, 5, spec_sex$pop[k]]) * pred.vals$mehg[j] * pred.vals$bd[j]
 )
  }
}
  
}

## Add in the details about the spec, sex, and pop 
if (n_specs > 1) {
pred.vals %<>% mutate(
  pop  = spec_sex$pop[k]
, spec = spec_sex$spec[k]
, sex  = spec_sex$sex[k]
)
} else {
pred.vals %<>% mutate(
  pop  = spec_sex$pop[k]
, sex  = spec_sex$sex[k]
)
}

pred.vals <- cbind(pred.vals, pred.est) 

if (n_specs > 1) {
pred.vals %<>% pivot_longer(., c(-bd, -len, -mehg, -pop, -spec, -sex), names_to = "iter", values_to = "est")
pred.vals.gg <- pred.vals %>% group_by(bd, len, mehg, spec, pop, sex)
} else {
pred.vals %<>% pivot_longer(., c(-bd, -len, -mehg, -pop, -sex), names_to = "iter", values_to = "est")  
pred.vals.gg <- pred.vals %>% group_by(bd, len, mehg, pop, sex)
}

pred.vals.gg %<>% summarize(
    lwr   = quantile(est, 0.025)
  , lwr_n = quantile(est, 0.200)
  , mid   = quantile(est, 0.500)
  , upr_n = quantile(est, 0.800)
  , upr   = quantile(est, 0.975)
  )

print(paste("Through", k, "population:sex", sep = " "))

## stick it all together
if (k == 1) {
pred.vals.gg.f <- pred.vals.gg
} else {
pred.vals.gg.f <- rbind(pred.vals.gg.f, pred.vals.gg)
}
}

pred.vals.gg.f %<>% ungroup() %>% mutate(
  pop  = plyr::mapvalues(pop, from = unique(pred.vals.gg.f$pop), to = spec_pop_plot_labels)
, sex  = plyr::mapvalues(sex, from = unique(pred.vals.gg.f$sex), to = these_sexes)
)

if (n_specs > 1) {
pred.vals.gg.f %<>% mutate(spec = plyr::mapvalues(spec, from = unique(pred.vals.gg.f$spec), to = spec_names))
pred.vals.gg.f %<>% mutate(pop =  as.factor(pop), spec = as.factor(spec))
} else {
pred.vals.gg.f %<>% mutate(pop =  as.factor(pop))
}

pred.vals.gg <- pred.vals.gg.f; rm(pred.vals.gg.f)

return(pred.vals.gg)
   
}
clean_CI          <- function(list_of_CI, capt_history.p, sampling, fit_mehg, fit_ind_mehg, spec_names, spec_in_pop, which_fit) {
  
int.est <- list_of_CI[[1]]
int.est <- reshape2::melt(int.est)
names(int.est) <- c("Sample", "Population", "Value")
int.est %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = spec_names)
)

bd.est <- list_of_CI[[2]]
bd.est <- reshape2::melt(bd.est)
names(bd.est) <- c("Sample", "Population", "Value")
bd.est %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = spec_names)
)

if (fit_mehg) {
  mehg.est <- list_of_CI[[3]]
  mehg.est <- reshape2::melt(mehg.est)
  names(mehg.est) <- c("Sample", "Population", "Value")
mehg.est %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = spec_names)
)
}

if (fit_ind_mehg) {
  bd.mehg.est <- list_of_CI[[4]]
  bd.mehg.est <- reshape2::melt(bd.mehg.est)
  names(bd.mehg.est) <- c("Sample", "Population", "Value")
bd.mehg.est %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = spec_names)
)
}

len.est <- list_of_CI[[5]]
len.est <- reshape2::melt(len.est)
names(len.est) <- c("Sample", "Population", "Value")
len.est %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = spec_names)
)
  
## A data frame for later to help scale the effect of Bd to the more interpretable probability scale
int.est2 <- int.est %>% group_by(Population, Species) %>% summarize(mid   = quantile(Value, 0.500))
## Second option of doing so
bd.est2  <- left_join(
  int.est %>% rename(int = Value)
, bd.est  %>% rename(bd = Value)
)

int.est %<>% mutate(
  Value = plogis(Value)
  ) %>% 
  group_by(Population, Species) %>% summarize(
    lwr   = quantile(Value, 0.025)
  , lwr_n = quantile(Value, 0.200)
  , mid   = quantile(Value, 0.500)
  , upr_n = quantile(Value, 0.800)
  , upr   = quantile(Value, 0.975)
) %>% ungroup() %>% mutate(
  pop_spec = unique(capt_history.p$pop_spec)
) %>% mutate(
  CI_width = upr - lwr
)

bd.est %<>% group_by(Population, Species) %>% summarize(
    lwr   = quantile(Value, 0.025)
  , lwr_n = quantile(Value, 0.200)
  , mid   = quantile(Value, 0.500)
  , upr_n = quantile(Value, 0.800)
  , upr   = quantile(Value, 0.975)
) %>% ungroup() %>% mutate(
  pop_spec = unique(capt_history.p$pop_spec)
) %>% mutate(
  CI_width = upr - lwr
)

len.est %<>% group_by(Population, Species) %>% summarize(
  lwr   = quantile(Value, 0.025)
  , lwr_n = quantile(Value, 0.200)
  , mid   = quantile(Value, 0.500)
  , upr_n = quantile(Value, 0.800)
  , upr   = quantile(Value, 0.975)
) %>% ungroup() %>% mutate(
  pop_spec = unique(capt_history.p$pop_spec)
) %>% mutate(
  CI_width = upr - lwr
)

if (fit_mehg) {
 
mehg.est %<>% group_by(Population, Species) %>% summarize(
    lwr   = quantile(Value, 0.025)
  , lwr_n = quantile(Value, 0.200)
  , mid   = quantile(Value, 0.500)
  , upr_n = quantile(Value, 0.800)
  , upr   = quantile(Value, 0.975)
) %>% ungroup() %>% mutate(
  pop_spec = unique(capt_history.p$pop_spec)
) %>% mutate(
  CI_width = upr - lwr
)
   
}

if (fit_ind_mehg) {
  
bd.mehg.est %<>% group_by(Population, Species) %>% summarize(
    lwr   = quantile(Value, 0.025)
  , lwr_n = quantile(Value, 0.200)
  , mid   = quantile(Value, 0.500)
  , upr_n = quantile(Value, 0.800)
  , upr   = quantile(Value, 0.975)
) %>% ungroup() %>% mutate(
  pop_spec = unique(capt_history.p$pop_spec)
) %>% mutate(
  CI_width = upr - lwr
)
  
}

bd.est2 %<>% mutate(
  int.p = plogis(int)
, bd.p  = plogis(int + bd)) %>% 
  mutate(
  rel_surv = bd.p / int.p
  )

bd.est2 %<>% group_by(Population, Species) %>% summarize(
    lwr   = quantile(rel_surv, 0.025)
  , lwr_n = quantile(rel_surv, 0.200)
  , mid   = quantile(rel_surv, 0.500)
  , upr_n = quantile(rel_surv, 0.800)
  , upr   = quantile(rel_surv, 0.975)
) %>% ungroup() %>% mutate(
  pop_spec = unique(capt_history.p$pop_spec)
) %>% mutate(
  CI_width = upr - lwr
)

## quick aside to calculate some metrics of effort vs CI
## Calculate a number of metrics of effort for some coarse correlations between effort and CI width
num_samples  <- sampling %>% group_by(pop_spec) %>% summarize(n_dates = n_distinct(CaptureDate)) %>%
  arrange(desc(n_dates)) %>% mutate(pop_spec = factor(pop_spec, levels = unique(pop_spec)))

recaps <- capt_history.p %>% group_by(pop_spec, Mark) %>% 
  filter(captured == 1) %>%
  summarize(recaps = sum(captured)) %>% 
  mutate(recaptured = ifelse(recaps > 1, 1, 0)) %>% 
  ungroup(Mark) %>%
  summarize(
    inds_capt  = n_distinct(Mark)
  , recapt_ind = sum(recaptured)
  , caps_per_ind = mean(recaps)) %>% 
  ungroup() 

n_swabs <- capt_history.p %>% group_by(pop_spec, Mark) %>% 
  filter(swabbed == 1) %>%
  summarize(swabbs = sum(swabbed)) %>% 
  mutate(reswabbed = ifelse(swabbs > 1, 1, 0)) %>% 
  ungroup(Mark) %>%
  summarize(
    inds_swabbed  = n_distinct(Mark)
  , reswabbed_ind = sum(reswabbed)
  , swabs_per_ind = mean(swabbs)) %>% 
  ungroup() 

int.est %<>% left_join(., num_samples) %>% left_join(., recaps)
bd.est  %<>% left_join(., num_samples) %>% left_join(., recaps) %>% 
  left_join(., n_swabs) %>% mutate(Population = spec_pop_plot_labels)

int.est.gg <- int.est %>% mutate(Population = plyr::mapvalues(Population
  , from = unique(Population)
  , to = spec_pop_plot_labels)
) %>% arrange(desc(mid)) %>% 
  mutate(pop_spec = factor(pop_spec, levels = pop_spec)) 

len.est.gg <- len.est %>%
  left_join(., num_samples) %>% left_join(., recaps) %>% 
  mutate(Population = plyr::mapvalues(Population
  , from = unique(Population)
  , to = spec_pop_plot_labels)
) %>% arrange(desc(mid)) %>% 
  mutate(pop_spec = factor(pop_spec, levels = pop_spec)) 

bd.est.gg <- bd.est %>% arrange(desc(mid)) %>% 
  mutate(pop_spec = factor(pop_spec, levels = pop_spec)) 

if (fit_mehg) {
  
mehg.est.gg <- mehg.est %>% mutate(Population = plyr::mapvalues(Population
  , from = unique(Population)
  , to = spec_pop_plot_labels)
) %>% arrange(desc(mid)) %>% 
  mutate(pop_spec = factor(pop_spec, levels = pop_spec))   
  
}

if (fit_ind_mehg) {
  
bd.mehg.est.gg <- bd.mehg.est %>% mutate(Population = plyr::mapvalues(Population
  , from = unique(Population)
  , to = spec_pop_plot_labels)
) %>% arrange(desc(mid)) %>% 
  mutate(pop_spec = factor(pop_spec, levels = pop_spec))   
  
}

bd.est2 <- int.est2 %>% ungroup() %>% 
  mutate(Population = plyr::mapvalues(Population
  , from = unique(Population)
  , to = spec_pop_plot_labels)
) %>% rename(mid_p = mid) %>% left_join(.
  , 
  bd.est # %>% mutate(Population = plyr::mapvalues(Population
#  , from = unique(Population)
#  , to = spec_pop_plot_labels
#)
# )
  ) %>% mutate(
    lwr   = plogis(mid_p) / plogis(mid_p - lwr)    
  , lwr_n = plogis(mid_p) / plogis(mid_p - lwr_n)
  , mid   = plogis(mid_p) / plogis(mid_p - mid)    
  , upr_n = plogis(mid_p) / plogis(mid_p - upr_n)  
  , upr   = plogis(mid_p) / plogis(mid_p - upr)    
  ) %>% arrange(desc(mid)) %>% 
  mutate(pop_spec = factor(pop_spec, levels = pop_spec)) 

bd.est2.gg <- bd.est2 %>% 
#  mutate(Population = plyr::mapvalues(Population
#  , from = unique(Population)
#  , to = spec_pop_plot_labels
#)) %>% 
  arrange(desc(mid)) %>% 
  mutate(pop_spec = factor(pop_spec, levels = pop_spec)) 

return(
list(
  int.est        = int.est %>% mutate(which_fit = which_fit)
, int.est2       = int.est2 %>% mutate(which_fit = which_fit)
, int.est.gg     = int.est.gg %>% mutate(which_fit = which_fit)
, bd.est         = bd.est %>% mutate(which_fit = which_fit)
, bd.est2        = bd.est2 %>% mutate(which_fit = which_fit)
, bd.est2.gg     = bd.est2.gg %>% mutate(which_fit = which_fit)
, bd.est.gg      = bd.est.gg %>% mutate(which_fit = which_fit)
, mehg.est       = {if(fit_mehg){mehg.est %>% mutate(which_fit = which_fit)}else{NULL}}
, mehg.est.gg    = {if(fit_mehg){mehg.est.gg %>% mutate(which_fit = which_fit)}else{NULL}}
, bd.mehg.est    = {if(fit_ind_mehg){bd.mehg.est %>% mutate(which_fit = which_fit)}else{NULL}}
, bd.mehg.est.gg = {if(fit_ind_mehg){bd.mehg.est.gg %>% mutate(which_fit = which_fit)}else{NULL}}
, len.est        = len.est %>% mutate(which_fit = which_fit)
, len.est.gg     = len.est.gg %>% mutate(which_fit = which_fit)
, recaps         = recaps
, nswabs         = n_swabs
)
)

}

####--------
## Fit with no ind-mehg
####

which_fit   <- "pop_mehg"
saved_model <- "fits/pop_mehg_full.Rds"

spec_pop_plot_labels <- c(
  "FL - SMNWR East"
, "FL - SMNWR West"
, "WY - Blackrock"
, "MT - Two Medicine"
, "WI - Kettle Moraine"
, "WI - Mud Lake"
, "PA - Scotia Barrens"
, "FL - SMNWR West"
, "MA - Springfield"
, "CA - Summit Meadow"
    )

spec_labs <- c(
  expression(italic("Ambystoma cingulatum")~"[Intercept]")
, expression(italic("Anaxyrus boreas"))
, expression(italic("Notophthalmus viridescens"))
, expression(italic("Rana spp."))
)

spec_names <- c(
  "Ambystoma cingulatum"
, "Anaxyrus boreas"
, "Notophthalmus viridescens"
, "Rana spp."
)

pop_mehg_preds <- extract_estimates(saved_model, NULL, spec_pop_plot_labels, spec_labs, spec_names, fit_mehg = T, fit_ind_mehg = F
                                  , multi_spec_red = F, which_fit)


####--------
## Fit with ind-mehg
####

which_fit   <- "ind_mehg"
saved_model <- "fits/ind_mehg_full.Rds"
#saved_model <- "fits/full_stan_fit.Rds"

spec_pop_plot_labels <- c(
  "MT - Jones Pond"
, "CA - Sonoma Mountain"
, "CO - Lily Pond"
, "CO - Matthews Pond"
, "OR - Dilman Meadows"
, "CA - Fox Creek"
, "MT - Jones Pond"
, "OR - Little Three Creeks"
, "MT - Lost Horse"
, "CA - San Francisquito"
    )

spec_labs <- c(
  expression(italic("Anaxyrus boreas"))
, expression(italic("Pseudacris maculata"))
, expression(italic("Rana")~"spp.")
)

spec_names <- c(
  "Anaxyrus boreas"
, "Pseudacris maculata"
, "Rana spp."
)

ind_mehg_preds <- extract_estimates(saved_model, sampling = sampling, spec_pop_plot_labels, spec_labs, spec_names, fit_mehg = T
                                  , multi_spec_red = F, fit_ind_mehg = T, which_fit)

####--------
## Plots with these fits
####

#### Figure 1 ------
 ## Figure 1: Coefficient estimates (on logit scale) for the effect of Bd on survival
  ## -- bd.est.gg

bd.est.gg <- rbind(
  ind_mehg_preds[[2]]$bd.est.gg
, pop_mehg_preds[[2]]$bd.est.gg
) %>% arrange(desc(mid)) %>% 
  mutate(
    pop_spec = factor(pop_spec, levels = unique(pop_spec))
  , Species  = factor(Species, levels = unique(Species) %>% as.character() %>% sort())
  )

spec_labs <- c(
  expression(italic("Ambystoma cingulatum"))
, expression(italic("Anaxyrus boreas"))
, expression(italic("Notophthalmus viridescens"))
, expression(italic("Pseudacris maculata"))
, expression(italic("Rana spp."))
)

bd.est.gg %>% {
  ggplot(., aes(mid, pop_spec, colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_point(size = 2.5) +
    scale_color_manual(
      values = RColorBrewer::brewer.pal(8,"Dark2")[c(6,7,3,4,1)]
    , labels  = spec_labs
    ) +
    scale_y_discrete(labels = bd.est.gg$Population) +
    xlab("Effect of Bd on Survival (Logit Scale)") +
    scale_x_continuous(breaks = c(-2.5, -1.5, -0.75, 0, 0.75, 1.5, 2.5)
                       , lim = c(-4.4, 4.9)) +
    ylab("Population") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , axis.text.x = element_text(size = 10)
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13)
      , legend.text.align = 0)
}

#### Figure 2 ------
 ## Figure 2: Survival (y-axis) as a function of Bd load (x-axis) across populations (faceted)
  ## -- ind_mehg_preds[[3]] (bd)

pred.vals.gg <- rbind(
  ind_mehg_preds[[3]] %>% mutate(
    bd_raw = (bd * 4.187317 + 5.332455) %>% exp()
  )
, pop_mehg_preds[[3]] %>% mutate(
    bd_raw = (bd * 4.289586 + 6.24991) %>% exp()
  )
) 

pred.vals.gg %<>% mutate(
  spec = factor(spec, levels = c(
    unique(pred.vals.gg$spec) %>% as.character() %>% sort()
  ))
)

pred.vals.gg %>% 
  filter(sex == "M", len == 0, mehg == 0) %>% {
  ggplot(., aes(bd_raw, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr
      , fill = spec, colour = spec
      ), alpha = 0.3) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n
      , fill = spec, colour = spec
      ), alpha = 0.3) +
    geom_line(aes(
      colour = spec
      ), size = 1) + 
    scale_colour_manual(
      name   = "Species"
    , values = brewer.pal(8,"Dark2")[c(6,7,3,4,1)]
    , labels = spec_labs
      ) +
    scale_fill_manual(
      name   = "Species"
    , values = brewer.pal(8,"Dark2")[c(6,7,3,4,1)]
    , labels = spec_labs) +
      scale_x_log10() +
    #scale_x_continuous(breaks = c(-1.40, -0.75, 0, 0.75, 1.40)) +
    facet_wrap(~pop) +
    theme(
    strip.text.x    = element_text(size = 11)
  , axis.text.y     = element_text(size = 12)
  , axis.text.x     = element_text(size = 11)
  , legend.key.size = unit(0.65, "cm")
  , legend.text     = element_text(size = 11)
  , legend.title    = element_text(size = 13)
  , legend.text.align = 0
    ) +
   # xlab("Bd Load (scaled)") +
    xlab("Bd Load (log10 Copy Number)") +
    ylab("Between-Season Survival")
}

#### Figure 3 ------
 ## Figure 3: Coefficient estimates (on logit scale) for the effect of MeHg on survival
  ## -- mehg.est
  ## -- bd.mehg.est

mehg.est.gg <- rbind(
  ind_mehg_preds[[2]]$mehg.est.gg
, pop_mehg_preds[[2]]$mehg.est.gg
) %>% mutate(
  coef = "MeHg Main Effect"
) %>% rbind(
  .
, ind_mehg_preds[[2]]$bd.mehg.est.gg %>% mutate(coef = "MeHg-Bd Interaction")
)

mehg.est.gg %<>% mutate(
  Species = factor(Species, levels = c(
    unique(mehg.est.gg$Species) %>% as.character() %>% sort()
  ))
, Population = factor(Population, levels = c(
  unique(bd.est.gg$Population)
))
)

mehg.est.gg %>% {
  ggplot(., aes(mid, pop_spec, colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_point(size = 2.5) +
    scale_color_manual(values = brewer.pal(8,"Dark2")[c(6,7,3,4,1)]) +
    scale_y_discrete(labels = mehg.est.gg$Population) +
    xlab("Effect of MeHg on Survival (Logit Scale)") +
    scale_x_continuous(breaks = c(-2.5, -1.5, -0.75, 0, 0.75, 1.5, 2.5)) +
    ylab("Population") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
    facet_wrap(~coef) +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , axis.text.x = element_text(size = 9)
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13)
      , legend.text.align = 0)
}


#### Figure S8 ------
 ## Figure S8: Proportional increase/decrease in survival given 1 SD change in Bd 
  ## -- bd.est2.gg

bd.est2.gg <- rbind(
  ind_mehg_preds[[2]]$bd.est2.gg
, pop_mehg_preds[[2]]$bd.est2.gg
) %>% arrange(desc(mid)) %>% 
  mutate(
    pop_spec = factor(pop_spec, levels = unique(pop_spec))
  , Species  = factor(Species, levels = unique(Species) %>% as.character() %>% sort())
  )

spec_labs <- c(
  expression(italic("Ambystoma cingulatum"))
, expression(italic("Anaxyrus boreas"))
, expression(italic("Notophthalmus viridescens"))
, expression(italic("Pseudacris maculata"))
, expression(italic("Rana spp."))
)

bd.est2.gg %>% {
  ggplot(., aes(mid, pop_spec, colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_point(size = 2.5) +
    scale_color_manual(
      values = brewer.pal(8,"Dark2")[c(6,7,3,4,1)]
    , labels  = spec_labs
    ) +
    scale_y_discrete(labels = bd.est2.gg$Population) +
    xlab("Effect of Bd on Survival (Logit Scale)") +
   # scale_x_continuous(breaks = c(-2.5, -1.5, -0.5, 0, 0.5, 1.5, 2.5)) +
    scale_x_log10(breaks = c(0.50, 0.66, 1, 2, 5, 10, 20, 30)) +
    ylab("Population") +
    geom_vline(xintercept = 1, linetype = "dashed", size = 0.4) +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , axis.text.x = element_text(size = 10)
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13)
      , legend.text.align = 0)
}

#### Figure S10 ------
 ## Figure S10: Average survival on probability scale
  ## -- int.est2 (I think)

int.est.gg <- rbind(
  ind_mehg_preds[[2]]$int.est.gg
, pop_mehg_preds[[2]]$int.est.gg
) %>% arrange(desc(mid)) %>% 
  mutate(
    pop_spec = factor(pop_spec, levels = unique(pop_spec))
  , Species  = factor(Species, levels = unique(Species) %>% as.character() %>% sort())
  )

spec_labs <- c(
  expression(italic("Ambystoma cingulatum"))
, expression(italic("Anaxyrus boreas"))
, expression(italic("Notophthalmus viridescens"))
, expression(italic("Pseudacris maculata"))
, expression(italic("Rana spp."))
)

int.est.gg %>% {
  ggplot(., aes(mid, pop_spec, colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_point(size = 2.5) +
    scale_color_manual(
      values = brewer.pal(8,"Dark2")[c(6,7,3,4,1)]
    , labels  = spec_labs
    ) +
    scale_y_discrete(labels = int.est.gg$Population) +
    xlab("Average Survival (Probability Scale)") +
    ylab("Population") +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , axis.text.x = element_text(size = 10)
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13)
      , legend.text.align = 0)
}


#### Figure S11 ------
 ## Figure S11: Length effect (estimates on logit scale)
  ## -- ind_mehg_preds[[3]] (len)

pred.vals.gg <- rbind(
  ind_mehg_preds[[3]]
, pop_mehg_preds[[3]]
) 

pred.vals.gg %<>% mutate(
  spec = factor(spec, levels = c(
    unique(pred.vals.gg$spec) %>% as.character() %>% sort()
  ))
)

pred.vals.gg %>% filter(sex == "M", bd == 0, mehg == 0) %>% {
  ggplot(., aes(len, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr
      , fill = spec, colour = spec
      ), alpha = 0.3) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n
      , fill = spec, colour = spec
      ), alpha = 0.3) +
    geom_line(aes(
      colour = spec
      ), size = 1) + 
    scale_colour_manual(
      name   = "Species"
    , values = brewer.pal(8,"Dark2")[c(6,7,3,4,1)]
    , labels = spec_labs
      ) +
    scale_fill_manual(
      name   = "Species"
    , values = brewer.pal(8,"Dark2")[c(6,7,3,4,1)]
    , labels = spec_labs) +
    scale_x_continuous(breaks = c(-1.40, -0.75, 0, 0.75, 1.40)) +
    facet_wrap(~pop) +
    theme(
    strip.text.x    = element_text(size = 11)
  , axis.text.y     = element_text(size = 12)
  , axis.text.x     = element_text(size = 9)
  , legend.key.size = unit(0.65, "cm")
  , legend.text     = element_text(size = 11)
  , legend.title    = element_text(size = 13)
  , legend.text.align = 0
    ) +
    xlab("Length (scaled)") +
    ylab("Between-Season Survival")
}

#### Figure SX -----
## Duplicated from above but for length and not bd

len.est.gg <- rbind(
  ind_mehg_preds[[2]]$len.est.gg
  , pop_mehg_preds[[2]]$len.est.gg
) %>% arrange(desc(mid)) %>% 
  mutate(
    pop_spec = factor(pop_spec, levels = unique(pop_spec))
    , Species  = factor(Species, levels = unique(Species) %>% as.character() %>% sort())
  )

spec_labs <- c(
  expression(italic("Ambystoma cingulatum"))
  , expression(italic("Anaxyrus boreas"))
  , expression(italic("Notophthalmus viridescens"))
  , expression(italic("Pseudacris maculata"))
  , expression(italic("Rana spp."))
)

len.est.gg %>% {
  ggplot(., aes(mid, pop_spec, colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_point(size = 2.5) +
    scale_color_manual(
      values = RColorBrewer::brewer.pal(8,"Dark2")[c(6,7,3,4,1)]
      , labels  = spec_labs
    ) +
    scale_y_discrete(labels = bd.est.gg$Population) +
    xlab("Effect of SVL on Survival (Logit Scale)") +
    #    scale_x_continuous(breaks = c(-2.5, -1.5, -0.75, 0, 0.75, 1.5, 2.5)
    #                       , lim = c(-4.4, 4.9)) +
    ylab("Population") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
    theme(axis.text.y = element_text(size = 11)
          , legend.key.size = unit(0.6, "cm")
          , axis.text.x = element_text(size = 10)
          , legend.text = element_text(size = 11)
          , legend.title = element_text(size = 13)
          , legend.text.align = 0)
}

#### Figure S12 ------
 ## Figure S12: Within-year survival
  ## -- STILL TO DO


#### Figure SX-SY ------
  ## Population sizes

pop_size_ests <- capt_history.p %>% 
  ungroup() %>% 
  group_by(pop_spec, capture_date) %>%
  slice(1) %>% 
  dplyr::select(pop_spec, capture_date) %>% 
  ungroup() %>%
  mutate(
      lwr   = apply(pop_size, 2, FUN = function(x) quantile(x, 0.025))
    , lwr_n = apply(pop_size, 2, FUN = function(x) quantile(x, 0.200))
    , mid   = apply(pop_size, 2, FUN = function(x) quantile(x, 0.500))
    , upr_n = apply(pop_size, 2, FUN = function(x) quantile(x, 0.800))
    , upr   = apply(pop_size, 2, FUN = function(x) quantile(x, 0.975))
  ) %>% ungroup() %>% 
  group_by(pop_spec) %>% 
  mutate(ss = seq(n())) %>% 
  filter(ss != min(ss)) %>% ungroup()

#pop_size_ests[pop_size_ests$capt_per_day == 0, 5:9] <- NA

#x_labs <- pop_size_ests %>% 
#  group_by(pop_spec) %>% 
#  filter(ss %in% seq(min(ss), max(ss), by = 5))

pop_size_ests$spec <- apply(pop_size_ests$pop_spec %>% matrix, 1, FUN = function(x) strsplit(x, "[.]")[[1]][1])

pop_size_ests %<>% mutate(pop = pop_spec) 
pop_size_ests %<>% mutate(
  pop  = plyr::mapvalues(pop, from = unique(pop_size_ests$pop)
                         , to = spec_pop_plot_labels
  )
  , spec = plyr::mapvalues(spec, from = unique(pop_size_ests$spec)
                           , to = spec_names)
)

pop_size_ests %<>% unite(pop, pop, spec, sep = " - ")

pop_size_ests1 <- pop_size_ests
pop_size_ests2 <- pop_size_ests

pop_size_ests.f <- rbind(pop_size_ests1, pop_size_ests2)

pop_size_ests.f$spec <- apply(pop_size_ests.f$pop %>% matrix, 1, FUN = function(x) strsplit(x, " - ")[[1]][3])
pop_size_ests.f$ppop <- apply(pop_size_ests.f$pop %>% matrix, 1, FUN = function(x) {
  tt <- strsplit(x, " - ")[[1]][1:2]
  paste(tt, collapse = " - ")
  })

first_specs  <- unique(pop_size_ests.f$spec)[c(1, 2, 3)]
second_specs <- unique(pop_size_ests.f$spec)[c(4, 5)]

# "#E6AB02" "#A6761D" "#7570B3" "#E7298A" "#1B9E77"

pop_size_ests.f %>% 
  filter(spec %in% c(
    unique(pop_size_ests.f$spec)[c(4, 5)]
  )) %>% {
    ggplot(., aes(capture_date, mid)) + 
      geom_errorbar(aes(ymin = lwr, ymax = upr, colour = spec), width = 0.2, size = 0.3) +
      geom_errorbar(aes(ymin = lwr_n, ymax = upr_n, colour = spec), width = 0, size = 0.8) +
      scale_color_manual(
      # values  = RColorBrewer::brewer.pal(8,"Dark2")[c(6,7,3,4,1)]
        values  = RColorBrewer::brewer.pal(8,"Dark2")[c(4,1)]
      , labels  = spec_labs[c(4, 5)]
      , name    = "Species"
      ) +
      xlab("Date") +
      ylab("Population Estimate") +
      theme(
        axis.text.x = element_text(size = 10)
      , axis.text.y = element_text(size = 10)
      , strip.text.x = element_text(size = 10)
      , ) +
      facet_wrap(~ppop, scales = "free", ncol = 2) +
      scale_y_log10()
  }

#   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
#   geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.2) +
#   geom_line() +
#   geom_point(aes(capture_date, capt_per_day), colour = "firebrick3", size = 3) +
#      scale_x_continuous(
#        breaks = x_labs$date_fac
#      , labels = x_labs$capture_date
#        ) +
# axis.text.x = element_text(angle = 300, hjust = 0, size = 10)


#### Figure SX ------
## Figure SX: 
 ## Width of CI for Average Survival against Caps per ind
 ## Width of CI for Bd effect Swabs per individual

spec_labs <- c(
  expression(italic("Ambystoma cingulatum"))
  , expression(italic("Anaxyrus boreas"))
  , expression(italic("Notophthalmus viridescens"))
  , expression(italic("Pseudacris maculata"))
  , expression(italic("Rana spp."))
)

ggs.1 <- bd.est.gg %>% {
  ggplot(., aes(swabs_per_ind, CI_width, colour = Species)) + 
    geom_point(aes(size = caps_per_ind)) +
    scale_size_continuous(name = "Average Number
of Captures
Per Individual") +
    scale_color_manual(
      values = RColorBrewer::brewer.pal(8,"Dark2")[c(6,7,3,4,1)]
      , labels  = spec_labs
    ) +
    xlab("Average Number of Swabs Per Individual") +
    ylab("Width of 95% CI for 
Bd Effect (logit scale)") +
    theme(axis.text.y = element_text(size = 11)
          , legend.key.size = unit(0.3, "cm")
          , axis.text.x = element_text(size = 11)
          , axis.title.x = element_text(size = 13)
          , axis.title.y = element_text(size = 13)
          , legend.text = element_text(size = 9)
          , legend.title = element_text(size = 10)
          , legend.text.align = 0)
}

int.est.gg <- rbind(
  ind_mehg_preds[[2]]$int.est.gg
  , pop_mehg_preds[[2]]$int.est.gg
) %>% arrange(desc(mid)) %>% 
  mutate(
    pop_spec = factor(pop_spec, levels = unique(pop_spec))
    , Species  = factor(Species, levels = unique(Species) %>% as.character() %>% sort())
  )

spec_labs <- c(
  expression(italic("Ambystoma cingulatum"))
  , expression(italic("Anaxyrus boreas"))
  , expression(italic("Notophthalmus viridescens"))
  , expression(italic("Pseudacris maculata"))
  , expression(italic("Rana spp."))
)

ggs.2 <- int.est.gg %>% {
  ggplot(., aes(caps_per_ind, CI_width, colour = Species)) + 
    geom_point(size = 2.5) +
    scale_color_manual(
      values = RColorBrewer::brewer.pal(8,"Dark2")[c(6,7,3,4,1)]
      , labels  = spec_labs
    ) +
    xlab("Average Number of Captures Per Individual") +
    ylab("Width of 95% CI for Average 
Survival (probability scale)") +
    theme(axis.text.y = element_text(size = 11)
          , legend.key.size = unit(0.3, "cm")
          , axis.text.x = element_text(size = 11)
          , axis.title.x = element_text(size = 13)
          , axis.title.y = element_text(size = 13)
          , legend.text = element_text(size = 9)
          , legend.title = element_text(size = 10)
          , legend.text.align = 0)
}

grid.arrange(ggs.1, ggs.2, ncol = 1)

