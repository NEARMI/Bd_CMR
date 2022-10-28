###################################################
## Plot diagnostics for a joint population model ##
###################################################

####
## Data loading, sample and summary extraction
####

## stan.fit.samples <- readRDS("fits/CMR_multiple_populations_mehg_ssp_alt_p_len_RANA_2022-10-26_samples.Rds")

if (plot_from == "fit" | plot_from == "saved_model") {

 if (plot_from == "fit") {
  stan.fit.summary <- summary(stan.fit)[[1]]
  stan.fit.samples <- extract(stan.fit)
 } else if (plot_from == "saved_model") {
  stan.fit         <- readRDS(saved_model)
  stan.fit.summary <- summary(stan.fit[[1]])[[1]]
  stan.fit.samples <- extract(stan.fit[[1]])
 }
  
## Can easily have problems with memory, so subset to just the parameters needed for plotting
needed_entries <- c(
#   "beta_offseason_int"
# , "beta_offseason_bd"
# , "beta_offseason_len"
# , "beta_offseason_mehg"
# , "beta_offseason_mehg_bd"
   "beta_phi"
 , "z_r" 
 , "beta_inseason_int"
 , "inseason_pop"
 , "bd_ind"
 , "beta_p_int"
 , "p_pop"
 , "p_day_dev"
 , "beta_p_slope"
 , "pop_size"
)

stan.fit.samples <- stan.fit.samples[needed_entries]
  
} else if (plot_from == "saved_samples") {
  
 cleaned_output.temp <- readRDS(saved_samples)
 stan.fit.summary    <- cleaned_output.temp[[1]]
 stan.fit.samples    <- cleaned_output.temp[[2]]
  
}

## Add some list entries so I don't have to change all of the code below
stan.fit.samples$beta_offseason_int     <- stan.fit.samples$beta_phi[, 1:3]
stan.fit.samples$beta_offseason_bd      <- stan.fit.samples$beta_phi[, 4]
stan.fit.samples$beta_offseason_len     <- stan.fit.samples$beta_phi[, 5]
stan.fit.samples$beta_offseason_mehg    <- stan.fit.samples$beta_phi[, 6]
stan.fit.samples$beta_offseason_mehg_bd <- stan.fit.samples$beta_phi[, 7]

## extract the species and population names for this fit
these_specs <- unique(capt_history.phi$Species)
n_specs     <- length(these_specs)
these_pops  <- unique(capt_history.phi$pop_spec)
these_sexes <- c("F", "M")    ## skipping U for now

## Very non-dynamic... Write out the Species and Populations in a full way for beautified plots
 ## Use either just the population, or Species and Population
plot_labs <- "Pop" # "Spec-Pop" 
# source("plotting_facet_names.R")
# source("scr_facet_names.R")

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
      (sweep(stan.fit.samples$beta_offseason_bd, 2, spec_sex_mm.t[1:n_spec], `*`) %>% rowSums()) + 
        stan.fit.samples$z_r[, 2, spec_sex$pop[k]]
      ) * pred.vals$bd[j] +
    (
      (sweep(stan.fit.samples$beta_offseason_len, 2, spec_sex_mm.t[1:n_spec], `*`) %>% rowSums()) +
        stan.fit.samples$z_r[, 3, spec_sex$pop[k]]
      ) * pred.vals$len[j] +
    (
      sweep(stan.fit.samples$beta_offseason_mehg, 2, spec_sex_mm.t[1:n_spec], `*`) %>% rowSums()
      ) * pred.vals$mehg[j]
 )
    } else {
 pred.est[j, ] <- plogis(
    (sweep(stan.fit.samples$beta_offseason_int, 2, spec_sex_mm.t, `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 1, spec_sex$pop[k]] + 
    (
      (sweep(stan.fit.samples$beta_offseason_bd, 2, spec_sex_mm.t[1:n_spec], `*`) %>% rowSums()) + 
        stan.fit.samples$z_r[, 2, spec_sex$pop[k]]
      ) * pred.vals$bd[j] +
    (
      (sweep(stan.fit.samples$beta_offseason_len, 2, spec_sex_mm.t[1:n_spec], `*`) %>% rowSums()) +
        stan.fit.samples$z_r[, 3, spec_sex$pop[k]]
      ) * pred.vals$len[j] +
    (
      sweep(stan.fit.samples$beta_offseason_mehg, 2, spec_sex_mm.t[1:n_spec], `*`) %>% rowSums() +
        stan.fit.samples$z_r[, 4, spec_sex$pop[k]]
      ) * pred.vals$mehg[j] +
    (
      sweep(stan.fit.samples$beta_offseason_mehg_bd, 2, spec_sex_mm.t[1:n_spec], `*`) %>% rowSums() +
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
    (stan.fit.samples$beta_offseason_len + stan.fit.samples$z_r[, 3, spec_sex$pop[k]]) * pred.vals$len[j] +
    (stan.fit.samples$beta_offseason_mehg + stan.fit.samples$z_r[, 4, spec_sex$pop[k]]) * pred.vals$mehg[j]
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
  pop  = plyr::mapvalues(pop , from = unique(pred.vals.gg.f$pop) , to = as.character(these_pops))
, sex  = plyr::mapvalues(sex , from = unique(pred.vals.gg.f$sex) , to = these_sexes)
)

if (n_specs > 1) {
pred.vals.gg.f %<>% mutate(spec = plyr::mapvalues(spec, from = unique(pred.vals.gg.f$spec)
  , to = as.character(these_specs)))
pred.vals.gg.f %<>% mutate(pop =  as.factor(pop), spec = as.factor(spec))
} else {
pred.vals.gg.f %<>% mutate(pop =  as.factor(pop))
}

pred.vals.gg <- pred.vals.gg.f; rm(pred.vals.gg.f)

####
## Plotting
####

gg1 <- pred.vals.gg %>% mutate(
  pop = plyr::mapvalues(pop, from = unique(pred.vals.gg$pop)
  , to = spec_pop_plot_labels
)) %>% filter(sex == "M", len == 0, mehg == 0) %>% {
  ggplot(., aes(bd, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr
    #  , fill = spec, colour = spec
      ), alpha = 0.3) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n
    #  , fill = spec, colour = spec
      ), alpha = 0.3) +
    geom_line(aes(
    #  colour = spec
      ), size = 1) + 
#    scale_colour_brewer(name = "Species", palette = "Dark2", labels = spec_labs) +
 #   scale_fill_brewer(name = "Species", palette = "Dark2", labels = spec_labs) +
   # scale_x_continuous(breaks = c(0, 3, 6, 9, 12)) +
    scale_x_continuous(breaks = c(-1.40, -0.75, 0, 0.75, 1.40)) +
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
    xlab("Bd Load (scaled)") +
    ylab("Between-Season Survival")
}

gg2 <- pred.vals.gg %>% mutate(
  pop = plyr::mapvalues(pop, from = unique(pred.vals.gg$pop)
  , to = spec_pop_plot_labels
)) %>% filter(sex == "M", bd == 0, mehg == 0) %>% {
  ggplot(., aes(len, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr
    #  , fill = spec, colour = spec
      ), alpha = 0.3) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n
     # , fill = spec, colour = spec
      ), alpha = 0.3) +
    geom_line(aes(
    #  colour = spec
      ), size = 1) + 
    scale_colour_brewer(name = "Species", palette = "Dark2"
      , labels = spec_labs) +
    scale_fill_brewer(name = "Species", palette = "Dark2"
      , labels = spec_labs) +
    scale_x_continuous(breaks = c(-1.75, -0.75, 0, 0.75, 1.75)) +
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
    xlab("Length (scaled)") +
    ylab("Between-Season Survival")
}

gg3 <- pred.vals.gg %>% mutate(
  pop = plyr::mapvalues(pop, from = unique(pred.vals.gg$pop)
  , to = spec_pop_plot_labels
)) %>% filter(sex == "M", bd == 0, len == 0) %>% {
  ggplot(., aes(mehg, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr
     # , fill = spec, colour = spec
      ), alpha = 0.3) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n
     # , fill = spec, colour = spec
      ), alpha = 0.3) +
    geom_line(aes(
    #  colour = spec
      ), size = 1) + 
    scale_colour_brewer(name = "Species", palette = "Dark2"
      , labels = spec_labs) +
    scale_fill_brewer(name = "Species", palette = "Dark2"
      , labels = spec_labs) +
    scale_x_continuous(breaks = c(-1.75, -0.75, 0, 0.75, 1.75)) +
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
    xlab("MeHg (scaled)") +
    ylab("Between-Season Survival")
}

if (fit_ind_mehg) {
  
pred.vals.gg %>% mutate(
  pop = plyr::mapvalues(pop, from = unique(pred.vals.gg$pop)
  , to = spec_pop_plot_labels
)) %>% filter(sex == "M", len == 0, mehg %in% unique(pred.vals.gg$mehg)[c(2, 4, 6, 8, 10)]) %>% 
 # filter(pop == "MT - Jones Pond", spec == "ANBO") %>% 
    {
  ggplot(., aes(bd, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr
      , fill = pop
      ), alpha = 0.3) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n
      , fill = pop
      ), alpha = 0.3) +
    geom_line(size = 1) + 
    scale_x_continuous(breaks = c(-1.75, -0.75, 0, 0.75, 1.75)) +
    facet_wrap(~mehg) +
    theme(
    strip.text.x    = element_text(size = 11)
  , axis.text.y     = element_text(size = 12)
  , axis.text.x     = element_text(size = 11)
  , legend.key.size = unit(0.65, "cm")
  , legend.text     = element_text(size = 11)
  , legend.title    = element_text(size = 13)
  , legend.text.align = 0
    ) +
    xlab("MeHg (scaled)") +
    ylab("Between-Season Survival")
}
  
}

### CI on Intercept and Bd effect

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

bd.mehg <- matrix(
  data = 0
, nrow = dim(stan.fit.samples$z_r)[1]
, ncol = dim(stan.fit.samples$z_r)[3]
)

## Derived quantities here being average survival at the mean of all continuous predictors and
 ## the effect of Bd
for (j in 1:ncol(bd.est)) {
  if (n_spec > 1) {
  if (multi_spec_red) {
   int.est[, j] <- stan.fit.samples$beta_offseason_int[, 1] + stan.fit.samples$z_r[, 1, j]
   bd.est[, j]  <- stan.fit.samples$beta_offseason_bd + stan.fit.samples$z_r[, 2, j]
  } else {
  int.est[, j] <- (sweep(stan.fit.samples$beta_offseason_int[, 1:n_spec], 2
    , spec_sex_mm[spec_sex_mm$spec == spec_in_pop[j] & spec_sex_mm$sex == "M", 1:n_spec] %>% as.matrix() %>% c(), `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 1, j]
   bd.est[, j] <- (sweep(stan.fit.samples$beta_offseason_bd, 2
    , spec_sex_mm[spec_sex_mm$spec == spec_in_pop[j] & spec_sex_mm$sex == "M", 1:n_spec] %>% as.matrix() %>% c(), `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 2, j]
  }
  if (fit_ind_mehg) {
    bd.mehg[, j] <- (sweep(stan.fit.samples$beta_offseason_mehg_bd[, 1:n_spec], 2
    , spec_sex_mm[spec_sex_mm$spec == spec_in_pop[j] & spec_sex_mm$sex == "M", 1:n_spec] %>% as.matrix() %>% c(), `*`) %>% rowSums()) +
      stan.fit.samples$z_r[, 5, j]
  }
  } else {
   int.est[, j] <- stan.fit.samples$beta_offseason_int[, 1] + stan.fit.samples$z_r[, 1, j]
   bd.est[, j]  <- stan.fit.samples$beta_offseason_bd + stan.fit.samples$z_r[, 2, j] 
  if (fit_ind_mehg) {
   bd.mehg[, j] <- stan.fit.samples$beta_offseason_mehg_bd + stan.fit.samples$z_r[, 5, j]
  }
    
  }
}

int.est        <- reshape2::melt(int.est)
names(int.est) <- c("Sample", "Population", "Value")
bd.est        <- reshape2::melt(bd.est)
names(bd.est) <- c("Sample", "Population", "Value")

if (!fit_ind_mehg) {

## Not dynamic, needs to get manually updated if the species change
int.est %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = c(
 #   "Ambystoma cingulatum"
    "Anaxyrus boreas"
  , "Pseudacris maculata"
  #, "Notophthalmus viridescens"
  , "Rana spp."
  ))
)

bd.est %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = c(
#    "Ambystoma cingulatum"
 #   "Anaxyrus boreas"
#  , "Pseudacris maculata"
#  , "Notophthalmus viridescens"
    "Rana spp."
  ))
)

} else {
  
## Not dynamic, needs to get manually updated if the species change
int.est %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = c(
 #   "Anaxyrus boreas"
 # , "Pseudacris maculata"
    "Rana spp."
  ))
)

bd.est %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = c(
#    "Anaxyrus boreas"
#  , "Pseudacris maculata"
    "Rana spp."
  ))
) 
  
}

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
  pop_spec = unique(capt_history$pop_spec)
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
  pop_spec = unique(capt_history$pop_spec)
) %>% mutate(
  CI_width = upr - lwr
)

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
  pop_spec = unique(capt_history$pop_spec)
) %>% mutate(
  CI_width = upr - lwr
)

if (fit_ind_mehg) {
bd.mehg        <- reshape2::melt(bd.mehg)
names(bd.mehg) <- c("Sample", "Population", "Value")
bd.mehg        %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = c(
 #   "Anaxyrus boreas"
 # , "Pseudacris maculata" 
    "Rana spp."
  ))
)

## save full un-summarized for some text description
bd.mehg.full <- bd.mehg

bd.mehg %<>% group_by(Population, Species) %>% summarize(
    lwr   = quantile(Value, 0.025)
  , lwr_n = quantile(Value, 0.200)
  , mid   = quantile(Value, 0.500)
  , upr_n = quantile(Value, 0.800)
  , upr   = quantile(Value, 0.975)
) %>% ungroup() %>% mutate(
  pop_spec = unique(capt_history$pop_spec)
) %>% mutate(
  CI_width = upr - lwr
)
}

## Calculate a number of metrics of effort for some coarse correlations between effort and CI width
num_samples  <- sampling %>% group_by(pop_spec) %>% summarize(n_dates = n_distinct(CaptureDate)) %>%
  arrange(desc(n_dates)) %>% mutate(pop_spec = factor(pop_spec, levels = unique(pop_spec)))

recaps <- capt_history %>% group_by(pop_spec, Mark) %>% 
  filter(captured == 1) %>%
  summarize(recaps = sum(captured)) %>% 
  mutate(recaptured = ifelse(recaps > 1, 1, 0)) %>% 
  ungroup(Mark) %>%
  summarize(
    inds_capt  = n_distinct(Mark)
  , recapt_ind = sum(recaptured)
  , caps_per_ind = mean(recaps)) %>% 
  ungroup() 

n_swabs <- capt_history %>% group_by(pop_spec, Mark) %>% 
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

int.est %>% {
  ggplot(., aes(caps_per_ind, CI_width)) + geom_point() + 
    scale_x_log10() +
    xlab("Number of Captures Per Individual") +
    ylab("95% CI Width")
}

bd.est %<>% left_join(., num_samples) %>% left_join(., recaps) %>% left_join(., n_swabs)

bd.est %<>% mutate(Population = spec_pop_plot_labels)

gg.cor.1 <- ggplot(bd.est, aes(caps_per_ind, CI_width)) + 
    geom_point(size = 3, aes(colour = pop_spec, shape = Species)) + 
    xlab("Average Number of Captures Per Individual") +
    ylab("Width of 95% CI for Bd Survival Effect") +
  scale_color_manual(name = "Population"
    , labels = bd.est$Population
    , values = colorRampPalette(RColorBrewer::brewer.pal(name="Paired", n = 12))(20)) +
  scale_shape_manual(values = c(2, 15, 16, 0, 17, 1)) +
theme(legend.key.size = unit(0.65, "cm")
 # , legend.position = "none"
  )

gg.cor.2 <- ggplot(bd.est, aes(swabs_per_ind, CI_width)) + 
    geom_point(size = 3, aes(colour = pop_spec, shape = Species)) + 
    xlab("Average Number of Bd Measures Per Individual") +
    ylab("Width of 95% CI for Bd Survival Effect") +
  scale_color_manual(name = "Population"
    , labels = bd.est$Population
    , values = colorRampPalette(RColorBrewer::brewer.pal(name="Paired", n = 12))(20)) +
  scale_shape_manual(values = c(2, 15, 16, 0, 17, 1)) +
theme(legend.position = "none")

gridExtra::grid.arrange(gg.cor.1, gg.cor.2, ncol = 1)

int.est.gg <- int.est %>% mutate(Population = plyr::mapvalues(Population
  , from = unique(Population)
  , to = spec_pop_plot_labels)
) %>% arrange(desc(mid)) %>% 
  mutate(pop_spec = factor(pop_spec, levels = pop_spec)) 

int.est.gg %>% {
  ggplot(., aes(mid, pop_spec, colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_point(size = 2.5) +
    scale_color_brewer(
      palette = "Dark2"
    , labels = spec_labs) +
    scale_y_discrete(labels = int.est.gg$Population) +
    xlab("Survival at mean conditions (Bd load, Length, MeHg)") +
    ylab("Population") +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13)
      , legend.text.align = 0)
}

bd.est.gg <- bd.est %>% mutate(Population = plyr::mapvalues(Population
  , from = unique(Population)
  , to = spec_pop_plot_labels
)) %>% arrange(desc(mid)) %>% 
  mutate(pop_spec = factor(pop_spec, levels = pop_spec)) 

bd.est.gg %>% {
  ggplot(., aes(mid, pop_spec, colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_point(size = 2.5) +
    scale_color_brewer(
      palette = "Dark2"
    , labels = spec_labs) +
    scale_y_discrete(labels = bd.est.gg$Population) +
    xlab("Effect of Bd on Survival (Logit Scale)") +
    scale_x_continuous(breaks = c(-2.5, -1.5, -0.5, 0, 0.5, 1.5, 2.5)) +
    ylab("Population") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , axis.text.x = element_text(size = 10)
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13)
      , legend.text.align = 0)
}

## ****** Manually combining the estimates from the MeHg and main fit for a supplemental figure

bd.est.gg.mehg_model %>% 
  rename(lwr.m = lwr, lwr_n.m = lwr_n, mid.m = mid, upr_n.m = upr_n, upr.m = upr) %>% 
  dplyr::select(-CI_width) %>% 
  left_join(.
  , bd.est.gg %>% dplyr::select(-CI_width)
  ) %>% {
  ggplot(., aes(mid, pop_spec, colour = Species)) + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3, size = 0.8, linetype = "dotted"
      , position = position_nudge(y = -0.15)) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5
      , position = position_nudge(y = -0.15)) +
    geom_point(size = 2
      , position = position_nudge(y = -0.15)) +
      
    geom_errorbarh(aes(xmin = lwr.m, xmax = upr.m), height = 0.3, size = 0.8
      , position = position_nudge(y = 0.15)) +
    geom_errorbarh(aes(xmin = lwr_n.m, xmax = upr_n.m), height = 0.0, size = 1.5
      , position = position_nudge(y = 0.15)) +
    geom_point(aes(mid.m, pop_spec), size = 2
      , position = position_nudge(y = 0.15)) +
      
    scale_y_discrete(labels = bd.est.gg.mehg_model$Population) +
    scale_color_manual(values = c("#D95F02", "#7570B3", "#66A61E")) +
    xlab("Bd-MeHg Interactive Effect (logit scale)") +
    ylab("Population") +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13))
  }

## ******

##########

## Sample plot as the above, but converting the effect to the probability scale

bd.est2 <- int.est2 %>% ungroup() %>% mutate(Population = plyr::mapvalues(Population
  , from = unique(Population)
  , to = spec_pop_plot_labels)
) %>% rename(mid_p = mid) %>% left_join(.
  , 
  bd.est %>% mutate(Population = plyr::mapvalues(Population
  , from = unique(Population)
  , to = spec_pop_plot_labels
))
  ) %>% mutate(
    lwr   = plogis(mid_p) / plogis(mid_p - lwr)    
  , lwr_n = plogis(mid_p) / plogis(mid_p - lwr_n)
  , mid   = plogis(mid_p) / plogis(mid_p - mid)    
  , upr_n = plogis(mid_p) / plogis(mid_p - upr_n)  
  , upr   = plogis(mid_p) / plogis(mid_p - upr)    
  ) %>% arrange(desc(mid)) %>% 
  mutate(pop_spec = factor(pop_spec, levels = pop_spec)) 

bd.est2.gg <- bd.est2 %>% mutate(Population = plyr::mapvalues(Population
  , from = unique(Population)
  , to = spec_pop_plot_labels
)) %>% arrange(desc(mid)) %>% 
  mutate(pop_spec = factor(pop_spec, levels = pop_spec)) 

bd.est2.gg %>% {
  ggplot(., aes(mid, pop_spec, colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_point(size = 2.5) +
    scale_color_brewer(
      palette = "Dark2"
    , labels = spec_labs
      ) +
    scale_y_discrete(labels = bd.est2.gg$Population) +
    xlab("Surival probability at one SD above the mean Bd load 
relative to survival at  the mean Bd load") +
    scale_x_log10(breaks = c(0.1, 0.40, 0.70, 1.0, 1.5, 3.0)) +
    ylab("Population") +
    geom_vline(xintercept = 1, linetype = "dashed", size = 0.4) +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , axis.text.x = element_text(size = 10)
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13)
      , legend.text.align = 0)
}

if (fit_ind_mehg) {

## Different populations so the previously stored vector of pop spec names won't work
bd.mehg %<>% mutate(Population = plyr::mapvalues(Population, from = unique(Population), to = spec_pop_plot_labels)) %>% 
    arrange(desc(mid)) %>% mutate(pop_spec = factor(pop_spec, levels = pop_spec)) 
  
bd.mehg %>% {
  ggplot(., aes(mid, pop_spec, colour = Species)) + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_point(size = 2) +
    scale_y_discrete(labels = bd.mehg$Population) +
    scale_color_manual(values = c("#D95F02", "#7570B3", "#66A61E")) +
    xlab("Bd-MeHg Interactive Effect (logit scale)") +
    ylab("Population") +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13))
}

## Two side by side
bd.mehg.main %>% rename(lwr.m = lwr, lwr_n.m = lwr_n, mid.m = mid, upr_n.m = upr_n, upr.m = upr) %>% dplyr::select(-CI_width) %>% 
  left_join(.
  , bd.mehg %>% dplyr::select(-CI_width)
  ) %>% {
  ggplot(., aes(mid, pop_spec, colour = Species)) + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3, size = 0.8, linetype = "dotted"
      , position = position_nudge(y = -0.15)) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5
      , position = position_nudge(y = -0.15)) +
    geom_point(size = 2
      , position = position_nudge(y = -0.15)) +
      
    geom_errorbarh(aes(xmin = lwr.m, xmax = upr.m), height = 0.3, size = 0.8
      , position = position_nudge(y = 0.15)) +
    geom_errorbarh(aes(xmin = lwr_n.m, xmax = upr_n.m), height = 0.0, size = 1.5
      , position = position_nudge(y = 0.15)) +
    geom_point(aes(mid.m, pop_spec), size = 2
      , position = position_nudge(y = 0.15)) +
      
    scale_y_discrete(labels = bd.mehg.main$Population) +
    scale_color_manual(values = c("#D95F02", "#7570B3", "#66A61E")) +
    xlab("Bd-MeHg Interactive Effect (logit scale)") +
    ylab("Population") +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13))
  }
  
## quick check on posterior less than 0
  length(which(
    (bd.mehg.full %>% filter(Population == 8))$Value < 0
  )) / 3000

}

## Continuing with all of the other beta estimates

beta_est <- stan.fit.summary[grep("beta", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%')

param_names <- apply(
  matrix(beta_est$Var1 %>% as.character())
, 1
, FUN = function(x) strsplit(x, "[[]")[[1]][1]
)

beta_est %<>% 
  mutate(Var1 = as.character(Var1)) %>% 
  mutate(Var1 = plyr::mapvalues(Var1, from = unique(beta_est$Var1), to = param_names)) %>%
  rename(params = Var1) %>%
  group_by(params) %>%
  mutate(param_lev = seq(n())) %>% 
  relocate(param_lev, .after = params) %>%
  mutate(param_lev = as.character(param_lev))

## Also pretty non-dynamic, not too sure what to do here to make this better.
 ## I guess name my parameters better?
if (n_specs > 1) {
beta_est.int <- beta_est %>% filter(
  params %in% c("beta_offseason_int", "beta_p_int", "beta_len")
) 
beta_est.spec <- beta_est %>% filter(
  params %in% c("beta_bd_spec", "beta_inseason", "beta_mehg_spec")
)
beta_est.spec %<>% mutate(
  param_lev = plyr::mapvalues(param_lev, from = unique(beta_est.spec$param_lev)
    , to = c(these_specs, "F", "U")
    )
)
beta_est.slopes <- beta_est %>% filter(
  params %in% c("beta_bd_temp", "beta_bd_len", "beta_p_slope", "beta_mehg_drawdown")
)

## Very specific few estimates for a manuscript coefficient plot
beta_est.len_mehg <- beta_est %>% filter(
  params %in% c("beta_offseason_len", "beta_offseason_mehg")
)
beta_est.mehg <- beta_est %>% filter(
  params %in% c("beta_offseason_mehg")
)

} else {
beta_est.int <- beta_est %>% filter(
  params %in% c("beta_offseason_int", "beta_inseason_int", "beta_p_int"
  #  , "beta_len"
    , "beta_mehg_int")
) 
beta_est.slopes <- beta_est %>% filter(
  params %in% c("beta_bd_temp", "beta_bd_len", "beta_p_slope", "beta_mehg_drawdown")
)
}

gg3 <- beta_est.int %>% 
  mutate(params = plyr::mapvalues(params
    , from = unique(beta_est.int$params)
    , to = c("Offseason Average Survival", "Average Detection Probability")
    )) %>% {
    ggplot(., aes(mid, param_lev)) + geom_point() +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3) +
      facet_wrap(~params) +
      ylab("Parameter Level") + 
      scale_y_discrete(labels = c(spec_labs, "Sex: Female", "Sex: Unknown")) +
      geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
      theme(
        axis.text.y  = element_text(size = 12)
      , strip.text.x   = element_text(size = 12)) +
      xlab("Coefficient Estimate")
}

if (n_specs > 1) {
gg4 <- beta_est.spec %>%
  mutate(params = plyr::mapvalues(params
    , from = unique(beta_est.spec$params)
    , to = c("Species Average Survival Between
Primary Periods Within A Year")
    )) %>% {
    ggplot(., aes(mid, param_lev)) + geom_point() +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3) +
      facet_wrap(~params) +
      ylab("Parameter Level") +
      scale_y_discrete(labels = c(spec_labs)) +
      geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
      theme(axis.text.y = element_text(size = 12)) +
      xlab("Coefficient Estimate")
}
}

gg5 <- beta_est.slopes %>% 
  mutate(params = plyr::mapvalues(params
    , from = unique(beta_est.slopes$params)
    , to = c("Coefficients Affecting Detection")
    )) %>% {
    ggplot(., aes(mid, param_lev)) + geom_point() +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3) +
      facet_wrap(~params) +
      ylab("Parameter") +
      scale_y_discrete(labels = c("Drawdown", "Vegetation")) +
      geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
      theme(axis.text.y = element_text(size = 12)) +
      xlab("Coefficient Estimate")
}

beta_est.len_mehg %>% 
  mutate(params = plyr::mapvalues(params
    , from = unique(beta_est.len_mehg$params)
    , to = c("Length", "MeHg Concentration")
    )) %>% {
    ggplot(., aes(mid, param_lev)) + 
      geom_point() +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2) +
      facet_wrap(~params) +
      ylab("Parameter") +
      scale_y_discrete(labels = spec_labs) +
      geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
      theme(axis.text.y = element_text(size = 12)
      , axis.text.x = element_text(size = 12)) +
      xlab("Coefficient Estimate")
}

## individual bd deviates

stan.ind_pred_var <- stan.fit.samples$bd_ind %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value) %>%
  group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% arrange(mid) 

capt_history.temp <- capt_history.p %>% group_by(Mark) %>% slice(1) %>% dplyr::select(pop_spec, Mark) %>% rename(ind = Mark)

stan.ind_pred_var %<>% left_join(., capt_history.temp)

stan.ind_pred_var %<>% ungroup() %>% group_by(pop_spec) %>% arrange(mid) %>%
  mutate(est_rank = seq(n()))

ind_bd_meas <- capt_history.p %>%
  filter(swabbed == 1) %>%
  group_by(Mark, pop_spec) %>% 
  summarize(
   mean_bd = mean(log_bd_load)
  ) %>% ungroup() %>% group_by(pop_spec) %>% 
  arrange(mean_bd) %>% mutate(real_rank = seq(n())) %>% 
  rename(ind = Mark)

stan.ind_pred_var <- left_join(stan.ind_pred_var, ind_bd_meas)

gg6 <- stan.ind_pred_var %>% group_by(pop_spec) %>% mutate(max_rank = max(est_rank)) %>%
  filter(est_rank < 30 | est_rank > (max_rank - 30)) %>% arrange(desc(mid)) %>% 
  mutate(ind = factor(ind, levels = ind)) %>% {
  ggplot(., aes(mid, ind)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3) + 
        geom_point() +
        theme(axis.text.y = element_text(size = 8)) +
        xlab("") +
        ylab("Individual") +
        xlab("Individual bd deviate") +
      facet_wrap(~pop_spec, scales = "free")
}

gg7 <- stan.ind_pred_var %>% {
  ggplot(., aes(real_rank, est_rank)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    xlab("Real Bd Rank") +
    ylab("Estimated Bd Rank") +
          facet_wrap(~pop_spec, scales = "free")
}

## detection 

each_pop_day <- capt_history.p %>% group_by(date_fac) %>% slice(1) %>%
  dplyr::select(pop_spec, capture_date, date_fac, veg_cont, drawdown_cont, Species) %>%
  mutate(pop_spec = as.numeric(pop_spec))

spec_sex_mm.p <- spec_sex_mm %>% filter(sex == "M")

pred.est.p <- matrix(data = 0, nrow = nrow(each_pop_day), ncol = dim(stan.fit.samples[[1]])[1])

for (i in 1:nrow(each_pop_day)) {

if (n_specs > 1) {
  
  which_temp_spec <- which(these_specs == each_pop_day$Species[i])
  spec_sex_mm.p.t <- spec_sex_mm.p %>% filter(spec == which_temp_spec) %>%
    dplyr::select(-spec, -sex) %>% as.matrix()
  
  pred.est.p[i, ] <- plogis(
    (sweep(stan.fit.samples$beta_p_int, 2, spec_sex_mm.p.t, `*`) %>% rowSums()) +
     stan.fit.samples$p_pop[, each_pop_day$pop_spec[i]] +
     stan.fit.samples$p_day_dev[, each_pop_day$date_fac[i]] +
     stan.fit.samples$beta_p_slope[, 1] * each_pop_day$drawdown_cont[i] +
     stan.fit.samples$beta_p_slope[, 2] * each_pop_day$veg_cont[i]
  )
  
} else {
  
  spec_sex_mm.p.t <- spec_sex_mm.p[1, 1:3] %>% as.matrix()
  
  pred.est.p[i, ] <- plogis(
    (sweep(stan.fit.samples$beta_p_int, 2, spec_sex_mm.p.t, `*`) %>% rowSums()) +
     stan.fit.samples$p_pop[, each_pop_day$pop_spec[i]] +
     stan.fit.samples$p_day_dev[, each_pop_day$date_fac[i]] +
     stan.fit.samples$beta_p_slope[, 1] * each_pop_day$drawdown_cont[i] +
     stan.fit.samples$beta_p_slope[, 2] * each_pop_day$veg_cont[i]
  )
  
}

}

each_pop_day %<>% cbind(., as.data.frame(pred.est.p))

each_pop_day %<>% ungroup() %>% 
  pivot_longer(., c(-pop_spec, -capture_date, -date_fac, -veg_cont, -drawdown_cont, -Species)
  , names_to = "iter", values_to = "est")
  
each_pop_day.gg <- each_pop_day %>%
  group_by(pop_spec, capture_date, date_fac, veg_cont, drawdown_cont, Species) %>%
  summarize(
    lwr   = quantile(est, 0.025)
  , lwr_n = quantile(est, 0.200)
  , mid   = quantile(est, 0.500)
  , upr_n = quantile(est, 0.800)
  , upr   = quantile(est, 0.975)
  ) %>% ungroup() 

each_pop_day.gg %<>% mutate(pop_spec = plyr::mapvalues(pop_spec, from = unique(each_pop_day.gg$pop_spec), to = as.character(these_pops)))

each_pop_day.gg %<>%
  ungroup() %>%
  group_by(pop_spec) %>%
  arrange(mid) %>%
  mutate(est_val = seq(n())) %>% 
  ungroup() 

day_capts <- capt_history.p %>% ungroup() %>% group_by(pop_spec, date_fac) %>% 
  summarize(tot_capt = sum(captured)) %>% ungroup() %>%
  group_by(pop_spec) %>%
  arrange(tot_capt) %>%
  mutate(ncapt_rank = seq(n()))

each_pop_day.gg %<>% left_join(., day_capts)

gg8 <- each_pop_day.gg %>% 
    group_by(pop_spec) %>%
    filter(capture_date != min(capture_date)) %>%
    mutate(capture_date = as.factor(capture_date)) %>% {
  ggplot(., aes(mid, capture_date)) +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), size = 0.75, height = 0.3) +
    geom_point() +
    xlab("Daily Detection Probability (for a Male)") +
    ylab("Capture outing") + 
    facet_wrap(~pop_spec, scales = "free")
}
  
gg9 <- each_pop_day.gg %>% {
  ggplot(., aes(tot_capt, est_val)) +
    geom_point() +
    xlab("Number of Individuals Captured") +
    ylab("Daily detection estimate (rank)") +
    facet_wrap(~pop_spec, scales = "free")
}

pop_size_ests <- capt_history.p %>% ungroup() %>% 
  group_by(pop_spec, capture_date) %>%
#  group_by(date_fac) %>%
  slice(1) %>% 
  dplyr::select(pop_spec, capture_date) %>% 
  ungroup() %>%
  mutate(capt_per_day = n_capt_per_day_sex %>% rowSums()) %>%
  mutate(
    lwr   = apply(stan.fit.samples$pop_size, 2, FUN = function(x) quantile(x, 0.025))
  , lwr_n = apply(stan.fit.samples$pop_size, 2, FUN = function(x) quantile(x, 0.200))
  , mid   = apply(stan.fit.samples$pop_size, 2, FUN = function(x) quantile(x, 0.500))
  , upr_n = apply(stan.fit.samples$pop_size, 2, FUN = function(x) quantile(x, 0.800))
  , upr   = apply(stan.fit.samples$pop_size, 2, FUN = function(x) quantile(x, 0.975))
    ) %>% ungroup() %>% group_by(pop_spec) %>% 
  mutate(ss = seq(n())) %>% filter(ss != min(ss)) %>% ungroup()

pop_size_ests[pop_size_ests$capt_per_day == 0, 5:9] <- NA

#x_labs <- pop_size_ests %>% 
#  group_by(pop_spec) %>% 
#  filter(ss %in% seq(min(ss), max(ss), by = 5))

pop_size_ests$spec <- apply(pop_size_ests$pop_spec %>% matrix, 1, FUN = function(x) strsplit(x, "[.]")[[1]][1])

pop_size_ests %<>% mutate(
  pop = pop_spec
) 
pop_size_ests %<>% mutate(
  pop  = plyr::mapvalues(pop, from = unique(pop_size_ests$pop)
    , to = spec_pop_plot_labels
    )
, spec = plyr::mapvalues(spec, from = unique(pop_size_ests$spec)
      , to = spec_names)
)

pop_size_ests %<>% unite(pop, pop, spec, sep = " - ")


for (i in 1:(these_pops %>% length())) {

gg.temp <- pop_size_ests %>% 
  filter(pop == unique(pop)[i]) %>% {
    ggplot(., aes(capture_date, mid)) + 
   #   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
   #   geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.2) +
   #   geom_line() +
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, size = 0.3) +
      geom_errorbar(aes(ymin = lwr_n, ymax = upr_n), width = 0, size = 0.8, colour = "dodgerblue4") +
      geom_point(aes(capture_date, capt_per_day), colour = "firebrick3", size = 3) +
      xlab("Date") +
      ylab("Population Estimate") +
#      scale_x_continuous(
#        breaks = x_labs$date_fac
#      , labels = x_labs$capture_date
#        ) +
      theme(
       # axis.text.x = element_text(angle = 300, hjust = 0, size = 10)
        axis.text.x = element_text(size = 12)
        ) +
      facet_wrap(~pop, scales = "free") +
      scale_y_log10() +
      ggtitle("Red Points Show Number of Captures - Lines and Ribbons Show Population Estimates")
  }
  
pdf(paste("plots/pop_sizes/", paste("pop_size", these_pops[i], sep = "_"), ".pdf", sep = ""), onefile = TRUE, width = 9, height = 8)
get("gg.temp") %>% print()
dev.off()

}
  
if (n_specs > 1) {
gglist    <- c("gg1", "gg2", "gg3", "gg4", "gg5", "gg6", "gg7", "gg8", "gg9", "gg10")
} else {
gglist    <- c("gg1", "gg2", "gg3", "gg5", "gg6", "gg7", "gg8", "gg9", "gg10")
}

pdf(paste("plots/", paste("stan_fit_multipop", Sys.Date(), sep = "_"), ".pdf", sep = ""), onefile = TRUE)

for (i in seq(length(gglist))) {
   get(gglist[i]) %>% print()
}
dev.off()
