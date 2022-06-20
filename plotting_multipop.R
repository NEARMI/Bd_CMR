###################################################
## Plot diagnostics for a joint population model ##
###################################################

# model_name <- "fits/stan_fit_multipop_mehg_red_2022-06-14.Rds"

stan.fit         <- readRDS(model_name)
stan.fit.summary <- summary(stan.fit[[1]])[[1]]
stan.fit.samples <- extract(stan.fit[[1]])

## For running this code on local laptop with memory issues
 ## WITH MEMORY ISSUES RUN LINE BY LINE WITH rm() and gc(), do not just execute the script
mem.ish <- TRUE

## Problems with memory, so subset
needed_entries <- c(
  "beta_offseason_int"
, "beta_offseason_bd"
, "beta_offseason_len"
, "beta_offseason_mehg"
, "beta_offseason_mehg_bd"
, "z_r" 
, "beta_inseason_int"
, "bd_ind"
, "beta_p_int"
, "p_pop"
, "p_day_dev"
, "beta_p_slope"
, "pop_size"
)

stan.fit.samples <- stan.fit.samples[needed_entries]

## stan.fit.samples <- readRDS("chains/stan_multipop_samples.Rds")

capt_history.phi <- stan.fit$capt_history.phi
capt_history.p   <- stan.fit$capt_history.p

## extract the species and population names for this fit
these_specs <- unique(capt_history.phi$Species)
n_specs     <- length(these_specs)
these_pops  <- unique(capt_history.phi$pop_spec)
these_sexes <- c("F", "M")    ## skipping U for now

## Re-establish n_capt_per_day_sex for plotting
## number of each sex captured each day 
n_capt_per_day_sex <- capt_history.p %>% group_by(date_fac, Sex) %>% summarize(num_capt = sum(captured)) %>%
    pivot_wider(., date_fac, values_from = num_capt, names_from = Sex) %>% ungroup() %>% dplyr::select(-date_fac) %>%
    as.matrix()
n_capt_per_day_sex[is.na(n_capt_per_day_sex)] <- 0

####
## Plotting Setup
####

spec_in_pop <- (capt_history.p %>% group_by(pop_spec) %>% slice(1))$Species %>% as.numeric()
spec_sex    <- data.frame(
  spec = rep(spec_in_pop, 2)
, sex  = rep(c("F", "M"), each = length(spec_in_pop))
, pop  = rep(seq(length(spec_in_pop)), 2)
)

## rather non-dynamic, can come back and clean this up later possibly
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
spec_sex_mm <- model.matrix(~Sex + pop_spec, capt_history.phi)[, ] %>% as.data.frame() %>% distinct()
spec_sex_mm %<>% mutate(sex = ifelse(SexF == 1 & SexU != 1, "F", (ifelse(SexF != 1 & SexU == 1, "U", "M"))))
}

for (k in 1:nrow(spec_sex)) {
  
pred.vals <- expand.grid(
  bd   = scale(seq(0, 14, by = 1))[, 1]
, len  = seq(-2, 2, by = 0.4)
, mehg = seq(-2, 2, by = 0.4)
)

pred.est <- matrix(data = 0, nrow = nrow(pred.vals), ncol = dim(stan.fit.samples[[1]])[1])

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
 pred.est[j, ] <- plogis(
    (sweep(stan.fit.samples$beta_offseason_int, 2, spec_sex_mm.t, `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 1, spec_sex$pop[k]] + 
    (
      (sweep(stan.fit.samples$beta_offseason_bd, 2, spec_sex_mm.t[1:5], `*`) %>% rowSums()) + 
        stan.fit.samples$z_r[, 2, spec_sex$pop[k]]) * pred.vals$bd[j] +
    (
      (sweep(stan.fit.samples$beta_offseason_len, 2, spec_sex_mm.t[1:5], `*`) %>% rowSums()) +
        stan.fit.samples$z_r[, 3, spec_sex$pop[k]]) * pred.vals$len[j] +
    (sweep(stan.fit.samples$beta_offseason_mehg, 2, spec_sex_mm.t[1:5], `*`) %>% rowSums()) * pred.vals$mehg[j]
 )
  }
} else {
 pred.est[j, ] <- plogis(
    (sweep(stan.fit.samples$beta_offseason_int, 2, spec_sex_mm.t, `*`) %>% rowSums()) +
     stan.fit.samples$offseason_pop[, spec_sex$pop[k]] + 
    (stan.fit.samples$beta_offseason_bd + stan.fit.samples$offseason_pop_bd[, spec_sex$pop[k]]) * pred.vals$bd[j] +
    (stan.fit.samples$beta_offseason_len + stan.fit.samples$offseason_pop_len[, spec_sex$pop[k]]) * pred.vals$len[j] +
     stan.fit.samples$beta_offseason_mehg * pred.vals$mehg[j]
 )
}
  
}

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
print(paste("Through", k, "population:sex", sep = " "))

if (k == 1) {
pred.vals.f <- pred.vals
} else {
pred.vals.f <- rbind(pred.vals.f, pred.vals)
}
}

if (mem.ish) { 

## More memory problems
pred.vals.f.1 <- pred.vals.f[, c(1:5, 6:1005)]
pred.vals.f.2 <- pred.vals.f[, c(1:5, 1006:2005)]
pred.vals.f.3 <- pred.vals.f[, c(1:5, 2006:3005)]

pred.vals.f.1 %<>% pivot_longer(., c(-bd, -len, -mehg, -pop, -sex), names_to = "iter", values_to = "est") 
pred.vals.f.2 %<>% pivot_longer(., c(-bd, -len, -mehg, -pop, -sex), names_to = "iter", values_to = "est")  
pred.vals.f.3 %<>% pivot_longer(., c(-bd, -len, -mehg, -pop, -sex), names_to = "iter", values_to = "est")

pred.vals.f.1 <- rbind(pred.vals.f.1, pred.vals.f.2)
pred.vals.f.1 <- rbind(pred.vals.f.1, pred.vals.f.3)

pred.vals.f   <- pred.vals.f.1

} else {

if (n_specs > 1) {
pred.vals.f %<>% pivot_longer(., c(-bd, -len, -mehg, -pop, -spec, -sex), names_to = "iter", values_to = "est")
} else {
pred.vals.f %<>% pivot_longer(., c(-bd, -len, -mehg, -pop, -sex), names_to = "iter", values_to = "est")  
}

}

if (n_specs > 1) {
pred.vals.gg <- pred.vals.f %>% group_by(bd, len, mehg, spec, pop, sex)
} else {
pred.vals.gg <- pred.vals.f %>% group_by(bd, len, mehg, pop, sex)
}

pred.vals.gg %<>% summarize(
    lwr   = quantile(est, 0.025)
  , lwr_n = quantile(est, 0.200)
  , mid   = quantile(est, 0.500)
  , upr_n = quantile(est, 0.800)
  , upr   = quantile(est, 0.975)
  )

pred.vals.gg %<>% ungroup() %>% mutate(
  pop  = plyr::mapvalues(pop , from = unique(pred.vals.gg$pop) , to = as.character(these_pops))
, sex  = plyr::mapvalues(sex , from = unique(pred.vals.gg$sex) , to = these_sexes)
)

if (n_specs > 1) {
pred.vals.gg %<>% mutate(spec = plyr::mapvalues(spec, from = unique(pred.vals.gg$spec), to = as.character(these_specs)))
pred.vals.gg %<>% mutate(pop =  as.factor(pop), spec = as.factor(spec))
} else {
pred.vals.gg %<>% mutate(pop =  as.factor(pop))
}

gg1 <- pred.vals.gg %>% filter(sex == "M", len == 0, mehg == 0) %>% {
  ggplot(., aes(bd, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = pop), alpha = 0.3) +
    geom_line(aes(colour = pop), size = 1) + 
    scale_colour_discrete() +
    scale_fill_discrete() + {
      if (n_specs > 1) {
        facet_grid(~spec*pop)
      } else {
        facet_grid(~pop)
      }
    } +
    xlab("Bd Load") +
    ylab("Between-Season Survival")
}

gg2 <- pred.vals.gg %>% filter(sex == "M", mehg == 0, bd == 7) %>% {
  ggplot(., aes(len, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = pop), alpha = 0.3) +
    geom_line(aes(colour = pop), size = 1) + 
    scale_colour_discrete() +
    scale_fill_discrete() + {
      if (n_specs > 1) {
        facet_grid(~spec*pop)
      } else {
        facet_grid(~pop)
      }
    } +
    xlab("Length") +
    ylab("Between-Season Survival")
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

for (j in 1:ncol(bd.est)) {
  if (multi_spec_red) {
   int.est[, j] <- stan.fit.samples$beta_offseason_int[, 1] + stan.fit.samples$z_r[, 1, j]
   bd.est[, j]  <- stan.fit.samples$beta_offseason_bd + stan.fit.samples$z_r[, 2, j]
  } else {
  int.est[, j] <- (sweep(stan.fit.samples$beta_offseason_int[, 1:5], 2
    , spec_sex_mm[spec_sex_mm$spec == spec_in_pop[j] & spec_sex_mm$sex == "M", 1:5] %>% as.matrix() %>% c(), `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 1, j]
   bd.est[, j] <- (sweep(stan.fit.samples$beta_offseason_bd, 2
    , spec_sex_mm[spec_sex_mm$spec == spec_in_pop[j] & spec_sex_mm$sex == "M", 1:5] %>% as.matrix() %>% c(), `*`) %>% rowSums()) +
    stan.fit.samples$z_r[, 2, j]
  }
  if (fit_ind_mehg) {
  bd.mehg[, j] <- stan.fit.samples$beta_offseason_mehg_bd + stan.fit.samples$z_r[, 5, j]
  }
}

int.est        <- reshape2::melt(int.est)
names(int.est) <- c("Sample", "Population", "Value")
int.est %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = c(
    "Ambystoma cingulatum", "Anaxyrus boreas", "Pseudacris maculata"
  , "Notophthalmus viridescens", "Rana spp."
  ))
)

bd.est        <- reshape2::melt(bd.est)
names(bd.est) <- c("Sample", "Population", "Value")
bd.est %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = c(
    "Ambystoma cingulatum", "Anaxyrus boreas", "Pseudacris maculata"
  , "Notophthalmus viridescens", "Rana spp."
  ))
)

bd.mehg        <- reshape2::melt(bd.mehg)
names(bd.mehg) <- c("Sample", "Population", "Value")
bd.mehg        %<>% mutate(Species = as.factor(spec_in_pop[Population])) %>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = c(
    "Anaxyrus boreas", "Rana spp."
  ))
)

int.est %<>% mutate(Value = plogis(Value)) %>% group_by(Population, Species) %>% summarize(
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

ggplot(bd.est, aes(recapt_ind, CI_width)) + 
    geom_point(size = 2, aes(colour = reswabbed_ind)) + 
    scale_x_log10() +
    xlab("Number of Recaptured Individuals") +
    ylab("95% CI Width") +
  scale_color_continuous(name = "Total
Individuals
Reswabbed") +
theme(legend.key.size = unit(0.65, "cm"), legend.position = c(0.85, 0.7))

int.est %>% mutate(Population = plyr::mapvalues(Population, from = unique(Population), to = c(
"Ambystoma cingulatum
SMNWR East"
, "Ambystoma cingulatum
SMNWR West"
, "Anaxyrus boreas
Blackrock Complex"
, "Anaxyrus boreas
Blackrock H"
, "Anaxyrus boreas
Jones Pond"
, "Anaxyrus boreas
Sonoma Mountain"
, "Anaxyrus boreas
Two Medicine"
, "Pseudacris maculata
Lily Pond"
, "Pseudacris maculata
Matthews Pond"
, "Notophthalmus viridescens
Mud Lake"
, "Notophthalmus viridescens
SMNWR West"
, "Rana pretiosa
Dilman Meadows"
, "Rana boylii
Fox Creek"
, "Rana luteiventris
Jones Pond"
, "Rana luteiventris
Lost Horse"
, "Rana draytonii
San Francisquito"
, "Rana sierrae
Summit Meadow"
, "Rana cascadae
Three Creeks"
    ))
) %>% arrange(desc(mid)) %>% mutate(Population = factor(Population, levels = Population)) %>% {
  ggplot(., aes(mid, Population, colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "Dark2") +
    xlab("Survival at mean Bd load") +
    ylab("Population") +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13))
}

bd.est %>% mutate(Population = plyr::mapvalues(Population, from = unique(Population), to = c(
"Ambystoma cingulatum
SMNWR East"
, "Ambystoma cingulatum
SMNWR West"
, "Anaxyrus boreas
Blackrock Complex"
, "Anaxyrus boreas
Blackrock H"
, "Anaxyrus boreas
Jones Pond"
, "Anaxyrus boreas
Sonoma Mountain"
, "Anaxyrus boreas
Two Medicine"
, "Pseudacris maculata
Lily Pond"
, "Pseudacris maculata
Matthews Pond"
, "Notophthalmus viridescens
Mud Lake"
, "Notophthalmus viridescens
SMNWR West"
, "Rana pretiosa
Dilman Meadows"
, "Rana boylii
Fox Creek"
, "Rana luteiventris
Jones Pond"
, "Rana luteiventris
Lost Horse"
, "Rana draytonii
San Francisquito"
, "Rana sierrae
Summit Meadow"
, "Rana cascadae
Three Creeks"
    ))
) %>% arrange(desc(mid)) %>% mutate(Population = factor(Population, levels = Population)) %>% {
  ggplot(., aes(mid, Population, colour = Species)) + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "Dark2") +
    xlab("Bd-Survival Effect") +
    ylab("Population") +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13))
}

bd.mehg %>% mutate(Population = plyr::mapvalues(Population, from = unique(Population), to = c(
  "Anaxyrus boreas
Blackrock H"
, "Anaxyrus boreas
Jones Pond"
, "Anaxyrus boreas
Sonoma Mountain"
, "Rana pretiosa
Dilman Meadows"
, "Rana boylii
Fox Creek"
, "Rana luteiventris
Jones Pond"
, "Rana luteiventris
Lost Horse"
, "Rana draytonii
San Francisquito"
, "Rana cascadae
Three Creeks"
    ))
) %>% arrange(desc(mid)) %>% mutate(Population = factor(Population, levels = Population)) %>% {
  ggplot(., aes(mid, Population, colour = Species)) + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_point(size = 2) +
    scale_color_brewer(palette = "Dark2") +
    xlab("Bd-Survival Effect") +
    ylab("Population") +
    theme(axis.text.y = element_text(size = 11)
      , legend.key.size = unit(0.6, "cm")
      , legend.text = element_text(size = 11)
      , legend.title = element_text(size = 13))
}

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
  param_lev = plyr::mapvalues(param_lev, from = unique(beta_est.spec$param_lev), to = these_specs)
)
beta_est.slopes <- beta_est %>% filter(
  params %in% c("beta_bd_temp", "beta_bd_len", "beta_p_slope", "beta_mehg_drawdown")
)
} else {
beta_est.int <- beta_est %>% filter(
  params %in% c("beta_offseason_int", "beta_inseason_int", "beta_p_int", "beta_len", "beta_mehg_int")
) 
beta_est.slopes <- beta_est %>% filter(
  params %in% c("beta_bd_temp", "beta_bd_len", "beta_p_slope", "beta_mehg_drawdown")
)
}

gg3 <- beta_est.int %>% {
    ggplot(., aes(mid, param_lev)) + geom_point() +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3) +
      facet_wrap(~params, scales = "free") +
      ylab("Species")
}

if (n_specs > 1) {
gg4 <- beta_est.spec %>% {
    ggplot(., aes(mid, param_lev)) + geom_point() +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3) +
      facet_wrap(~params, scales = "free") +
      ylab("Species")
}
}

gg5 <- beta_est.slopes %>% {
    ggplot(., aes(mid, param_lev)) + geom_point() +
      geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.3) +
      facet_wrap(~params, scales = "free") +
      ylab("Species")
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

pop_size_ests <- capt_history.p %>% ungroup() %>% group_by(date_fac) %>% slice(1) %>% 
  dplyr::select(pop_spec, capture_date) %>% ungroup() %>%
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

x_labs <- pop_size_ests %>% group_by(pop_spec) %>% filter(
  date_fac %in% seq(min(date_fac), max(date_fac), by = 5)
)

gg10 <- pop_size_ests %>% {
    ggplot(., aes(date_fac, mid)) + 
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
      geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.2) +
      geom_line() +
      geom_point(aes(date_fac, capt_per_day), colour = "firebrick3", size = 3) +
      xlab("Date") +
      ylab("Population Estimate") +
      scale_x_continuous(
        breaks = x_labs$date_fac
      , labels = x_labs$capture_date
        ) +
      theme(axis.text.x = element_text(angle = 300, hjust = 0, size = 10)) +
      facet_wrap(~pop_spec, scales = "free") +
      scale_y_log10() +
      ggtitle("Red Points Show Number of Captures - Lines and Ribbons Show Population Estimates")
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
