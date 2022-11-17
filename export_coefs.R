##### Load and clean fits -----

all_pops <- c(
  "AMCI.SMNWR_E"
, "AMCI.SMNWR_W"
, 'ANBO.Blackrock'  
, "ANBO.JonesPond"
, "ANBO.SonomaMountain"
, "ANBO.TwoMedicine"
, "NOVI.KettleMoraine"
, "NOVI.MudLake"
, "NOVI.ScotiaBarrens"
, "NOVI.SMNWR_W"
, "NOVI.Springfield"
, "PSMA.LilyPond"
, "PSMA.MatthewsPond"
, "RANA.DilmanMeadows"
, "RANA.FoxCreek"
, "RANA.JonesPond"
, "RANA.LittleThreeCreeks"
, "RANA.LostHorse"
, "RANA.SanFrancisquito"
, "RANA.SummitMeadow"
)

file_samps <- list.files("fits/exported_samples")

## list of the pops in each fit
pops_in_fit <- list(
  c(4, 5, 12, 13, 14, 15, 16, 17, 18, 19)
, 14:20
, 3:6
, 1
, 2
, 11
, 4
, 10
, 9
, 7
, 8
)

## list of the specs of each pop in the fit
spec_in_fit <- list(
  c(1, 1, 2, 2, 3, 3, 3, 3, 3, 3)
, 1
, 1
, 1
, 1
, 1
, 1
, 1
, 1
, 1
, 1
)

for (i in seq_along(file_samps)) {

stan.fit.samples <- readRDS(
  paste("fits/exported_samples/", file_samps[i], sep = "")
  )

## quick change for inseason name
if (i == 1) {
  names(stan.fit.samples)[9] <- "beta_inseason_int"
}

print(paste("Starting", file_samps[i], sep = " ---- "))

samps.ss <- ifelse(length(grep("single", file_samps[i])) > 0, 1, 0)

if (samps.ss == 0 & i %notin% c(1)) {
  
stan.fit.samples$beta_offseason_int     <- stan.fit.samples$beta_phi[, 1:3]
stan.fit.samples$beta_offseason_bd      <- stan.fit.samples$beta_phi[, 4]
stan.fit.samples$beta_offseason_len     <- stan.fit.samples$beta_phi[, 5]
stan.fit.samples$beta_offseason_mehg    <- stan.fit.samples$beta_phi[, 6]

if (length(grep("mehg", file_samps[i])) == 1) {

stan.fit.samples$beta_offseason_mehg_bd    <- stan.fit.samples$beta_phi[, 7]
    
}
  
}
    
temp_array.int <- matrix(
  data = 0
, nrow = ifelse(samps.ss == 1, dim(stan.fit.samples$beta_offseason)[1], dim(stan.fit.samples$pop_size)[1])
, ncol = ifelse(samps.ss == 1, 1, dim(stan.fit.samples$z_r)[3])
, dimnames = list(
  NULL
, all_pops[pops_in_fit[[i]]]
)
)

temp_array.len <- matrix(
  data = 0
, nrow = ifelse(samps.ss == 1, dim(stan.fit.samples$beta_offseason)[1], dim(stan.fit.samples$pop_size)[1])
, ncol = ifelse(samps.ss == 1, 1, dim(stan.fit.samples$z_r)[3])
, dimnames = list(
  NULL
, all_pops[pops_in_fit[[i]]]
)
)

temp_array.bd <- matrix(
  data = 0
, nrow = ifelse(samps.ss == 1, dim(stan.fit.samples$beta_offseason)[1], dim(stan.fit.samples$pop_size)[1])
, ncol = ifelse(samps.ss == 1, 1, dim(stan.fit.samples$z_r)[3])
, dimnames = list(
  NULL
, all_pops[pops_in_fit[[i]]]
)
)

if (samps.ss == 0 | length(grep("mehg", file_samps[i])) == 1) {
  
temp_array.mehg <- matrix(
  data = 0
, nrow = ifelse(samps.ss == 1, dim(stan.fit.samples$beta_offseason)[1], dim(stan.fit.samples$pop_size)[1])
, ncol = ifelse(samps.ss == 1, 1, dim(stan.fit.samples$z_r)[3])
, dimnames = list(
  NULL
, all_pops[pops_in_fit[[i]]]
)
)
  
}

if (length(grep("mehg", file_samps[i])) == 1) {
  
temp_array.mehg_int <- matrix(
  data = 0
, nrow = ifelse(samps.ss == 1, dim(stan.fit.samples$beta_offseason)[1], dim(stan.fit.samples$pop_size)[1])
, ncol = ifelse(samps.ss == 1, 1, dim(stan.fit.samples$z_r)[3])
, dimnames = list(
  NULL
, all_pops[pops_in_fit[[i]]]
)
)

}

temp_array.within <- matrix(
  data = 0
, nrow = ifelse(samps.ss == 1, dim(stan.fit.samples$beta_offseason)[1], dim(stan.fit.samples$pop_size)[1])
, ncol = ifelse(samps.ss == 1, 1, dim(stan.fit.samples$z_r)[3])
, dimnames = list(
  NULL
, all_pops[pops_in_fit[[i]]]
)
)

if (samps.ss == 0) {

for (j in 1:length(pops_in_fit[[i]])) {
  
## the mehg model is a problem... Need to extract the correct intercepts and which Bd to
 ## use given that it varies by species
if (i == 1) {
  ## species than sex for the model matrix
  which_spec_this_ind <- spec_in_fit[[i]][j]
  
  if (which_spec_this_ind == 1) {
  beta_offseason_int <- c(stan.fit.samples$beta_offseason_int[, which_spec_this_ind])
  beta_offseason_bd  <- c(stan.fit.samples$beta_offseason_bd[, which_spec_this_ind])
  beta_offseason_len <- c(stan.fit.samples$beta_offseason_len[, which_spec_this_ind])
  beta_offseason_mehg    <- c(stan.fit.samples$beta_offseason_mehg[, which_spec_this_ind])
  beta_offseason_mehg_bd <- c(stan.fit.samples$beta_offseason_mehg_bd[, which_spec_this_ind])
  } else {
  beta_offseason_int <- c(stan.fit.samples$beta_offseason_int[, 1]) +
      c(stan.fit.samples$beta_offseason_int[, spec_in_fit[[i]][j]])
  beta_offseason_bd  <- c(stan.fit.samples$beta_offseason_bd[, 1]) + 
      c(stan.fit.samples$beta_offseason_bd[, spec_in_fit[[i]][j]])
  beta_offseason_len <- c(stan.fit.samples$beta_offseason_len[, 1]) + 
    c(stan.fit.samples$beta_offseason_len[, spec_in_fit[[i]][j]])
  beta_offseason_mehg    <- c(stan.fit.samples$beta_offseason_mehg[, 1]) + 
    c(stan.fit.samples$beta_offseason_mehg[, spec_in_fit[[i]][j]])
  beta_offseason_mehg_bd <- c(stan.fit.samples$beta_offseason_mehg_bd[, 1]) +
    c(stan.fit.samples$beta_offseason_mehg_bd[, spec_in_fit[[i]][j]])
  }
  
} else {
  beta_offseason_int  <- stan.fit.samples$beta_offseason_int[, 1]
  beta_offseason_bd   <- stan.fit.samples$beta_offseason_bd
  beta_offseason_len  <- stan.fit.samples$beta_offseason_len
  beta_offseason_mehg <- stan.fit.samples$beta_offseason_mehg
  
if (length(grep("mehg", file_samps[i])) == 1) {
  beta_offseason_mehg_bd <- stan.fit.samples$beta_offseason_mehg_bd
  }
}
  
  temp_array.int[, j] <- (beta_offseason_int  + stan.fit.samples$z_r[, 1, j]) %>% plogis()
  temp_array.bd[, j]  <- beta_offseason_bd    + stan.fit.samples$z_r[, 2, j]
  temp_array.len[, j] <- beta_offseason_len   + stan.fit.samples$z_r[, 3, j]
  if (i == 1) {
   temp_array.mehg[, j] <- beta_offseason_mehg + stan.fit.samples$z_r[, 4, j]
  } else {
   temp_array.mehg[, j] <- beta_offseason_mehg
  }
  
if (length(grep("mehg", file_samps[i])) == 1) {
  temp_array.mehg_int[, j]   <- beta_offseason_mehg_bd + stan.fit.samples$z_r[, 5, j]
}
   
  if (i == 1) {
    which_spec_this_ind <- spec_in_fit[[i]][j]
    if (which_spec_this_ind == 1) {
  temp_array.within[, j] <- (stan.fit.samples$beta_inseason_int[, 1] + 
    stan.fit.samples$inseason_pop[, j]) %>% plogis()
    } else {
  temp_array.within[, j] <- (stan.fit.samples$beta_inseason_int[, 1] + 
    stan.fit.samples$beta_inseason_int[, which_spec_this_ind] +
    stan.fit.samples$inseason_pop[, j]) %>% plogis()
    }
  } else {
  temp_array.within[, j] <- (stan.fit.samples$beta_inseason_int + 
    stan.fit.samples$inseason_pop[, j]) %>% plogis()
  }
  
}
  
temp_array.int <- reshape2::melt(temp_array.int) %>% rename(samp = Var1, pop = Var2) %>% 
  mutate(model = file_samps[i], parameter = "intercept")   
temp_array.bd  <- reshape2::melt(temp_array.bd) %>% rename(samp = Var1, pop = Var2) %>% 
  mutate(model = file_samps[i], parameter = "bd")
temp_array.len <- reshape2::melt(temp_array.len) %>% rename(samp = Var1, pop = Var2) %>% 
  mutate(model = file_samps[i], parameter = "length")

if (samps.ss == 0 | length(grep("mehg", file_samps[i])) == 1) {

temp_array.mehg  <- reshape2::melt(temp_array.mehg) %>% rename(samp = Var1, pop = Var2) %>% 
  mutate(model = file_samps[i], parameter = "mehg")  

temp_array.main <- rbind(temp_array.int, temp_array.bd, temp_array.len, temp_array.mehg)
  
} else {
  
temp_array.main <- rbind(temp_array.int, temp_array.bd, temp_array.len)
  
}

if (length(grep("mehg", file_samps[i])) == 1) {
  
temp_array.mehg_int <- reshape2::melt(temp_array.mehg_int) %>% rename(samp = Var1, pop = Var2) %>% 
  mutate(model = file_samps[i], parameter = "mehg_bd")

temp_array.main <- rbind(temp_array.main, temp_array.mehg_int)
  
}

temp_array.within <- reshape2::melt(temp_array.within) %>% 
  rename(samp = Var1, pop = Var2) %>% 
  mutate(model = file_samps[i], parameter = "within")
  
} else {
  
  temp_array.int[, 1] <- stan.fit.samples$beta_offseason_sex[, 1] %>% plogis()
  temp_array.bd[, 1]  <- stan.fit.samples$beta_offseason[, 1]
  temp_array.len[, 1] <- stan.fit.samples$beta_offseason[, 2]
  
if (length(grep("mehg", file_samps[i])) == 1) {
  
  temp_array.mehg[, 1]     <- stan.fit.samples$beta_offseason[, 3]
  temp_array.mehg_int[, 1] <- stan.fit.samples$beta_offseason[, 4]
  
}
  
  temp_array.within[, 1] <- stan.fit.samples$beta_phi %>% plogis()
  
temp_array.int <- reshape2::melt(temp_array.int) %>% rename(samp = Var1, pop = Var2) %>% 
  mutate(model = file_samps[i], parameter = "intercept")   
temp_array.bd  <- reshape2::melt(temp_array.bd) %>% rename(samp = Var1, pop = Var2) %>% 
  mutate(model = file_samps[i], parameter = "bd")
temp_array.len <- reshape2::melt(temp_array.len) %>% rename(samp = Var1, pop = Var2) %>% 
  mutate(model = file_samps[i], parameter = "length")

temp_array.main <- rbind(temp_array.int, temp_array.bd, temp_array.len)  

if (length(grep("mehg", file_samps[i])) == 1) {
  
temp_array.mehg     <- reshape2::melt(temp_array.mehg) %>% rename(samp = Var1, pop = Var2) %>% 
  mutate(model = file_samps[i], parameter = "mehg")
temp_array.mehg_int <- reshape2::melt(temp_array.mehg_int) %>% rename(samp = Var1, pop = Var2) %>% 
  mutate(model = file_samps[i], parameter = "mehg_bd")
  
temp_array.main <- rbind(temp_array.main, temp_array.mehg, temp_array.mehg_int) 

}

temp_array.within <- reshape2::melt(temp_array.within) %>% 
  rename(samp = Var1, pop = Var2) %>% 
  mutate(model = file_samps[i], parameter = "within")
  
}

if (i == 1) {
  f_array.main   <- temp_array.main
  f_array.within <- temp_array.within
} else {
  f_array.main   <- rbind(f_array.main, temp_array.main)
  f_array.within <- rbind(f_array.within, temp_array.within)
}

print(paste("Finished", file_samps[i], sep = " ---- "))
  
}

##### Calculate Bd effect on a probability scale -----

f_array.main.gg <- f_array.main %>% pivot_wider(., c(pop, model, samp)
  , values_from = "value", names_from = "parameter") %>% 
  mutate(
    avg_bd = plogis(intercept)
  , sd_bd  = plogis(intercept + bd * 1)
    ) %>% mutate(
    prop_bd = sd_bd / avg_bd
    ) %>% dplyr::select(
      -c(intercept, bd, length, mehg, mehg_bd, avg_bd, sd_bd)
    ) %>% 
  group_by(pop, model) %>%
  rename(value = prop_bd) %>%
  summarize(
    lwr   = quantile(value, 0.025)
  , lwr_n = quantile(value, 0.200)
  , mid   = quantile(value, 0.500)
  , upr_n = quantile(value, 0.800)
  , upr   = quantile(value, 0.975)
) 

##### Summarize the coefs -----

f_array.main <- f_array.within

f_array.main.gg <- f_array.main %>% group_by(pop, model, parameter) %>% summarize(
    lwr   = quantile(value, 0.025)
  , lwr_n = quantile(value, 0.200)
  , mid   = quantile(value, 0.500)
  , upr_n = quantile(value, 0.800)
  , upr   = quantile(value, 0.975)
) 

f_array.main.gg %<>% mutate(
  Species    = strsplit(pop %>% as.character(), "[.]")[[1]][1]
, Population = strsplit(pop %>% as.character(), "[.]")[[1]][2]
) 

f_array.main.gg %<>% ungroup() %>% mutate(
  Population = plyr::mapvalues(Population, from = unique(Population), to = c(
  "MT - Jones Pond"
, "CA - Sonoma Mountain"
, "CO - Lily Pond"
, "CO - Matthews Pond"
, "OR - Dilman Meadows"
, "CA - Fox Creek"
, "OR - Little Three Creeks"
, "MT - Lost Horse"
, "CA - San Francisquito"
, "CA - Summit Meadow"
, "WY - Blackrock"
, "MT - Two Medicine"
, "FL - SMNWR East"
, "FL - SMNWR West"
, "MA - Springfield"
, "PA - Scotia Barrens"
, "WI - Kettle Moraine"
, "WI - Mud Lake"
  ))
)

f_array.main.gg %<>% mutate(
  model = plyr::mapvalues(model, from = unique(model), to = c(
    "All MeHg pops"
  , "All Western Toad Populations"
  , "Single"
  , "All Rana Species"
  , "Single"
  , "Single"
  , "Single"
  , "Single"
  , "Single"
  , "Single"
  , "Single"
  ))
)

f_array.main.gg %<>% mutate(
  Species = plyr::mapvalues(Species, from = unique(Species), to = c(
    "Western toad"
  , "Boreal chorus frog"
  , "Rana sp."
  , "Frosted flatwoods salamander"
  , "Eastern newt"
)
  )
)

f_array.main.gg %<>% mutate(Pop_Mod_Spec = interaction(Population, model, Species)) %>%
  filter(Pop_Mod_Spec %notin% c(
#    "MT - Jones Pond.All Western Toad Populations.Western toad"
    "CA - Sonoma Mountain.All Western Toad Populations.Western toad"
  , "OR - Dilman Meadows.All Rana Species.Rana sp."
  , "CA - Fox Creek.All Rana Species.Rana sp."
  , "MT - Jones Pond.All Rana Species.Rana sp."
  , "OR - Little Three Creeks.All Rana Species.Rana sp."
  , "MT - Lost Horse.All Rana Species.Rana sp."
  , "CA - San Francisquito.All Rana Species.Rana sp."
  )) %>% mutate(
    Pop_Spec = interaction(Species, Population, sep = " : ")
  ) %>% mutate(
    Pop_Spec = plyr::mapvalues(Pop_Spec, from = c(
      "Rana sp. : OR - Dilman Meadows"
    , "Rana sp. : MT - Jones Pond"
    , "Rana sp. : MT - Lost Horse"
    , "Rana sp. : CA - Summit Meadow"
    , "Rana sp. : CA - Fox Creek"
    , "Rana sp. : OR - Little Three Creeks"
    , "Rana sp. : CA - San Francisquito"
    )
      , to = c(
      "Oregon spotted frog : OR - Dilman Meadows"
    , "Columbia spotted frog : MT - Jones Pond"
    , "Columbia spotted frog : MT - Lost Horse"
    , "Sierra Nevada yellow-legged frog : CA - Summit Meadow"
    , "Foothill yellow-legged frog : CA - Fox Creek"
    , "Cascades frog : OR - Little Three Creeks"
    , "California red-legged frog : CA - San Francisquito"
      ))
  )

## Just the three Jones Pond ANBO populations first to compare the estimates from fitting in three different ways
f_array.main.gg %>% filter(pop == "ANBO.JonesPond") %>% 
  mutate(parameter = factor(parameter, levels = c(
    "intercept", "length", "bd"
  , "mehg", "mehg_bd"
  ))) %>%
  mutate(parameter = plyr::mapvalues(parameter
  , from = c(
    "intercept", "length", "bd"
  , "mehg", "mehg_bd"
  )
  , to = c(
    "Intercept", "Length", "Bd"
  , "MeHg Main 
Effect", "MeHg-Bd
Interaction"
  )
  )) %>% arrange(desc(mid)) %>%
  mutate(Pop_Spec = factor(Pop_Spec, levels = unique(Pop_Spec))) %>% {
  ggplot(., aes(mid, model)) + 
    geom_point() + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
    ylab("Model Fit") +
    xlab("Estimate") +
    theme(axis.text.y = element_text(size = 10)) +
    facet_wrap(~parameter, nrow = 1)
  }

f_array.main.gg %<>% filter(Pop_Mod_Spec %notin% c(
    "MT - Jones Pond.All Western Toad Populations.Western toad"
,   "MT - Jones Pond.Single.Western toad"
))

##### Plots -----

f_array.main.gg %>% 
  filter(parameter == "intercept") %>% 
  arrange(desc(mid)) %>%
  mutate(Pop_Spec = factor(Pop_Spec, levels = unique(Pop_Spec))) %>% {
  ggplot(., aes(mid, Pop_Spec)) + 
    geom_point(aes(colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr, colour = Species), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n, colour = Species), height = 0.0, size = 1.5) +
  # geom_vline(xintercept = 1, linetype = "dashed", size = 0.5) +
 #   geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
    scale_colour_brewer(palette = "Dark2") +
    ylab("Population") +
    xlab("Estimate") +
  #  scale_x_continuous(breaks = c(0, 0.5, 0.75, 0.9, 1.0, 1.5, 2, 3)) +
    theme(
      axis.text.y = element_text(size = 10)
  #  , axis.text.x = element_text(size = 10, angle = 300, hjust = 0)
    ) 
}

f_array.main.gg %>% filter(parameter == "intercept") %>% arrange(desc(mid)) %>%
  mutate(Pop_Spec = factor(Pop_Spec, levels = unique(Pop_Spec))) %>% {
  ggplot(., aes(mid, Pop_Spec)) + 
    geom_point(aes(colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr, colour = Species), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n, colour = Species), height = 0.0, size = 1.5) +
    scale_colour_brewer(palette = "Dark2") +
    ylab("Population") +
    xlab("Estimate") +
    theme(
      axis.text.y = element_text(size = 10)
    ) 
  }

f_array.main.gg %>% 
  filter(parameter == "length") %>% arrange(desc(mid)) %>%
  mutate(Pop_Spec = factor(Pop_Spec, levels = unique(Pop_Spec))) %>% {
  ggplot(., aes(mid, Pop_Spec)) + 
    geom_point(aes(colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr, colour = Species), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n, colour = Species), height = 0.0, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
    scale_colour_brewer(palette = "Dark2") +
    ylab("Population") +
    xlab("Estimate") +
    theme(
      axis.text.y = element_text(size = 10)
    ) 
  }

f_array.main.gg %>% filter(parameter %in% c("mehg", "mehg_bd")) %>% 
  mutate(parameter = plyr::mapvalues(parameter, from = c("mehg", "mehg_bd"), to = c("MeHg Main
Effect", "MeHg Bd
Interaction"))) %>%
  mutate(parameter = factor(parameter, levels = c(
    "MeHg Main
Effect", "MeHg Bd
Interaction"
  ))) %>% arrange(Species) %>%
  mutate(Pop_Spec = factor(Pop_Spec, levels = unique(Pop_Spec))) %>% 
  {
  ggplot(., aes(mid, Pop_Spec)) + 
    geom_point(aes(colour = Species)) + 
    geom_errorbarh(aes(xmin = lwr, xmax = upr, colour = Species), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n, colour = Species), height = 0.0, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
    scale_colour_brewer(palette = "Dark2") +
    ylab("Population") +
    xlab("Estimate") +
    theme(axis.text.y = element_text(size = 10)) + 
      facet_wrap(~parameter, nrow = 1)
  }

##### Some figures for effort vs CI -----

### Also need to run manuscript_data_details2.R

# inds_capt
# prop_recap
# avg_caps
# avg_swabs
# prop_merc

gg.eff.1 <- f_array.main.gg %>% ungroup() %>% 
  left_join(., pops_info %>% ungroup() %>% rename(pop = pop_spec)) %>% 
  mutate(CI_width = upr - lwr) %>% 
  filter(parameter == "intercept") %>% 
  dplyr::select(Species, CI_width, inds_capt, prop_recap, avg_caps, avg_swabs) %>% {
   ggplot(., aes(avg_caps, CI_width)) + 
    geom_point(aes(colour = Species), size = 3) +
    scale_colour_brewer(palette = "Dark2") +
    xlab("Average Number of Captures per Individual") +
    ylab("Width of CI for Average Survival") +
      theme(
        legend.key.size = unit(.45, "cm")
      , plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
      )
  } 

gg.eff.2 <- f_array.main.gg %>% ungroup() %>% 
  left_join(., pops_info %>% ungroup() %>% rename(pop = pop_spec)) %>% 
  mutate(CI_width = upr - lwr) %>% 
  filter(parameter == "bd") %>% 
  dplyr::select(Species, CI_width, inds_capt, prop_recap, avg_caps, avg_swabs) %>% {
   ggplot(., aes(avg_swabs, CI_width)) + 
    geom_point(aes(colour = Species, size = avg_caps)) +
    scale_colour_brewer(palette = "Dark2", guide = "none") +
    xlab("Average Number of Swabs per Individual") +
    ylab("Width of CI for Bd Effect") +
      theme(
        legend.key.size = unit(.45, "cm")
      , plot.margin = unit(c(0.2,1.67,0.2,0.65), "cm")
      ) +
      scale_size_continuous(
  
       name   = "Average Number 
of Captures
per Individual"
      ) 
  }

gridExtra::grid.arrange(gg.eff.1, gg.eff.2, ncol = 1)


##### load and combine all predictions to plot survival over Bd load and length -----

file_samps <- list.files("fits/cleaned")

for (i in seq_along(file_samps)) {

pred_vals.gg <- readRDS(
 paste("fits/cleaned/", file_samps[i], sep = "")
)

if (class(pred_vals.gg)[1] == "list") {
  pred_vals.gg <- pred_vals.gg[[1]]
}

pred_vals.gg %<>% mutate(
  model = file_samps[i]
)

if (i %in% c(1,6,8)) {
  pred_vals.gg %<>% filter(sex == "M") %>% dplyr::select(-sex) 
  if (i == 6) {
  pred_vals.gg %<>% dplyr::select(-spec)
  }
}

if (i %notin% c(1,6,8)) {
  pred_vals.gg %<>% mutate(mehg = 0) %>% rename(pop = population) %>% 
    dplyr::select(-c(when, location, species)) %>% relocate(bd, len, mehg, pop, lwr, lwr_n, mid, upr_n, upr, model)
}

if (i == 1) {
  pred_vals.gg.f <- pred_vals.gg
} else {
  pred_vals.gg.f <- rbind(pred_vals.gg.f, pred_vals.gg)
}

}

ordered_pops <- unique(pred_vals.gg.f$pop) %>% as.character() %>% sort()
pred_vals.gg.f %<>% mutate(pop = factor(pop, levels = ordered_pops))
pred_vals.gg.f %<>% arrange(pop)

spp <- apply(pred_vals.gg.f$pop %>% matrix(), 1, FUN = function(x) strsplit(x, "[.]")[[1]][1])
pred_vals.gg.f %<>% mutate(spec = spp) %>%
  mutate(spec = plyr::mapvalues(spec, from = unique(spec), to = c(
    "Ambystoma cingulatum"
  , "Anaxyrus boreas"
  , "Notophthalmus viridescens"
  , "Pseudacris maculata"
  , "Rana spp"
  )))

pred_vals.gg.f %<>% mutate(
  pop_mod = interaction(pop, model)
)

pred_vals.gg.f %>% filter(
  pop_mod %notin% c(
    "ANBO.JonesPond.ANBO.Rds"
  , "ANBO.SonomaMountain.ANBO.Rds"
  , "RANA.FoxCreek.RANA.Rds"
  , "RANA.DilmanMeadows.RANA.Rds"
  , "RANA.JonesPond.RANA.Rds"
  , "RANA.LittleThreeCreeks.RANA.Rds"
  , "RANA.LostHorse.RANA.Rds"
  , "RANA.SanFrancisquito.RANA.Rds"
  )
) %>% mutate(
  pop = plyr::mapvalues(pop, from = unique(pop)
  , to = c(
"Ambystoma cingulatum
SMNWR East"

, "Ambystoma cingulatum
SMNWR West"

, "Anaxyrus boreas
Blackrock"

, "Anaxyrus boreas
Jones Pond"

, "Anaxyrus boreas
Sonoma Mountain"
  
, "Anaxyrus boreas
Two Medicine"
  
, "Notophthalmus viridescens
Kettle Moraine"   
    
, "Notophthalmus viridescens
Mud Lake"
  
, "Notophthalmus viridescens
Scotia Barrens"
  
, "Notophthalmus viridescens
SMNWR West"
  
, "Notophthalmus viridescens
Springfield"
    
, "Pseudacris maculata
Lily Pond"
  
, "Pseudacris maculata
Matthews Pond"
  
, "Rana pretiosa
Dilman Meadows"
  
, "Rana boylii
Fox Creek"
  
, "Rana luteiventris
Jones Pond"
    
, "Rana cascadae
Little Three Creeks" 
  
, "Rana luteiventris
Lost Horse"
  
, "Rana draytonii
San Francisquito"
  
, "Rana sierrae
Summit Meadow"
  
    )
)) %>% filter(len == 0, mehg == 0) %>% {
  ggplot(., aes(bd, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr
      , fill = spec, colour = spec
      ), alpha = 0.3) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n
      , fill = spec, colour = spec
      ), alpha = 0.3) +
    geom_line(aes(
      colour = spec
      ), size = 1) + 
    scale_colour_brewer(name = "Species", palette = "Dark2") +
    scale_fill_brewer(name = "Species", palette = "Dark2") +
    scale_x_continuous(breaks = c(-1.40, -0.75, 0, 0.75, 1.40)) +
    facet_wrap(~pop) +
    theme(
    strip.text.x    = element_text(size = 10)
  , axis.text.y     = element_text(size = 12)
  , axis.text.x     = element_text(size = 11)
  , legend.key.size = unit(0.75, "cm")
  , legend.spacing.y = unit(0.25, "cm")
  , legend.text     = element_text(size = 11)
  , legend.title    = element_text(size = 13)
  , legend.text.align = 0
    ) +
    xlab("Bd Load (scaled)") +
    ylab("Between-Season Survival")
}

pred_vals.gg.f %>% filter(
  pop_mod %notin% c(
    "ANBO.JonesPond.ANBO.Rds"
  , "ANBO.SonomaMountain.ANBO.Rds"
  , "RANA.FoxCreek.RANA.Rds"
  , "RANA.DilmanMeadows.RANA.Rds"
  , "RANA.JonesPond.RANA.Rds"
  , "RANA.LittleThreeCreeks.RANA.Rds"
  , "RANA.LostHorse.RANA.Rds"
  , "RANA.SanFrancisquito.RANA.Rds"
  )
) %>% mutate(
  pop = plyr::mapvalues(pop, from = unique(pop)
  , to = c(
"Ambystoma cingulatum
SMNWR East"

, "Ambystoma cingulatum
SMNWR West"

, "Anaxyrus boreas
Blackrock"

, "Anaxyrus boreas
Jones Pond"

, "Anaxyrus boreas
Sonoma Mountain"
  
, "Anaxyrus boreas
Two Medicine"
  
, "Notophthalmus viridescens
Kettle Moraine"   
    
, "Notophthalmus viridescens
Mud Lake"
  
, "Notophthalmus viridescens
Scotia Barrens"
  
, "Notophthalmus viridescens
SMNWR West"
  
, "Notophthalmus viridescens
Springfield"
    
, "Pseudacris maculata
Lily Pond"
  
, "Pseudacris maculata
Matthews Pond"
  
, "Rana pretiosa
Dilman Meadows"
  
, "Rana boylii
Fox Creek"
  
, "Rana luteiventris
Jones Pond"
    
, "Rana cascadae
Little Three Creeks" 
  
, "Rana luteiventris
Lost Horse"
  
, "Rana draytonii
San Francisquito"
  
, "Rana sierrae
Summit Meadow"
  
    )
)) %>% filter(bd == 0, mehg == 0) %>% {
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
    scale_colour_brewer(name = "Species", palette = "Dark2") +
    scale_fill_brewer(name = "Species", palette = "Dark2") +
    #scale_x_continuous(breaks = c(-1.40, -0.75, 0, 0.75, 1.40)) +
    xlim(c(-2, 2)) +
    facet_wrap(~pop) +
    theme(
    strip.text.x    = element_text(size = 10)
  , axis.text.y     = element_text(size = 12)
  , axis.text.x     = element_text(size = 11)
  , legend.key.size = unit(0.75, "cm")
  , legend.spacing.y = unit(0.25, "cm")
  , legend.text     = element_text(size = 11)
  , legend.title    = element_text(size = 13)
  , legend.text.align = 0
    ) +
    xlab("Individual length (SVL) (scaled)") +
    ylab("Between-Season Survival")
}


##### Detection -----
 ## (just a placeholder for now, will duplicate this for when the other fit eventually happens)

file_samps

stan.fit.samples <- readRDS(
  paste("fits/exported_samples/", file_samps[1], sep = "")
 )
hist(stan.fit.samples$beta_p_sex[, 2])

int.p           <- stan.fit.samples$beta_p_int 
dimnames(int.p) <- list(
  NULL
, c("ANBO", "PSMA", "RANA", "F", "U")
)
int.p[, "PSMA"] <- int.p[, "ANBO"] + int.p[, "PSMA"] 
int.p[, "RANA"] <- int.p[, "ANBO"] + int.p[, "RANA"] 
int.p[, "F"]    <- int.p[, "ANBO"] + int.p[, "F"] 

int.p           <- reshape2::melt(int.p)
names(int.p)    <- c("samp", "Species", "value")
int.p           %<>% mutate(value = plogis(value))

int.p           %<>% group_by(Species) %>% summarize(
    lwr   = quantile(value, 0.025)
  , lwr_n = quantile(value, 0.200)
  , mid   = quantile(value, 0.500)
  , upr_n = quantile(value, 0.800)
  , upr   = quantile(value, 0.975)
)
  
slope.p           <- stan.fit.samples$beta_p_slope
dimnames(slope.p) <- list(
  NULL
, c("Drawdown", "Vegetation")
)

slope.p           <- reshape2::melt(slope.p)
names(slope.p)    <- c("samp", "Covariate", "value")

slope.p           %<>% group_by(Covariate) %>% summarize(
    lwr   = quantile(value, 0.025)
  , lwr_n = quantile(value, 0.200)
  , mid   = quantile(value, 0.500)
  , upr_n = quantile(value, 0.800)
  , upr   = quantile(value, 0.975)
)


pop.p           <- stan.fit.samples$p_pop
dimnames(pop.p) <- list(
  NULL
, all_pops[pops_in_fit[[1]]]
)

pop.p           <- reshape2::melt(pop.p)
names(pop.p)    <- c("samp", "Population", "value")

pop.p           %<>% group_by(Population) %>% summarize(
    lwr   = quantile(value, 0.025)
  , lwr_n = quantile(value, 0.200)
  , mid   = quantile(value, 0.500)
  , upr_n = quantile(value, 0.800)
  , upr   = quantile(value, 0.975)
)

day.p <- stan.fit.samples$p_day_dev

day.p           <- reshape2::melt(day.p)
names(day.p)    <- c("samp", "dev", "value")

day.p           %<>% group_by(dev) %>% summarize(
    lwr   = quantile(value, 0.025, na.rm = T)
  , lwr_n = quantile(value, 0.200, na.rm = T)
  , mid   = quantile(value, 0.500, na.rm = T)
  , upr_n = quantile(value, 0.800, na.rm = T)
  , upr   = quantile(value, 0.975, na.rm = T)
)

day.p %<>% arrange(desc(mid)) %>% mutate(dev = factor(dev, levels = unique(dev)))

int.p <- rbind(int.p, slope.p %>% rename(Species = Covariate))

gg.p.int <- int.p %>% {
  ggplot(., aes(mid, Species)) + geom_point() +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
    ylab("Covariate") +
    xlab("Estimate") +
    scale_y_discrete(
      labels = rev(c(
        "Vegetation"
      , "Drawdown"
      , "Sex: Unknown"
      , "Sex: Female"
      , "Species: Rana"
      , "Species: PSMA"
      , "Species:ANBO
(intercept)")
    ))
}

gg.p.slope <- slope.p %>% {
  ggplot(., aes(mid, Covariate)) + geom_point() +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
    ylab("Covariate") +
    xlab("Estimate")
}

gg.p.pop <- pop.p %>% arrange(desc(mid)) %>% mutate(Population = factor(Population, levels = Population)) %>% {
  ggplot(., aes(mid, Population)) + geom_point() +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
    ylab("Population") +
    xlab("Estimate") 
}

gg.p.day <- day.p %>%slice(1:87) %>% {
  ggplot(., aes(mid, dev)) + geom_point() +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0.2, size = 0.8) +
    geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n), height = 0.0, size = 1.5) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
    ylab("PP-Subsite Deviates") +
    xlab("Estimate") +
    theme(
      axis.text.y = element_blank()
    )
}

gridExtra::grid.arrange(gg.p.pop, gg.p.day, nrow = 1)


##### Populations -----

for (i in seq_along(file_samps)) {

stan.fit.samples <- readRDS(
  paste("fits/exported_samples/", file_samps[i], sep = "")
 )

## Data frame of unique dates
uni_dates <- capt_history.p %>% group_by(pop_spec, capture_date) %>% slice(1) %>%
  dplyr::select(pop_spec, capture_date) %>% ungroup() %>%
  group_by(pop_spec) %>% mutate(
    entry = seq(n())
  )

pop_size <- stan.fit.samples$pop_size
#pop_size <- pop_size[, 1:((capt_history.p %>% group_by(pop_spec) %>% summarize(date_fac = length(unique(date_fac))) %>% 
#  filter(pop_spec %in% all_pops[pops_in_fit[[i]]]))$date_fac %>% sum())]

pop_entries <- rep(
  all_pops[pops_in_fit[[i]]]
, (capt_history.p %>% group_by(pop_spec) %>% 
   # summarize(date_fac = length(unique(date_fac))) %>% 
  summarize(capt_date = length(unique(capture_date))) %>% 
  filter(pop_spec %in% all_pops[pops_in_fit[[i]]]))$capt_date
   )

pop_entries <- (data.frame(
  pop = pop_entries
) %>% group_by(pop) %>%
  mutate(pop_entry = seq(n())) %>% 
  mutate(
    pop = interaction(pop, pop_entry)
  ))$pop

dimnames(pop_size) <- list(
  NULL
, pop_entries
)

pop_size <- reshape2::melt(pop_size) 
pop_size %<>% group_by(Var2) %>% summarize(
    lwr   = quantile(value, 0.025, na.rm = T)
  , lwr_n = quantile(value, 0.200, na.rm = T)
  , mid   = quantile(value, 0.500, na.rm = T)
  , upr_n = quantile(value, 0.800, na.rm = T)
  , upr   = quantile(value, 0.975, na.rm = T)
) %>% rename(pop_spec = Var2) 

pop_entries <- apply(pop_size$pop_spec %>% matrix, 1, FUN = function(x) strsplit(x, "[.]")[[1]][3]) %>% as.numeric()
pops        <- apply(pop_size$pop_spec %>% matrix, 1, FUN = function(x) strsplit(x, "[.]")[[1]][c(1,2)] %>% paste(collapse = "."))

pop_size %<>% mutate(
  pop_spec = pops
, entry    = pop_entries
  ) 

pop_size %<>% left_join(., uni_dates)

if (i == 1) {
  pop_size.f <- pop_size
} else {
  pop_size.f <- rbind(pop_size.f, pop_size)
}

}

pop_size %>% {
  ggplot(., aes(capture_date, mid)) + geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), height = 0.2, size = 0.8) +
    geom_errorbar(aes(ymin = lwr_n, ymax = upr_n), height = 0.0, size = 1.5) +
    scale_y_log10() +
    facet_wrap(~pop_spec, ncol = 1, scales = "free_y") +
    xlab("Date") + ylab("Population Size") 
}



