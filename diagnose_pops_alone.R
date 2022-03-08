#####################################################
## Diagnose and compare each single population fit ##
#####################################################

fits.files <- list.files("fits")

source("packages_functions.R")
red_ind    <- FALSE
single_pop <- FALSE
source("data_load.R")
source("data_manip.R")
source("data_stan.R")
source("data_covariates.R")

for (i in seq_along(fits.files)) {
  
this_file <- fits.files[i]
this_pop  <- strsplit(this_file, "_")[[1]]
this_pop  <- this_pop[length(this_pop)]
this_pop  <- strsplit(this_pop, "[.]")[[1]]
this_pop  <- this_pop[-length(this_pop)]
this_pop  <- paste(this_pop, collapse = ".")

this_loc  <- strsplit(this_pop, "[.]")[[1]][1]
this_spec <- strsplit(this_pop, "[.]")[[1]][2]
  
temp_fit  <- readRDS(paste("fits", this_file, sep = "/"))

stan.fit.summary <- summary(temp_fit)[[1]]
stan.fit.samples <- extract(temp_fit)

print("--------------")
print(paste("Samples Loaded for Population", i, sep = " "))
print("--------------")

beta_est <- stan.fit.summary[grep("beta", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>% 
  mutate(
    population = this_pop
  , location   = this_loc
  , species    = this_spec
    )

outval   <- seq(1, 10, by = 1)
out.pred <- matrix(nrow = dim(stan.fit.samples[[1]])[1], ncol = length(outval))

for (j in 1:ncol(out.pred)) {
  out.pred[, j] <- plogis(
    stan.fit.samples$beta_inseason_year[, 2] + 
    stan.fit.samples$beta_phi * outval[j] 
    )
}

out.pred <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "day_gap")

out.pred.in <- out.pred %>%
  group_by(gap, variable) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% mutate(
    when       = "inseason"
  , population = this_pop
  , location   = this_loc
  , species    = this_spec  
  )

outval   <- matrix(seq(1, 14, by = 1))
out.pred <- matrix(nrow = dim(stan.fit.samples[[1]])[1], ncol = length(outval))

for (j in 1:ncol(out.pred)) {
   out.pred[, j] <- plogis(
    stan.fit.samples$beta_offseason_year[, 2] + 
    stan.fit.samples$beta_offseason[, 1] * outval[j] +
    stan.fit.samples$beta_offseason[, 2] * 0 +
    stan.fit.samples$beta_offseason[, 3] * 0
    )
}

out.pred <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "Bd")

out.pred.off <- out.pred %>%
  group_by(gap, variable) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% mutate(
    when       = "offseason"
  , population = this_pop
  , location   = this_loc
  , species    = this_spec  
  )

outval   <- matrix(seq(-2, 2, by = .1))
out.pred <- matrix(nrow = dim(stan.fit.samples[[1]])[1], ncol = length(outval))

for (j in 1:ncol(out.pred)) {
   out.pred[, j] <- plogis(
    stan.fit.samples$beta_offseason_year[, 2] + 
    stan.fit.samples$beta_offseason[, 1] * 5 +
    stan.fit.samples$beta_offseason[, 2] * 0 +
    stan.fit.samples$beta_offseason[, 3] * outval[j]
    )
}

out.pred <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "MeHg")

out.pred.off.2 <- out.pred %>%
  group_by(gap, variable) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% mutate(
    when       = "offseason"
  , population = this_pop
  , location   = this_loc
  , species    = this_spec  
  )

out.pred.bd <- rbind(
  out.pred.in
, out.pred.off
, out.pred.off.2
)

outval   <- matrix(seq(1, 4, by = 1))
out.pred <- matrix(nrow = dim(stan.fit.samples[[1]])[1], ncol = length(outval))

for (j in 1:ncol(out.pred)) {
   out.pred[, j] <- plogis(stan.fit.samples$beta_p_year[, outval[j]])
}

out.pred <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "detection_year")

out.pred.p <- out.pred %>%
  group_by(gap, variable) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% mutate(
    population = this_pop
  , location   = this_loc
  , species    = this_spec  
  )

stan.ind_pred_var <- stan.fit.samples$X %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value) %>%
  group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% arrange(mid) %>%
  mutate(ind = factor(ind, levels = ind))

capt_history.temp <- capt_history %>% 
  filter(pop_spec == this_pop)
capt_history.slice <- capt_history.temp %>% group_by(X_stat_index) %>% slice(1)

stan.ind_pred_var <- cbind(
  Year = capt_history.slice$Year
, ind  = capt_history.slice$Mark
, Rep  = capt_history.slice$SecNumConsec
, t(stan.fit.samples$X)
) %>% as.data.frame() %>% 
  reshape2::melt(c("Year", "ind", "Rep")) %>% 
  dplyr::select(-Rep) %>%
  rename(iter = variable, eps = value) %>% 
  distinct() %>%
  group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% 
  arrange(mid) %>%
  mutate(ind = factor(ind, levels = ind))

stan.ind_pred_var %<>% mutate(
    population = this_pop
  , location   = this_loc
  , species    = this_spec  
)

ind_order.r <- capt_history.temp %>% 
  group_by(Mark) %>% 
  summarize(
   n_swabs = sum(swabbed)
) %>% left_join(
  ., capt_history.temp %>%
  filter(swabbed == 1) %>%
  group_by(Mark) %>% 
  summarize(
   tot_bd = mean(log_bd_load)
)) %>% mutate(
  given_mean = ifelse(
    is.na(tot_bd)
  , 1
  , 0
  ))

no_swabs <- (ind_order.r %>% filter(given_mean == 0))$Mark

ind_order.r %<>% 
  filter(given_mean == 0) %>%
  arrange(desc(tot_bd)) %>% 
  mutate(order_real = seq(n()), Mark = as.character(Mark))

ind_order.p <- stan.ind_pred_var %>% 
  filter(ind %in% no_swabs) %>%
  arrange(desc(mid)) %>% 
  mutate(order_pred = seq(n()), ind = as.character(ind)) %>% 
  rename(Mark = ind)

ind_order <- left_join(ind_order.r, ind_order.p)

ind_order %<>% mutate(
    population = this_pop
  , location   = this_loc
  , species    = this_spec 
)

if (i == 1) {
  beta_est.all          <- beta_est
  out.pred.all          <- out.pred.bd
  stan.ind_pred_var.all <- stan.ind_pred_var
  ind_order.all         <- ind_order
  out.pred.p.all        <- out.pred.p
} else {
  beta_est.all          <- rbind(beta_est.all, beta_est)
  out.pred.all          <- rbind(out.pred.all, out.pred.bd)
  stan.ind_pred_var.all <- rbind(stan.ind_pred_var.all, stan.ind_pred_var)
  ind_order.all         <- rbind(ind_order.all, ind_order)
  out.pred.p.all        <- rbind(out.pred.p.all, out.pred.p)
}

print("--------------")
print(paste("Through", i, "of", length(seq_along(fits.files)), "sites", sep = " "))
print("--------------")

}

## Some final cleanup

param_names <- apply(
  matrix((beta_est.all %>% filter(population == unique(population)[1]))$Var1 %>% as.character())
, 1
, FUN = function(x) strsplit(x, "[[]")[[1]][1]
)

beta_est.all %<>% 
  group_by(population) %>%
  mutate(param_names = param_names) %>%
  group_by(param_names, population) %>%
  mutate(param_lev = seq(n()))

####
## And plotting
####

beta_est.all %>% 
  filter((param_names == "beta_phi" & param_lev == 2) | (param_names == "beta_offseason" & param_lev == 1)) %>% {
    ggplot(., aes(population, mid)) +
      geom_errorbar(aes(ymin = lwr, ymax = upr, colour = species), width = 0.3) +
      geom_point(aes(colour = species)) +
      scale_color_brewer(palette = "Dark2", name = "Species") +
      xlab("Population") +
      ylab("Estimate") +
      facet_wrap(~param_names, scales = "free") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme(
        axis.text.x = element_text(angle = 300, hjust = 0)
      )
}

out.pred.all %>% filter(when == "offseason", variable == "Bd") %>% {
  ggplot(., aes(gap, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = location), alpha = 0.3) +
    geom_line(aes(colour = location), size = 1) +
    xlab("Bd Copies (log)") +
    ylab("Apparent Survival Between Seasons") +
    facet_wrap(~species)
}

out.pred.all %>% filter(when == "offseason", variable == "MeHg") %>% {
  ggplot(., aes(gap, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = location), alpha = 0.3) +
    geom_line(aes(colour = location), size = 1) +
    xlab("Bd Copies (log)") +
    ylab("Apparent Survival Between Seasons") +
    facet_wrap(~species)
}

out.pred.p.all %>% {
  ggplot(., aes(gap, mid)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, colour = location), width = 0.3, size = 1) +
    xlab("Size") +
    ylab("Detection Probability") +
    facet_wrap(~species)
}

ind_order.all %>% {
  ggplot(., aes(order_real, order_pred)) +
    geom_point(aes(colour = species, size = n_swabs)) +
    xlab("Real Bd Rank") +
    ylab("Estimated Bd Rank") +
    facet_wrap(~population, scales = "free")
}

capt_history.phi %>% 
  filter(pop_spec == "MatthewsPond.BCF") %>% 
  group_by(Mark) %>% 
  summarize(nswabs = sum(swabbed)) %>% 
  arrange(nswabs)
