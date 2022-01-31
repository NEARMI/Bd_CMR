#############################################################
## Extract the fit for each individual population from the ##
## joint model that fit all populations at once            ##
#############################################################

the_specs <- unique(capt_history.phi$Species)
the_sites <- unique(capt_history.phi$Site)
the_pops  <- unique(capt_history.phi$pop_spec)

for (i in seq_along(the_pops)) {

this_pop   <- capt_history %>% filter(pop_spec == the_pops[i])
which_spec <- which(the_specs == this_pop$Species[1])
which_site <- which(the_sites == this_pop$Site[1])
which_pop  <- i  
which_X    <- unique(this_pop$X_stat_index)


outval   <- matrix(seq(1, 10, by = 1))
out.pred <- matrix(nrow = dim(stan.fit.samples[[1]])[1], ncol = length(outval))

for (j in 1:ncol(out.pred)) {
   out.pred[, j] <- plogis(
    stan.fit.samples$beta_offseason[, 1] +
    stan.fit.samples$offseason_pop[, which_pop] +
    (stan.fit.samples$beta_offseason[, 2] + stan.fit.samples$beta_spec[, which_spec]) * outval[j] +
    stan.fit.samples$beta_offseason[, 3] * 0
    )
}

out.pred <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "bd")

out.pred.phi <- out.pred %>%
  group_by(gap, variable) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% mutate(
    when       = "offseason"
  , population = the_pops[which_pop]
  , location   = the_sites[which_site]
  , species    = the_specs[which_spec]
  )

outval   <- matrix(seq(1, 10, by = 1))
out.pred <- matrix(nrow = stan.length, ncol = length(outval))

for (j in 1:ncol(out.pred)) {
   out.pred[, j] <- plogis(
    stan.fit.samples$beta_p[, 1]  + stan.fit.samples$p_pop[, which_pop] +
    (stan.fit.samples$beta_p[, 2] + stan.fit.samples$p_pop_bd[, which_pop]) * outval[j] +
    stan.fit.samples$beta_p[, 3] * 0 # + stan.fit.samples$beta_p[, 4] * 0
    )
}

out.pred <- reshape2::melt(out.pred) %>% 
  rename(iter = Var1, gap = Var2) %>% 
  mutate(gap = plyr::mapvalues(gap
    , from = unique(gap), to = unique(outval))) %>%
  mutate(variable = "bd")

out.pred.p <- out.pred %>%
  group_by(gap, variable) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% mutate(
    when       = "offseason"
  , population = the_pops[which_pop]
  , location   = the_sites[which_site]
  , species    = the_specs[which_spec]
  )

stan.ind_pred_var <- stan.fit.samples$X[, which_X] %>%
  reshape2::melt(.) %>% 
  rename(ind = Var2, eps = value) %>%
  group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% arrange(mid) %>%
  mutate(ind = factor(ind, levels = ind))

capt_history.slice <- this_pop %>% group_by(X_stat_index) %>% slice(1)

stan.ind_pred_var <- cbind(
  Year = capt_history.slice$Year
, ind  = capt_history.slice$Mark
, Rep  = capt_history.slice$SecNumConsec
, t(stan.fit.samples$X[, which_X])
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
    population = the_pops[which_pop]
  , location   = the_sites[which_site]
  , species    = the_specs[which_spec]
)

ind_order.r <- this_pop %>% 
  group_by(Mark) %>% 
  summarize(
   n_swabs = sum(swabbed)
) %>% left_join(
  ., this_pop %>%
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
    population = the_pops[which_pop]
  , location   = the_sites[which_site]
  , species    = the_specs[which_spec]
)

if (i == 1) {
  out.pred.all.phi      <- out.pred.phi
  out.pred.all.p        <- out.pred.p
  stan.ind_pred_var.all <- stan.ind_pred_var
  ind_order.all         <- ind_order
} else {
  out.pred.all.phi      <- rbind(out.pred.all.phi, out.pred.phi)
  out.pred.all.p        <- rbind(out.pred.all.p, out.pred.p)
  stan.ind_pred_var.all <- rbind(stan.ind_pred_var.all, stan.ind_pred_var)
  ind_order.all         <- rbind(ind_order.all, ind_order)
}

print("--------------")
print(paste("Through", i, "of", n_sites, "sites", sep = " "))
print("--------------")

}

####
## And plotting
####

stan.fit.summary[grep("beta", dimnames(stan.fit.summary)[[1]]), ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>% {
    ggplot(., aes(Var1, mid)) + 
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3) +
      geom_point() +
      geom_hline(yintercept = 0) +
      xlab("Parameter") +
      ylab("Estimate") +
      theme(axis.text.x = element_text(angle = 300, hjust = 0))
  }

stan.fit.summary[grep("offseason_pop", dimnames(stan.fit.summary)[[1]]), ][15:27, ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>% {
    ggplot(., aes(Var1, mid)) + 
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3) +
      geom_point() +
      geom_hline(yintercept = 0) +
      xlab("Parameter") +
      ylab("Estimate") +
      theme(axis.text.x = element_text(angle = 300, hjust = 0)) +
      scale_x_discrete(
        labels = the_pops
      )
  }

stan.fit.summary[grep("p_pop", dimnames(stan.fit.summary)[[1]]), ][29:41, ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>% {
    ggplot(., aes(Var1, mid)) + 
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3) +
      geom_point() +
      geom_hline(yintercept = 0) +
      xlab("Parameter") +
      ylab("Estimate") +
      theme(axis.text.x = element_text(angle = 300, hjust = 0)) +
      scale_x_discrete(
        labels = the_pops
      )
  }

stan.fit.summary[grep("p_pop", dimnames(stan.fit.summary)[[1]]), ][42:54, ] %>% 
  reshape2::melt() %>%
  filter(Var2 %in% c('2.5%', '50%', '97.5%')) %>% 
  pivot_wider(names_from = "Var2", values_from = "value") %>% 
  rename(lwr = '2.5%', mid = '50%', upr = '97.5%') %>% {
    ggplot(., aes(Var1, mid)) + 
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3) +
      geom_point() +
      geom_hline(yintercept = 0) +
      xlab("Parameter") +
      ylab("Estimate") +
      theme(axis.text.x = element_text(angle = 300, hjust = 0)) +
      scale_x_discrete(
        labels = the_pops
      )
  }

out.pred.all.phi %>% filter(when == "offseason") %>% 
  mutate(species = factor(species, levels = c(
    "ANBO", "BCF", "CF", "FYLF", "RALU", "SNYLF"
    ))) %>% {
  ggplot(., aes(gap, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = location), alpha = 0.3) +
    geom_line(aes(colour = location), size = 1) +
    xlab("Bd Copies (log)") +
    ylab("Apparent Survival Between Seasons") +
    facet_wrap(~species)
    }

out.pred.all.p %>% filter(when == "offseason") %>% 
  mutate(species = factor(species, levels = c(
    "ANBO", "BCF", "CF", "FYLF", "RALU", "SNYLF"
    ))) %>% {
  ggplot(., aes(gap, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = location), alpha = 0.3) +
    geom_line(aes(colour = location), size = 1) +
    xlab("Bd Copies (log)") +
    ylab("Detection") +
    facet_wrap(~species)
}

ind_order.all %>% {
  ggplot(., aes(order_real, order_pred)) +
    geom_point(aes(colour = species, size = n_swabs)) +
    xlab("Real Bd Rank") +
    ylab("Estimated Bd Rank") +
    facet_wrap(~population, scales = "free")
}

####
## Some exploration about the odd populations
####

pop_plot <- "MatthewsPond.BCF"

pop_plot.d <- capt_history %>% 
  filter(pop_spec == pop_plot) %>% 
  group_by(Mark, Year) %>%
  filter(captured == 1) %>% 
  summarize(
    capt    = sum(captured)
  , mean_bd = mean(log_bd_load)
    ) %>%
  left_join(., capt_history) %>%
  filter(capt > 0) %>% 
  group_by(Mark, Year) %>%
  slice(1) %>% 
  ungroup(Year) %>%
  mutate(capt_number = seq(n())) %>% 
  mutate(recapture  = ifelse(capt_number > 1, 1, 0)) %>% 
  mutate(recaptured = ifelse(sum(recapture) > 0, 1, 0)) %>% 
  mutate(captured_again = c(recapture[-1], NA)) %>% 
  mutate(captured_again = ifelse(is.na(captured_again), 0, 1)) 

glm(
   captured_again ~ mean_bd
 , family = "binomial"
 , data = pop_plot.d
) %>% summary()

pop_plot.d %>% {
    ggplot(., aes(mean_bd, recaptured)) + stat_sum() 
  }

capt_history %>% 
  filter(pop_spec == "MatthewsPond.BCF") %>% 
  group_by(Mark) %>%
  mutate() %>% 
  summarize(total_caps = sum(captured)) %>%
  left_join(., capt_history) %>% 
  


capt_history %>% 
  group_by(Mark, pop_spec) %>% 
  summarize(total_caps = sum(captured)) %>%
  left_join(., capt_history) %>% {
    ggplot(., aes(log_bd_load, total_caps)) + 
      geom_point() +
      ylab("Total Times Captured") +
      xlab("Bd Load (log10)") +
      facet_wrap(~pop_spec, scales = "free")
  }
