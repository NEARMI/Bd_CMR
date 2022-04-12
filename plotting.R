######################################################
## Plot diagnostics for each model fit individually ##
######################################################

# stan.fit <- readRDS(paste(paste("fits/stan_fit", which.dataset, sep = "_"), "Rds", sep = "."))
# stan.fit <- readRDS("fits/stan_fit_RANA.JonesPond_adj.Rds")

print("---------------------")
print("Model Finished and Saved, Extracting samples and starting plotting")
print("---------------------")

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

this_pop  <- capt_history$pop_spec[1] %>% as.character()
this_loc  <- capt_history$Site[1]     %>% as.character()
this_spec <- capt_history$Species[1]  %>% as.character()

####
## Plotting Setup
####

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

outval   <- seq(1, by = 1)
out.pred <- matrix(nrow = dim(stan.fit.samples[[1]])[1], ncol = length(outval))

for (j in 1:ncol(out.pred)) {
  out.pred[, j] <- plogis(stan.fit.samples$beta_phi)
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

# stan.fit.samples$beta_offseason <- cbind(stan.fit.samples$beta_offseason, rep(0, dim(stan.fit.samples$beta_offseason)[1]))

for (j in 1:ncol(out.pred)) {
   out.pred[, j] <- plogis(
    stan.fit.samples$beta_offseason[, 1] +
    stan.fit.samples$beta_offseason[, 2] * outval[j] +
    stan.fit.samples$beta_offseason[, 3] * 0 +
    stan.fit.samples$beta_offseason[, 4] * 0 +
    stan.fit.samples$beta_offseason[, 5] * 0 +
    stan.fit.samples$beta_offseason_sex[, 1]
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
    stan.fit.samples$beta_offseason[, 1] +
    stan.fit.samples$beta_offseason[, 2] * 7 +
    stan.fit.samples$beta_offseason[, 3] * 0 +
    stan.fit.samples$beta_offseason[, 4] * outval[j] +
    stan.fit.samples$beta_offseason[, 5] * outval[j] * 7 +
    stan.fit.samples$beta_offseason_sex[, 1]
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

capt_history.temp  <- capt_history %>% filter(pop_spec == this_pop)
capt_history.slice <- capt_history.temp %>% 
  group_by(capture_date) %>% 
  filter(captured == 1) %>% summarize(n_capts = n()) %>%
  arrange(n_capts) %>%
  mutate(rough_real_order = seq(n())) %>% 
  arrange(capture_date)

stan.p_pred_baseline <- stan.fit.samples$beta_p %>% reshape2::melt(.)

stan.p_pred_var <- stan.fit.samples$p_day_dev %>%
  reshape2::melt(.) %>% 
  rename(day = Var2, eps = value) %>% 
  left_join(., stan.p_pred_baseline) %>%
  mutate(pred_p = plogis(eps + value)) %>%
  group_by(day) %>%
  summarize(
    mid = quantile(pred_p, 0.50)
  , lwr = quantile(pred_p, 0.025)
  , upr = quantile(pred_p, 0.975)
  ) %>% arrange(day) %>%
  mutate(day = plyr::mapvalues(day, from = unique(day), to = as.character(unique(capt_history.temp$capture_date))))

stan.p_pred_var %<>% arrange(mid) %>% mutate(
  capture_date = as.Date(day)
, pred_order = seq(n())
  ) %>% dplyr::select(-day) %>% arrange(capture_date)

stan.p_pred_var %<>% left_join(., capt_history.slice)

stan.p_pred_var %<>% mutate(
    population = this_pop
  , location   = this_loc
  , species    = this_spec 
)

pop_size_est <- stan.fit.samples$pop_size %>% 
  reshape2::melt() %>% 
  mutate(value = ifelse(value == 0, NA, value)) %>%
  group_by(Var2) %>% 
  summarize(
    lwr   = quantile(value, 0.025, na.rm = T)
  , lwr_n = quantile(value, 0.200, na.rm = T)
  , mid   = quantile(value, 0.500, na.rm = T)
  , upr_n = quantile(value, 0.800, na.rm = T)
  , upr   = quantile(value, 0.975, na.rm = T)
  ) %>% 
  rename(Sample_Date = Var2) %>% 
  mutate(Sample_Date = as.character(Sample_Date)) %>%
  mutate(
    Sample_Date = plyr::mapvalues(
      Sample_Date
    , from = Sample_Date
    , to   = as.character(unique(capt_history$capture_date)))) %>%
      mutate(sdate = seq(1, n())) %>% filter(sdate > 1)

#ind_p_est <- stan.fit.samples$p_ind_dev %>% reshape2::melt() %>%
#  group_by(Var2) %>%
#  summarize(
#    lwr_p = quantile(value, 0.025)
#  , mid_p = quantile(value, 0.500)
#  , upr_p = quantile(value, 0.975)
#  ) %>% 
#  rename(ind = Var2) %>%
#  mutate(ind = factor(ind, levels = ind)) %>% 
#  arrange(mid_p) %>% mutate(
#   ord_p  = seq(n())
#  )

ind_bd_est <- stan.fit.samples$bd_delta_eps %>% reshape2::melt() %>%
  group_by(Var2) %>%
  summarize(
    lwr_bd = quantile(value, 0.025)
  , mid_bd = quantile(value, 0.500)
  , upr_bd = quantile(value, 0.975)
  ) %>% 
  rename(ind = Var2) %>%
  mutate(ind = factor(ind, levels = ind)) %>% 
  arrange(mid_bd) %>% mutate(
   ord_bd  = seq(n())
  )

## A bit of extra cleanup
param_names <- apply(
  matrix((beta_est %>% filter(population == unique(population)[1]))$Var1 %>% as.character())
, 1
, FUN = function(x) strsplit(x, "[[]")[[1]][1]
)

beta_est %<>% 
  group_by(population) %>%
  mutate(Var1 = as.character(Var1)) %>% 
  mutate(Var1 = plyr::mapvalues(Var1, from = unique(beta_est$Var1), to = param_names)) %>%
  rename(params = Var1) %>%
  group_by(params, population) %>%
  mutate(param_lev = seq(n())) %>% 
  relocate(param_lev, .after = params) %>%
  mutate(param_lev = as.character(param_lev))

beta_est[beta_est$params == "beta_offseason", ]$param_lev <- c("Int", "Bd", "Size", "MeHg", "Bd-MeHg")
beta_est %<>% mutate(params = plyr::mapvalues(params, from = "beta_phi", to = "beta_inseason"))
beta_est[beta_est$params == "beta_inseason", ]$param_lev <- c("Int")
beta_est[beta_est$params == "beta_p", ]$param_lev <- c("Int")

beta_est %<>% mutate(param = interaction(params, param_lev))


####
## Actual Plotting
####

gg.1 <- beta_est %>% 
    filter(params != "beta_bd_year", params != "beta_mehg_have", params != "ind_len_beta") %>% {
 #  filter(params == "beta_p") %>% {
    ggplot(., aes(param, mid)) +
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.3) +
      geom_point() +
      xlab("Parameter") +
      ylab("Estimate") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme(
        axis.text.x = element_text(angle = 300, hjust = 0)
      )
    }

gg.2a <- out.pred.bd %>% filter(when == "offseason", variable == "Bd") %>% {
  ggplot(., aes(gap, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line(size = 1) +
    xlab("") +
    ylab("Apparent Survival Between Seasons") +
    theme(
      axis.text.x = element_blank()
    , axis.ticks.x = element_blank()
    , plot.margin = unit(c(.2,.2,0,.2), "cm")
      )
}

gg.2b <- capt_history.bd_load %>% {
  ggplot(., aes(x = log_bd_load)) +
    geom_histogram(bins = 30) +
    xlab("Bd Copies (log)") +
    ylab("Density") +
    theme(
      plot.margin = unit(c(0,.2,.2,.33), "cm")
    )
}

gg.2 <- gridExtra::arrangeGrob(gg.2a, gg.2b, layout_matrix = rbind(c(1, 1), c(1, 1), c(2, 2)))

scaled_hg <- data.frame(ind_hg = scale(ind.hg)[, 1])

gg.3a <- out.pred.bd %>% filter(when == "offseason", variable == "MeHg") %>% {
  ggplot(., aes(gap, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line(size = 1) +
    xlab("") +
    ylab("Apparent Survival Between Seasons") +
    scale_x_continuous(lim = c(min(scaled_hg, na.rm = T), max(scaled_hg, na.rm = T))) +
    theme(
      axis.text.x = element_blank()
    , axis.ticks.x = element_blank()
    , plot.margin = unit(c(.2,.2,0,.2), "cm")
      )
}

gg.3b <- scaled_hg %>% {
  ggplot(., aes(x = ind_hg)) +
    geom_histogram(bins = 30) +
    xlab("MeHg (scaled)") +
    ylab("Density") +
    theme(
      plot.margin = unit(c(0,.2,.2,.66), "cm")
    )
}

gg.3 <- gridExtra::arrangeGrob(gg.3a, gg.3b, layout_matrix = rbind(c(1, 1), c(1, 1), c(2, 2)))

gg.4 <- stan.p_pred_var %>% 
    group_by(population) %>%
    filter(capture_date != min(capture_date)) %>%
    mutate(capture_date = as.factor(capture_date)) %>% {
  ggplot(., aes(mid, capture_date)) +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), size = 0.75, height = 0.3) +
    geom_point() +
    xlab("Daily Detection Probability (for each individual)") +
    ylab("Capture outing") 
}

gg.5 <- ind_p_est %>% mutate(ind = factor(ind, levels = ind)) %>%
    filter(ind %in% sample(seq(1, n_distinct(ind)), min(250, n_distinct(ind)))) %>% {
  ggplot(., aes(mid_p, ind)) + 
    geom_errorbarh(aes(xmin = lwr_p, xmax = upr_p), height = 0.3) + 
        geom_point() +
        theme(axis.text.y = element_text(size = 8)) +
        xlab("") +
        ylab("Individual") +
        xlab("Individual detection deviate")
    }

gg.6 <- ind_bd_est %>% mutate(ind = factor(ind, levels = ind)) %>%
    filter(ind %in% sample(seq(1, n_distinct(ind)), min(250, n_distinct(ind)))) %>% {
  ggplot(., aes(mid_bd, ind)) + 
    geom_errorbarh(aes(xmin = lwr_bd, xmax = upr_bd), height = 0.3) + 
        geom_point() +
        theme(axis.text.y = element_text(size = 8)) +
        xlab("") +
        ylab("Individual") +
        xlab("Individual bd deviate")
    }

gg.7 <- ind_order %>% {
  ggplot(., aes(order_real, order_pred)) +
    geom_point(aes(size = n_swabs)) +
    geom_abline(intercept = 0, slope = 1) +
    xlab("Real Bd Rank") +
    ylab("Estimated Bd Rank")
}

n_cap <- capt_history %>% group_by(capture_date) %>% summarize(num_capt = sum(captured)) %>%
  mutate(capture_date = as.factor(capture_date)) %>% mutate(capture_date = as.numeric(capture_date))
n_cap <- n_cap[-1, ]

gg.8 <- pop_size_est %>% {
    ggplot(., aes(sdate, mid)) + 
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
      geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.2) +
      geom_line() +
      scale_x_continuous(
        breaks = seq(nrow(pop_size_est))
      , labels = as.character(pop_size_est$Sample_Date)
        ) +
      geom_point(data = n_cap, aes(capture_date, num_capt), colour = "firebrick3", size = 3) +
      xlab("Date") +
      ylab("Population Estimate") +
      theme(axis.text.x = element_text(angle = 300, hjust = 0)) +
      ggtitle("Red Points Show Number of Captures - Lines and Ribbons Show Population Estimates")
}

if (max(pop_size_est$upr, na.rm = T) > (max(n_cap$num_capt) * 100)) {
  gg.8 <- gg.8 + scale_y_log10()
}

gglist    <- c("gg.1", "gg.2", "gg.3", "gg.4"
 # , "gg.5"
  , "gg.6", "gg.7", "gg.8")
need_grob <- c(F, T, T, F
#  , F
  , F, F, F)

pdf(paste("plots/", this_pop,".pdf", sep = ""), onefile = TRUE)
for (i in seq(length(gglist))) {
  if (need_grob[i]) {
   gridExtra::grid.arrange(get(gglist[i]))
  } else {
   get(gglist[i]) %>% print()
  }
}
dev.off()

