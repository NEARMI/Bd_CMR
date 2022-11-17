####
## Script for temporary plots for first document for PIs
####

##### Single population fits -----------------------------------------------

# %% 0) sp_Bd_rep_meas
# %% 1) sp_Bd_modes
# %% 2) sp_Bd_rank
# %% 3) sp_neg_Bdsurv1
# %% 4) sp_neg_Bdsurv2
# %% 5) sp_pos_Bdsurv1
# %% 6) sp_pos_Bdsurv2
# %% 7) sp_detect
# %% 8) sp_pop_good
# %% 9) sp_pop_bad

# 0) measures of Bd through time

which_mark_sub <- (capt_history.bd_load %>% group_by(Mark) %>% 
  summarize(tot_caps = sum(swabbed)) %>% 
  arrange(desc(tot_caps)) %>% head(20))$Mark

capt_history.bd_load %>% filter(Mark %in% which_mark_sub) %>% {
  ggplot(., aes(capture_date, bd_load)) + 
    geom_point() +
    scale_y_log10() +
    facet_wrap(~Mark) +
    xlab("Capture Date") +
    ylab("Bd Load")
}

capt_history.bd_load %>%
  mutate(JD = as.POSIXlt(capture_date)$yday) %>%
  filter(Mark %in% which_mark_sub) %>% {
  ggplot(., aes(JD, bd_load)) +
    geom_point(aes(colour = as.factor(Year))) +
    geom_path(aes(colour = as.factor(Year))) +
    scale_y_log10() +
    facet_wrap(~Mark) +
    scale_colour_brewer(
      palette = "Dark2"
    , name = "Year"
      ) +
    xlab("Julian Day") +
    ylab("Bd Load") +
    theme(
      legend.key.size = unit(0.8, "cm")
    , strip.background = element_blank()
    , strip.text.x = element_blank()
    , axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
    , legend.text = element_text(size = 14)
    , legend.title = element_text(size = 16)
      )
}

# %% 1) sp_Bd_modes: 574 x 798; 5.98 x 8.49
ind_bd_est %>% 
  mutate(ind = factor(ind, levels = ind)) %>%{
  ggplot(., aes(mid_bd, ind)) + 
    geom_errorbarh(aes(xmin = lwr_bd, xmax = upr_bd), height = 0.2) + 
        geom_point() +
        theme(axis.text.y = element_text(size = 8)) +
        ylab("Individual") +
        xlab("Individual Bd Deviate") +
        theme(axis.text.y = element_blank())
    }

# %% 2) sp_Bd_rank: 6.03 x 8.49
ind_order %>% {
  ggplot(., aes(order_real, order_pred)) +
    geom_point(aes(size = n_swabs)) +
    geom_abline(intercept = 0, slope = 1) +
    scale_size_continuous(name = "Number
of Swabs") +
    xlab("Measured Bd Rank") +
    ylab("Estimated Bd Rank") 
}

# %% 3) sp_neg_Bdsurv: 574 x 798; 5.98 x 8.49

gg.2a <- pred.vals.gg %>% filter(len == 0
#  , mehg == 0
  ) %>% {
  ggplot(., aes(bd, mid)) + 
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.3) +
    geom_line(size = 1) + 
    xlab("Bd Load") + 
    ylab("Apparent Survival Between Seasons") +
    ggtitle(this_pop)
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

gridExtra::grid.arrange(gg.2)

# %% 4) sp_pos_Bdsurv ^^ Above but with a different population

# %% 7) sp_detect; 5.98 x 8.49

stan.p_pred_var %>% 
    group_by(population) %>%
    filter(capture_date != unique(capture_date)[c(1,2)]) %>%
    mutate(capture_date = as.factor(capture_date)) %>% {
  ggplot(., aes(mid, capture_date)) +
    geom_errorbarh(aes(xmin = lwr, xmax = upr), size = 0.75, height = 0.3) +
    geom_point() +
    scale_x_continuous(
   # breaks = c(0.0, 0.15, 0.30, 0.45)
     breaks = c(0.0, 0.025, 0.050, 0.075, 0.100)
      ) +
    xlab("Daily Detection Probability (for each individual)") +
    ylab("Capture outing") +
    geom_hline(aes(yintercept = 15.5), colour = "blue", linetype = "dashed", size = 0.25) +
    geom_hline(aes(yintercept = 39.5), colour = "blue", linetype = "dashed", size = 0.25) +
    theme(axis.text.y = element_text(size = 12))
}

# %% 8) sp_pop_good

x_axis_breaks <- seq(1, nrow(pop_size_est), length = 10) %>% round()

## 929 x 624, 9.68 x 6.60
pop_size_est %>% {
    ggplot(., aes(sdate, mid)) + 
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
      geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.2) +
      geom_line() +
      scale_x_continuous(
        breaks = x_axis_breaks
      , labels = as.character(pop_size_est$Sample_Date)[x_axis_breaks]
        ) +
  #     scale_y_continuous(
  #      breaks = c(50, 100, 150, 200, 250)
  #    , labels = c("   50", "  100", "  150", "  200", "  250")
  #     ) +
      scale_y_log10(
        breaks = c(1, 10, 100, 1000, 10000, 100000)
      ) +
      geom_point(data = n_cap, aes(capture_date, num_capt), colour = "firebrick3", size = 3) +
      geom_vline(xintercept = 16.5, linetype = "dashed", colour = "dodgerblue3", size = 0.5) +
      geom_vline(xintercept = 40.5, linetype = "dashed", colour = "dodgerblue3", size = 0.5) +
 #     geom_vline(xintercept = 17.5, linetype = "dashed", colour = "dodgerblue3", size = 0.5) +
      xlab("Date") +
      ylab("Population Estimate") +
      theme(axis.text.x = element_text(angle = 300, hjust = 0)) +
      ggtitle("Red Points Show Number of Captures - Lines and Ribbons Show Population Estimates")
}

# %% A few plots for the MA newts ^^ Just use the above +

capt_history %<>% cbind(., t(stan.fit.samples$X))
all_names  <- names(capt_history)[1:26]
names_need <- c("capture_date", "Mark", "swabbed", "log_bd_load", "X_stat_index")

capt_history %<>% dplyr::select(-c(all_names[all_names %notin% names_need]))
capt_history %<>% filter(swabbed == 1)

capt_history.2 <- capt_history %>% pivot_longer(-c(capture_date, Mark, swabbed, log_bd_load, X_stat_index))
capt_history.2 %<>% dplyr::select(-swabbed)
capt_history.2 %<>% group_by(Mark, X_stat_index, name) %>% summarize(max_val = max(value))

capt_history.2 %<>% ungroup() %>% group_by(Mark, X_stat_index) %>% summarize(
    lwr   = quantile(max_val, 0.025, na.rm = T)
  , lwr_n = quantile(max_val, 0.200, na.rm = T)
  , mid   = quantile(max_val, 0.500, na.rm = T)
  , upr_n = quantile(max_val, 0.800, na.rm = T)
  , upr   = quantile(max_val, 0.975, na.rm = T)
)

capt_history.r <- capt_history %>% dplyr::select(capture_date, Mark, swabbed, log_bd_load, X_stat_index)
capt_history.r %<>% droplevels()

JD  <- as.POSIXlt(capt_history.r$capture_date)$yday
yrr <- apply(matrix(as.character(capt_history.r$capture_date)), 1, FUN = function(x) strsplit(x, "-")[[1]][1] %>% as.numeric())

capt_history.r %<>% mutate(JD = JD, year = yrr)

capt_history.2 %<>% left_join(.
  , capt_history.r %>% dplyr::select(X_stat_index, year) %>% distinct()
  )

ind_order_for_plot <- (capt_history.r %>% group_by(Mark) %>% 
  summarize(mean_bd = mean(log_bd_load)) %>% arrange(desc(mean_bd)))$Mark

capt_history.r %<>% mutate(Mark = factor(Mark, levels = ind_order_for_plot)) 
capt_history.2 %<>% mutate(Mark = factor(Mark, levels = ind_order_for_plot))

##### Plot 1 being all individual data points and maxes
capt_history.r %>% mutate(Mark = factor(Mark, levels = unique(Mark))) %>% {
  ggplot(., aes(log_bd_load, Mark)) + 
    geom_point(aes(colour = as.factor(year)), alpha = 0.5) +
#    geom_point(data = capt_history.2
#      , aes(mid, Mark
#      , colour = as.factor(year))
#      , shape = 2) +
#     geom_errorbarh(data = capt_history.2 %>% rename(log_bd_load = mid)
#      , aes(xmin = lwr_n, xmax = upr_n, y = Mark
#      , colour = as.factor(year)), alpha = 0.3) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 15, 18)) +
    xlab("Bd load (Log)") +
    ylab("Individual") + 
    theme(
      axis.text.y = element_blank()
    , legend.position = c(0.9, 0.9))
}

##### Plot 2 being CI for a few specific individuals

inds_to_plot <- (capt_history.r %>% group_by(Mark) %>% summarize(num_swabs = n()) %>% arrange(desc(num_swabs)))$Mark[1:20]

capt_history.r %>% filter(Mark %in% inds_to_plot) %>% {
  ggplot(., aes(JD, log_bd_load)) + 
    geom_point(aes(colour = as.factor(year))) +
    geom_line(aes(colour = as.factor(year))) +
    geom_errorbar(data = 
        left_join(
          expand.grid(Mark = inds_to_plot %>% as.factor(), JD = c(100:200)) 
        , capt_history.2 %>% rename(log_bd_load = mid) %>% filter(Mark %in% inds_to_plot)
        )
     , aes(ymin = lwr, ymax = upr, x = JD
     , colour = as.factor(year)), alpha = 0.4) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    theme(axis.text.x = element_text(size = 10)) +
    xlab("Julian Day") +
    ylab("Bd load (Log)") +
    facet_wrap(~Mark)
}


##### Single Species Multiple population fits -----------------------------------------------

# %% 10) ssp_Bd_ANBO

pred.vals.gg %>% filter(sex == "M", len == 0, mehg == 0) %>% 
mutate(
  pop = plyr::mapvalues(pop, from = unique(pop), to = c(
    "Blackrock-C", "Blackrock-H", "Jones Pond", "Sonoma Mountain", "Two Medicine"
  ))) %>% {
  ggplot(., aes(bd, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.3) +
    geom_line(size = 1) + 
    scale_colour_discrete() +
    scale_fill_discrete() +
  #  scale_x_continuous(breaks = c(0, 3, 6, 9, 12)) +
    scale_x_continuous(breaks = c(-1.5, -0.75, 0, 0.75, 1.5)) +
    scale_y_continuous(lim = c(0, 1)) +
    facet_wrap(~pop, ncol = 1, strip.position = "right") +
    theme(
      strip.text = element_text(size = 12)
    ) +
    xlab("Bd Load (Scaled)") +
    ylab("Between-Season Survival")
  }

## Quick asside prior vs posterior

compare_prior <- data.frame(Posterior = stan.fit.samples$beta_offseason_int[, 1]
  , Prior = rnorm(3000, 0, 0.65), intt = seq(3000)) %>% pivot_longer(-intt)

compare_prior %>% {
  ggplot(., aes(x = value)) + 
    geom_density(aes(colour = name), size = 2) +
    scale_colour_brewer(name = "Distribution", palette = "Dark2") +
    theme(
      legend.key.size = unit(0.8, "cm")
    , legend.text = element_text(size = 12)
    , legend.title = element_text(size = 14)
    ) +
    xlab("Density") +
    ylab("Value")
}

# %% 11) ssp_Bd_RANA

pred.vals.gg %>% filter(sex == "M", len == 0, mehg == 0) %>% 
mutate(
  pop = plyr::mapvalues(pop, from = unique(pop), to = c(
    "Rana pretiosa
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
) %>% {
  ggplot(., aes(bd, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.3) +
    geom_line(size = 1) + 
    scale_colour_discrete() +
    scale_fill_discrete() +
    scale_x_continuous(breaks = c(0, 3, 6, 9, 12)) +
    facet_wrap(~pop, ncol = 1, strip.position = "right") +
    theme(
      strip.text = element_text(size = 11)
  ,   axis.text.y = element_text(size = 12)
    ) +
    xlab("Bd Load") +
    ylab("Between-Season Survival")
}

# %% 12) ssp_len_RANA

pred.vals.gg %>% filter(sex == "M", bd == 6, mehg == 0) %>% 
mutate(
  pop = plyr::mapvalues(pop, from = unique(pop), to = c(
    "Rana pretiosa
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
) %>% {
  ggplot(., aes(len, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.3) +
    geom_line(size = 1) + 
    scale_colour_discrete() +
    scale_fill_discrete() +
    scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
    facet_wrap(~pop, ncol = 1, strip.position = "right") +
    theme(
      strip.text = element_text(size = 11)
  ,   axis.text.y = element_text(size = 12)
    ) +
    xlab("Individual Length") +
    ylab("Between-Season Survival") 
}

# %% 13) ssp_pop

pop_size_ests %>% mutate(
  pop_spec = plyr::mapvalues(pop_spec, from = unique(pop_spec), to = c(
    "Blackrock-C", "Blackrock-H", "Jones Pond", "Sonoma Mountain", "Two Medicine"
  ))
) %>% {
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
      theme(axis.text.x = element_text(angle = 285, hjust = 0, size = 10)) +
      facet_wrap(~pop_spec, scales = "free", ncol = 2) +
      scale_y_log10() 
}

##### Multiple Species Multiple population fits -----------------------------------------------

# %% 14) msp_Bdmort

ord_of_effort <- sampling %>% group_by(pop_spec) %>% summarize(n_dates = n_distinct(CaptureDate)) %>%
  arrange(desc(n_dates)) %>% mutate(pop_spec = factor(pop_spec, levels = unique(pop_spec)))

pred.vals.gg %>%
  mutate(
  pop = plyr::mapvalues(pop, from = unique(pred.vals.gg$pop)
    , to = c(
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
Scotia Barrens"
, "Notophthalmus viridescens
SMNWR West"
, "Notophthalmus viridescens
SMNWR Springfield"
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
) %>% filter(sex == "M", len == 0, mehg == 0) %>% {
  ggplot(., aes(bd, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = spec, colour = spec), alpha = 0.3) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n, fill = spec, colour = spec), alpha = 0.3) +
    geom_line(aes(colour = spec), size = 1) + 
    scale_colour_brewer(name = "Species", palette = "Dark2") +
    scale_fill_brewer(name = "Species", palette = "Dark2") +
   # scale_x_continuous(breaks = c(0, 3, 6, 9, 12)) +
    scale_x_continuous(breaks = c(-1.5, -0.75, 0, 0.75, 1.5)) +
    facet_wrap(~pop) +
    theme(
      strip.text.x = element_text(size = 11)
  ,   axis.text.y = element_text(size = 12)
  ,   axis.text.x = element_text(size = 11)
    ) +
    xlab("Bd Load") +
    ylab("Between-Season Survival")
}

# %% 15) msp_len_RANA


# %% 16) msp_BCF_AMCI


# %% 17) msp_NOVI


# %% Some priors vs posteriors  
bd_int_prior <- reshape2::melt(stan.fit.samples$beta_offseason_bd)
bd_int_prior %<>% mutate(prior = rnorm(n(), 0, 0.65))
bd_int_prior %<>% mutate(Var2 = plyr::mapvalues(Var2, from = unique(Var2)
  , c("Ambystoma 
cingulatum", "Anaxyrus 
boreas", "Pseudacris 
maculata"
    , "Notophthalmus 
viridescens", "Rana 
spp.")
  )) %>% rename(Species = Var2)

bd_int_prior %<>% rename(Posterior = value, Prior = prior) %>% pivot_longer(-c(iterations, Species))

bd_int_prior %>% {
  ggplot(., aes(x = value)) + 
    geom_density(aes(colour = name), size = 2) +
    scale_colour_brewer(name = "Distribution", palette = "Dark2") +
    theme(
      legend.key.size = unit(0.8, "cm")
    , legend.text = element_text(size = 12)
    , legend.title = element_text(size = 14)
    , strip.text.y = element_text(size = 12)
    ) +
    xlab("Density") +
    ylab("Value") +
    facet_wrap(~Species, ncol = 1, strip.position = "right")
}

# %% Some other priors vs posteriors
bd_rand_prior <- reshape2::melt(stan.fit.samples$offseason_pop_bd)
bd_rand_prior %<>% mutate(prior = rnorm(n(), 0, 0.85))
bd_rand_prior %<>% mutate(Var2 = plyr::mapvalues(Var2, from = unique(Var2)
  , c(
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
    )
  )) %>% rename(Species = Var2)

bd_rand_prior %<>% rename(Posterior = value, Prior = prior) %>% pivot_longer(-c(iterations, Species))

bd_rand_prior %>% {
  ggplot(., aes(x = value)) + 
    geom_density(aes(y = ..scaled.., colour = name), size = 0.5) +
    geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
    scale_colour_brewer(name = "Distribution", palette = "Dark2") +
    theme(
      legend.key.size = unit(0.8, "cm")
    , legend.text = element_text(size = 12)
    , legend.title = element_text(size = 14)
    , strip.text.x = element_text(size = 12)
    ) +
    xlab("Density") +
    ylab("Value") +
    facet_wrap(~Species, scales = "free")
}



