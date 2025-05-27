####
## Model fit exploration and figure making for manuscript 
####

cort_fit <- readRDS("../fits/cmr2.Rds")
samps    <- extract(cort_fit$fitted_model)

cort_fit_summary <- cort_fit$fitted_model %>% summary()

stan_diag(cort_fit$fitted_model)
stan_rhat(cort_fit$fitted_model)
stan_ess(cort_fit$fitted_model)
stan_mcse(cort_fit$fitted_model)

samps$ind_cort_s[, 2] %>% hist(breaks = 200)
samps$beta_cort_r[, 3] %>% hist(breaks = 200)
samps$beta_offseason[, 4] %>% hist(breaks = 200)

color_vec <- c(
  "#7CBA96" ## SVL
, "#3B99B1" ## BL CORT
, "#F5191C" ## SI CORT
, "#E78F0A" ## MeHg
, "#EACB2B" ## Bd
)

samps$beta_offseason %>% reshape2::melt() %>% 
  group_by(Var2) %>% summarize(
    lwr   = quantile(value, 0.025)
    , lwr_n = quantile(value, 0.200)
    , mid   = quantile(value, 0.500)
    , upr_n = quantile(value, 0.800)
    , upr   = quantile(value, 0.975)
  ) %>% mutate(
    Var2 = plyr::mapvalues(Var2, from = c(1, 2, 3, 4, 5), to = c(
      "Bd Load", "Snout-vent length (mm)", "MeHg", "Baseline CORT", "Stress-induced CORT"
    ))
  ) %>% mutate(
    Var2 = factor(Var2, levels = c(
      "Snout-vent length (mm)", "Baseline CORT", "Stress-induced CORT", "MeHg", "Bd Load"
    ) %>% rev())
  ) %>% {
    ggplot(., aes(mid, Var2)) + 
      geom_point(aes(colour = Var2), size = 3) + 
      geom_errorbarh(aes(xmin = lwr, xmax = upr, colour = Var2), height = 0.3) +
      geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n, colour = Var2), height = 0, linewidth = 1.5) +
      scale_colour_manual(values = rev(color_vec)) +
      theme(legend.position="none") +
      ylab("Covariate") +
      xlab("Effect on between-year survival (log odds)") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      annotate("text", x = -0.14, y = 5, label = "96%") +
      annotate("text", x = -0.09, y = 4, label = "98%") +
      annotate("text", x = 0.32, y = 3, label = "84%") +
      annotate("text", x = 0.21, y = 2, label = "93%") +
      annotate("text", x = 0.52, y = 1, label = "51%")
  }

length(which(samps$beta_cort_r[, 1] < 0)) / 8000
length(which(samps$beta_offseason[, 3] < 0)) / 8000

samps$beta_cort_s %>% reshape2::melt() %>% 
  group_by(Var2) %>% summarize(
    lwr   = quantile(value, 0.025)
    , lwr_n = quantile(value, 0.200)
    , mid   = quantile(value, 0.500)
    , upr_n = quantile(value, 0.800)
    , upr   = quantile(value, 0.975)
  ) %>% mutate(
    Var2 = plyr::mapvalues(Var2, from = c(1, 2, 3), to = c(
      "Bd Load", "Snout-vent length (mm)", "MeHg"
    ))
  ) %>% 
  mutate(CORT = "Stress-induced CORT", .after = Var2) %>% 
  rbind(
    .
    , samps$beta_cort_r %>% reshape2::melt() %>% 
      group_by(Var2) %>% summarize(
        lwr   = quantile(value, 0.025)
        , lwr_n = quantile(value, 0.200)
        , mid   = quantile(value, 0.500)
        , upr_n = quantile(value, 0.800)
        , upr   = quantile(value, 0.975)
      ) %>% mutate(
        Var2 = plyr::mapvalues(Var2, from = c(1, 2, 3), to = c(
          "Bd Load", "Snout-vent length (mm)", "MeHg"
        ))
      ) %>% mutate(CORT = "Baseline CORT", .after = Var2)
  ) %>% {
    ggplot(., aes(mid, Var2)) + 
      geom_rect(fill = color_vec[5], xmin = -Inf, xmax = Inf,
                ymin = 0.70, ymax = 1.3, alpha = 0.05) +
      geom_rect(fill = color_vec[4], xmin = -Inf, xmax = Inf,
                ymin = 1.7, ymax = 2.3, alpha = 0.05) +
      geom_rect(fill = color_vec[1], xmin = -Inf, xmax = Inf,
                ymin = 2.7, ymax = 3.3, alpha = 0.05) +
      geom_point(aes(colour = CORT), size = 3, position = position_dodge(width = 0.6)) + 
      geom_errorbarh(aes(xmin = lwr, xmax = upr, colour = CORT)
                     , height = 0.2, position = position_dodge(width = 0.6)) +
      geom_errorbarh(aes(xmin = lwr_n, xmax = upr_n, colour = CORT)
                     , height = 0, linewidth = 1.5, position = position_dodge(width = 0.6)) +
      scale_colour_manual(values = c(color_vec[2], color_vec[3])) +
      ylab("Covariate") +
      xlab("Effect on scaled CORT") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_x_continuous(limits = c(-0.9, 0.7)) +
      annotate("text", x = -0.66, y = 2.85, label = "87%") +
      annotate("text", x = -0.55, y = 1.86, label = "97.2%")
  }

length(which(samps$beta_cort_r[, 3] < 0)) / length(samps$beta_cort_r[, 3])

#mehg_range  <- seq(-3, 3, by = 0.1)
mehg_range  <- seq(-1.7, 3, by = 0.1)
pred_cort.b <- matrix(
  ncol = length(mehg_range)
  , nrow = length(samps$beta_cort_r[, 3])
)
pred_cort.s <- matrix(
  ncol = length(mehg_range)
  , nrow = length(samps$beta_cort_s[, 3])
)

for (i in seq_along(mehg_range)) {
  pred_cort.b[, i] <- samps$beta_cort_r[, 3] * mehg_range[i] + samps$beta_cort_r_sex[, 1]
  pred_cort.s[, i] <- samps$beta_cort_s[, 3] * mehg_range[i] + samps$beta_cort_r_sex[, 1]
}

pred_cort.b.s <- reshape2::melt(pred_cort.b) %>%
  group_by(Var2) %>% summarize(
    lwr        = quantile(value, 0.025)
    , lwr_n      = quantile(value, 0.200)
    , mid        = quantile(value, 0.500)
    , upr_n      = quantile(value, 0.800)
    , upr        = quantile(value, 0.975)
  ) %>% mutate(
    Var2 = plyr::mapvalues(Var2, from = unique(Var2), to = mehg_range)
  ) %>% rename(MeHg = Var2) %>%
  mutate(
    CORT = "Baseline CORT", .after = MeHg
  )

pred_cort.s.s <- reshape2::melt(pred_cort.s) %>%
  group_by(Var2) %>% summarize(
    lwr        = quantile(value, 0.025)
    , lwr_n      = quantile(value, 0.200)
    , mid        = quantile(value, 0.500)
    , upr_n      = quantile(value, 0.800)
    , upr        = quantile(value, 0.975)
  ) %>% mutate(
    Var2 = plyr::mapvalues(Var2, from = unique(Var2), to = mehg_range)
  ) %>% rename(MeHg = Var2) %>%
  mutate(
    CORT = "Stress-induced CORT", .after = MeHg
  )

cort_preds <- rbind(
  pred_cort.b.s, pred_cort.s.s
)

cort_pred.gg <- capt_history %>% 
  group_by(Mark) %>%
  summarize(
    m_bd   = mean(bd_load, na.rm = T)
    , m_merc = mean(merc, na.rm = T)
    , len    = mean(len, na.rm = T)
    , len_raw = mean(len_raw, na.rm = T)
    , Sex    = Sex[1]
    , b_cort = mean(cort_base_conc, na.rm = T)
    , s_cort = mean(cort_stress_conc, na.rm = T)
  ) %>% 
  pivot_longer(
    -c(Mark, m_bd, m_merc, len, len_raw, Sex)
    , names_to = "CORT", values_to = "conc") %>% 
  mutate(CORT = plyr::mapvalues(
    CORT, from = c("b_cort", "s_cort"), to = c("Baseline CORT", "Stress-induced CORT")
  )) %>% group_by(CORT) %>%
  mutate(
    conc_s = scale(log(conc))[, 1]
  #, merc_s = scale(log(m_merc))[, 1]
   , merc_s = scale(m_merc)[, 1]
  ) 

scaled_conc         <- cort_pred.gg %>% filter(CORT == "Baseline CORT")
scaled_conc_attr.bc <- scale(scaled_conc$conc %>% log())
scaled_conc_attr.bc <- c(attr(scaled_conc_attr.bc, 'scaled:scale'), attr(scaled_conc_attr.bc, 'scaled:center'))

x.lab.b  <- c(100, 200, 400, 800, 1600, 3200, 6400)
x.lab.bc <- (log(x.lab.b) - scaled_conc_attr.bc[2]) / scaled_conc_attr.bc[1]

scaled_conc         <- cort_pred.gg %>% filter(CORT == "Stress-induced CORT")
scaled_conc_attr.sc <- scale(scaled_conc$conc %>% log())
scaled_conc_attr.sc <- c(attr(scaled_conc_attr.sc, 'scaled:scale'), attr(scaled_conc_attr.sc, 'scaled:center'))

x.lab.s  <- c(500, 1000, 2000, 4000, 8000, 16000, 32000, 64000)
x.lab.sc <- (log(x.lab.s) - scaled_conc_attr.sc[2]) / scaled_conc_attr.sc[1]

scaled_conc        <- cort_pred.gg %>% filter(CORT == "Stress-induced CORT")
#scaled_conc_attr.m <- scale(scaled_conc$m_merc %>% log())
scaled_conc_attr.m <- scale(scaled_conc$m_merc)
scaled_conc_attr.m <- c(attr(scaled_conc_attr.m, 'scaled:scale'), attr(scaled_conc_attr.m, 'scaled:center'))

#x.lab.m  <- c(40, 80, 160, 320, 640)
x.lab.m  <- c(0, 100, 200, 300, 400, 500, 600)
#x.lab.mc <- (log(x.lab.m) - scaled_conc_attr.m[2]) / scaled_conc_attr.m[1]
x.lab.mc <- (x.lab.m - scaled_conc_attr.m[2]) / scaled_conc_attr.m[1]

scaled_conc        <- cort_pred.gg %>% filter(CORT == "Stress-induced CORT")
#scaled_conc_attr.l <- scale(scaled_conc$len_raw %>% log())
scaled_conc_attr.l <- scale(scaled_conc$len_raw)
scaled_conc_attr.l <- c(attr(scaled_conc_attr.l, 'scaled:scale'), attr(scaled_conc_attr.l, 'scaled:center'))

x.lab.l  <- c(20, 40, 60, 80)
#x.lab.lc <- (log(x.lab.l) - scaled_conc_attr.l[2]) / scaled_conc_attr.l[1]
x.lab.lc <- (x.lab.l - scaled_conc_attr.l[2]) / scaled_conc_attr.l[1]

gg.1 <- cort_pred.gg %>% 
  filter(CORT == "Baseline CORT") %>% {
    ggplot(., aes(merc_s, conc_s)) + 
      geom_line(data = cort_preds %>% rename(conc_s = mid) %>% filter(CORT == "Baseline CORT")
                , aes(x = MeHg, y = conc_s), alpha = 0.5, colour = color_vec[2]) +
      geom_ribbon(data = cort_preds %>% rename(conc_s = mid) %>% filter(CORT == "Baseline CORT")
                  , aes(x = MeHg, ymin = lwr, ymax = upr), alpha = 0.3, fill = color_vec[2]) +
      geom_point(alpha = 0.4, colour = color_vec[2]) +
      scale_y_continuous(breaks = x.lab.bc, labels = x.lab.b) +
      scale_x_continuous(breaks = x.lab.mc, labels = x.lab.m
                         , limits = c(-1.7, 2.7)
                         ) +
      xlab("MeHg") +
      ylab("Baseline CORT (pg/ml)")
  }

gg.2 <- cort_pred.gg %>% 
  filter(CORT == "Stress-induced CORT") %>% {
    ggplot(., aes(merc_s, conc_s)) + 
      geom_line(data = cort_preds %>% rename(conc_s = mid) %>% filter(CORT == "Stress-induced CORT")
                , aes(x = MeHg, y = conc_s), alpha = 0.5, colour = color_vec[3]) +
      geom_ribbon(data = cort_preds %>% rename(conc_s = mid) %>% filter(CORT == "Stress-induced CORT")
                  , aes(x = MeHg, ymin = lwr, ymax = upr), alpha = 0.3, fill = color_vec[3]) +
      geom_point(alpha = 0.4, colour = color_vec[3]) +
      scale_y_continuous(breaks = x.lab.sc, labels = x.lab.s) +
      scale_x_continuous(breaks = x.lab.mc, labels = x.lab.m
                         , limits = c(-1.7, 2.7)
                         ) +
      xlab("MeHg") +
      ylab("Stress-induced CORT (pg/ml)") 
  } 

gridExtra::grid.arrange(gg.1, gg.2, ncol = 1)

ttt <- capt_history %>% group_by(Mark) %>%
  summarize(
    m_bd   = mean(bd_load, na.rm = T)
    , m_merc = mean(merc, na.rm = T)
    , len    = mean(len, na.rm = T)
    , Sex    = Sex[1]
    , m_cort = mean(cort_base_conc, na.rm = T)
  )

ttt$m_bd   <- scale(log(ttt$m_bd + 1))[, 1]
ttt$m_merc <- scale(ttt$m_merc)[, 1]
ttt$len    <- scale(ttt$len)
ttt$m_cort <- scale(ttt$m_cort)[, 1]

lm(m_cort ~ m_bd + m_merc + len + Sex, data = ttt) %>% summary()

## Figure for impact of range of cort on survival

#ind_sex[ind_occ_min1_rep[phi_off_index], ] * beta_offseason_sex + 
#beta_offseason[1] * X_scaled[phi_bd_index[phi_off_index]] +
#beta_offseason[2] * ind_len_have[ind_occ_min1_rep[phi_off_index]] +
#beta_offseason[3] * ind_mehg_scaled[ind_occ_min1_rep[phi_off_index]] +
#beta_offseason[4] * ind_cort_r[ind_occ_min1_rep[phi_off_index]] +
#beta_offseason[5] * ind_cort_s[ind_occ_min1_rep[phi_off_index]]

color_vec <- c(
  "#7CBA96" ## SVL
  , "#3B99B1" ## BL CORT
  , "#F5191C" ## SI CORT
  , "#E78F0A" ## MeHg
  , "#EACB2B" ## Bd
)

######
### Resting Cort
######

cort_range <- seq(-3, 3, by = 0.1)
cort_est   <- data.frame(
  cort_range = cort_range
  , lwr        = 0
  , lwr_n      = 0
  , mid        = 0
  , upr_n      = 0
  , upr        = 0
)
for (i in seq_along(cort_range)) {
  temp_hist        <- samps$beta_offseason_sex[, 1] + 
    samps$beta_offseason[, 2] * 0 +
    samps$beta_offseason[, 4] * cort_range[i]
  temp_quant       <- quantile(temp_hist, c(0.025, 0.200, 0.500, 0.800, 0.975))
  cort_est[i, 2:6] <- temp_quant %>% plogis()
}

cort_map <- capt_history %>% 
  dplyr::select(cort_base_conc) %>% 
  filter(!is.na(cort_base_conc)) %>% 
  mutate(cort_scaled = scale(log(cort_base_conc))[, 1] %>% round(1)) %>%
  group_by(cort_scaled) %>%
  summarize(cort_base_conc = mean(cort_base_conc)) %>%
  ungroup() %>%
  rename(cort_range = cort_scaled)

cort_est.gg <- cort_est %>% mutate(cort_range = round(cort_range, 1)) %>% left_join(., cort_map)

cort_hist <- capt_history %>% 
  group_by(Mark) %>%
  slice(1) %>%
  dplyr::select(cort_base_conc) %>% 
  filter(!is.na(cort_base_conc)) %>% 
  ungroup() %>%
  mutate(cort_scaled = scale(log(cort_base_conc))[, 1] %>% round(1))

gg.1 <- cort_est.gg %>% {
  ggplot(., aes(cort_range, mid)) +
    geom_dotplot(
      data = cort_hist %>% mutate(mid = 0), aes(x = cort_scaled)
      , binwidth = 0.18, method = 'histodot', alpha = 0.25, fill = color_vec[2], colour = NA
    ) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.3, fill = color_vec[2]) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, fill = color_vec[2]) +
    geom_line(colour = color_vec[2]) +
    xlab("Baseline CORT (pg/ml)") +
    ylab("Survival probability between years") +
    scale_x_continuous(
      breaks = x.lab.bc, labels = x.lab.b
    ) +
    scale_y_continuous(
      limits = c(0, 1)
    )
}

######
### Stress Cort
######

cort_range <- seq(-3, 3, by = 0.1)
cort_est   <- data.frame(
  cort_range = cort_range
  , lwr        = 0
  , lwr_n      = 0
  , mid        = 0
  , upr_n      = 0
  , upr        = 0
)
for (i in seq_along(cort_range)) {
  temp_hist        <- samps$beta_offseason_sex[, 1] + 
    samps$beta_offseason[, 2] * 0 +
    samps$beta_offseason[, 5] * cort_range[i]
  temp_quant       <- quantile(temp_hist, c(0.025, 0.200, 0.500, 0.800, 0.975))
  cort_est[i, 2:6] <- temp_quant %>% plogis()
}

cort_map <- capt_history %>% 
  dplyr::select(cort_stress_conc) %>% 
  filter(!is.na(cort_stress_conc)) %>% 
  mutate(cort_scaled = scale(log(cort_stress_conc))[, 1] %>% round(1)) %>%
  group_by(cort_scaled) %>%
  summarize(cort_stress_conc = mean(cort_stress_conc)) %>%
  ungroup() %>%
  rename(cort_range = cort_scaled)

cort_est.gg <- cort_est %>% mutate(cort_range = round(cort_range, 1)) %>% left_join(., cort_map)

cort_hist <- capt_history %>% 
  group_by(Mark) %>%
  slice(1) %>%
  dplyr::select(cort_stress_conc) %>% 
  filter(!is.na(cort_stress_conc)) %>% 
  ungroup() %>%
  mutate(cort_scaled = scale(log(cort_stress_conc))[, 1] %>% round(1))

gg.2 <- cort_est.gg %>% {
  ggplot(., aes(cort_range, mid)) +
    geom_dotplot(
      data = cort_hist %>% mutate(mid = 0), aes(x = cort_scaled)
      , binwidth = 0.15, method = 'histodot', alpha = 0.25, fill = color_vec[3], colour = NA
    ) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.3, fill = color_vec[3]) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, fill = color_vec[3]) +
    geom_line(colour = color_vec[3]) +
    xlab("Stress-induced CORT (pg/ml)") +
    ylab("") +
    theme(axis.text.x = element_text(size = 10)) +
    scale_x_continuous(
      breaks = x.lab.sc, labels = x.lab.s
    ) +
    scale_y_continuous(
      limits = c(0, 1)
    )
}

######
### MeHg
######

mehg_range <- seq(-1.7, 3, by = 0.01)
mehg_est   <- data.frame(
  merc_range = mehg_range
  , lwr        = 0
  , lwr_n      = 0
  , mid        = 0
  , upr_n      = 0
  , upr        = 0
)
for (i in seq_along(mehg_range)) {
  temp_hist        <- samps$beta_offseason_sex[, 1] + 
    samps$beta_offseason[, 2] * 0 +
    samps$beta_offseason[, 3] * mehg_range[i]
  temp_quant       <- quantile(temp_hist, c(0.025, 0.200, 0.500, 0.800, 0.975))
  mehg_est[i, 2:6] <- temp_quant %>% plogis()
}

mehg_map <- capt_history %>% 
  dplyr::select(merc) %>% 
  filter(!is.na(merc)) %>% 
#  mutate(merc_scaled = scale(log(merc))[, 1] %>% round(2)) %>%
  mutate(merc_scaled = scale(merc)[, 1] %>% round(2)) %>%
  group_by(merc_scaled) %>%
  summarize(merc = mean(merc)) %>%
  ungroup() %>%
  rename(merc_range = merc_scaled)

#scaled_mehg_attr <- capt_history %>% 
#  dplyr::select(merc) %>% 
#  filter(!is.na(merc)) %>% 
#  pull(merc) %>% scale()
#scaled_conc_attr.s <- c(attr(scaled_mehg_attr, 'scaled:scale'), attr(scaled_mehg_attr, 'scaled:center'))

mehg_est.gg <- mehg_est %>% 
  mutate(merc_range = round(merc_range, 3)) %>% 
  left_join(., mehg_map)# %>% 
 # mutate(
 #   merc2 = scaled_conc_attr.s[1] * merc_range + scaled_conc_attr.s[2]
 # ) %>% dplyr::select(-merc) %>% 
 # rename(MeHg_scaled = merc_range, MeHg = merc2) %>%
 # relocate(MeHg, .after = MeHg_scaled) %>%
 # filter(MeHg_scaled >= -1.5) %>% 
 # dplyr::select(-MeHg_scaled)

mehg_hist <- capt_history %>% 
  group_by(Mark) %>%
  slice(1) %>%
  dplyr::select(merc) %>% 
  filter(!is.na(merc)) %>% 
  ungroup() %>%
#  mutate(merc_scaled = scale(log(merc))[, 1] %>% round(1))
  mutate(merc_scaled = scale(merc)[, 1] %>% round(1))

#x.lab.m  <- c(80, 160, 320, 640, 1280, 2560)
#x.lab.mc <- (log(x.lab.m) - scaled_conc_attr.m[2]) / scaled_conc_attr.m[1]

gg.3 <- mehg_est.gg %>% {
  ggplot(., aes(mehg_range, mid)) +
    geom_dotplot(
      data = mehg_hist %>% mutate(mid = 0), aes(x = merc_scaled)
      , binwidth = 0.096, method = 'histodot'
      , alpha = 0.25, fill = color_vec[4], colour = NA
    ) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.3, fill = color_vec[4]) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3, fill = color_vec[4]) +
    geom_line(colour = color_vec[4]) +
    xlab("MeHg (ng/g dw)") +
    ylab("Survival probability between years") +
    scale_x_continuous(
      breaks = x.lab.mc, labels = x.lab.m, limits = c(-1.7, 2.7)
    ) +
    scale_y_continuous(
      limits = c(0, 1)
    )
}


######
### Length
######

len_range <- seq(-3, 3, by = 0.1)
len_est   <- data.frame(
  len_range = len_range
  , lwr        = 0
  , lwr_n      = 0
  , mid        = 0
  , upr_n      = 0
  , upr        = 0
)
for (i in seq_along(len_range)) {
  temp_hist        <- samps$beta_offseason_sex[, 1] + 
    samps$beta_offseason[, 2] * len_range[i] 
  temp_quant       <- quantile(temp_hist, c(0.025, 0.200, 0.500, 0.800, 0.975))
  len_est[i, 2:6] <- temp_quant %>% plogis()
}

len_map <- capt_history %>% 
  dplyr::select(len) %>% 
  filter(!is.na(len)) %>% 
  group_by(len) %>%
  summarize(len = mean(len)) %>%
  ungroup() %>%
  rename(len_range = len)

len_est.gg <- len_est %>% mutate(len_range = round(len_range, 1)) %>% left_join(., len_map)

len_hist <- capt_history %>% 
  group_by(Mark) %>%
  slice(1) %>%
  dplyr::select(len) %>% 
  filter(!is.na(len)) %>% 
  ungroup() %>%
  mutate(len_scaled = len %>% round(1))

#scaled_conc        <- cort_pred.gg %>% filter(CORT == "Stress-induced CORT")
#scaled_conc_attr.l <- scale(scaled_conc$len_raw)
#scaled_conc_attr.l <- c(attr(scaled_conc_attr.l, 'scaled:scale'), attr(scaled_conc_attr.l, 'scaled:center'))

#x.lab.l  <- c(40, 60, 80)
#x.lab.lc <- (x.lab.l - scaled_conc_attr.l[2]) / scaled_conc_attr.l[1]

gg.4 <- len_est.gg %>% {
  ggplot(., aes(len_range, mid)) +
    geom_dotplot(
      data = len_hist %>% mutate(mid = 0), aes(x = len_scaled)
      , binwidth = 0.090, method = 'histodot', alpha = 0.3, fill = color_vec[1], colour = NA
    ) +
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n), alpha = 0.3, fill = color_vec[1]) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4, fill = color_vec[1]) +
    geom_line(colour = color_vec[1]) +
    xlab("Snout-vent length (mm)") +
    ylab("") +
    scale_x_continuous(
      breaks = x.lab.lc, labels = x.lab.l
      , limits = c(-1.8, 2.2)
    ) +
    scale_y_continuous(
      limits = c(0, 1)
    )
}

gridExtra::grid.arrange(
  gg.1 + theme(
    axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
    , axis.title.x = element_text(size = 14)
    , axis.title.y = element_text(size = 12)
  )
  , gg.2 + theme(
    axis.text.x = element_text(size = 9)
    , axis.text.y = element_text(size = 12)
    , axis.title.x = element_text(size = 14)
    , axis.title.y = element_text(size = 14)
  )
  , gg.3 + theme(
    axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
    , axis.title.x = element_text(size = 14)
    , axis.title.y = element_text(size = 12)
  )
  , gg.4 + theme(
    axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
    , axis.title.x = element_text(size = 14)
    , axis.title.y = element_text(size = 14)
  )
)

### Supplemental Figure about repeated measures of cort

source("JP_data.R")

gg.s.1 <- data.all %>% 
  filter(!is.na(cort_base_conc)) %>% 
  group_by(Mark) %>% 
  mutate(n_entry = n()) %>%
  ungroup() %>% 
  filter(n_entry > 1) %>% 
  mutate(Individual = as.factor(Mark) %>% as.numeric() %>% as.factor()) %>%
  ungroup() %>%
  mutate(Year = as.factor(Year)) %>% {
    ggplot(., aes(cort_base_conc, Individual)) + 
      geom_point(aes(colour = Year, shape = Year)) +
      scale_colour_brewer(palette = "Dark2") +
      geom_line() +
      scale_x_log10() +
      xlab("Baseline CORT (pg/ml)") +
      theme(axis.text.y = element_blank())
  }

gg.s.2 <- data.all %>% 
  filter(!is.na(cort_stress_conc)) %>% 
  group_by(Mark) %>% 
  mutate(n_entry = n()) %>%
  ungroup() %>% 
  filter(n_entry > 1) %>% 
  mutate(Individual = as.factor(Mark) %>% as.numeric() %>% as.factor()) %>%
  ungroup() %>%
  mutate(Year = as.factor(Year)) %>% {
    ggplot(., aes(cort_stress_conc, Individual)) + 
      geom_point(aes(colour = Year, shape = Year)) +
      scale_colour_brewer(palette = "Dark2") +
      geom_line() +
      scale_x_log10() +
      xlab("Stress-induced CORT (pg/ml)") +
      theme(axis.text.y = element_blank())
  }

gridExtra::grid.arrange(gg.s.1, gg.s.2, ncol = 1)


### Misc Statements about captures

data.all %>% group_by(Mark) %>% summarize(n_year = n_distinct(Year)) %>% 
  ungroup() %>% group_by(n_year) %>% summarize(n_ent = n())


### Table of coefficients

## bd, len, mehg
quantile(samps$beta_cort_s[, 1], c(0.025, 0.5, 0.975))
length(which(samps$beta_cort_s[, 1] < 0)) / length(samps$beta_cort_s[, 1])

quantile(samps$beta_cort_s_sex[, 3], c(0.025, 0.5, 0.975))
length(which(samps$beta_cort_s_sex[, 3] > 0)) / length(samps$beta_cort_s_sex[, 3])

## bd, len, mehg, cort_r, cort_s
quantile(samps$beta_offseason[, 5], c(0.025, 0.5, 0.975))
length(which(samps$beta_offseason[, 5] < 0)) / length(samps$beta_offseason[, 5])

quantile(samps$beta_offseason_sex[, 3], c(0.025, 0.5, 0.975))
length(which(samps$beta_offseason_sex[, 3] < 0)) / length(samps$beta_offseason_sex[, 3])

## beta_p_sex, beta_p_len
quantile(samps$beta_p_sex[, 3], c(0.025, 0.5, 0.975))
length(which(samps$beta_p_sex[, 3] < 0)) / length(samps$beta_p_sex[, 3])

quantile(samps$beta_p_len, c(0.025, 0.5, 0.975))
length(which(samps$beta_p_len < 0)) / length(samps$beta_p_len)