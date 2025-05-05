######################################################
## A place to store code for all manuscript figures ##
######################################################

#### Supplemental Newt Quadratic Figure ----

capt_history.bd_load %>% filter(pop_spec == "NOVI.ScotiaBarrens") %>% {
ggplot(., aes(julian, bd_load)) + 
  geom_point(alpha = 0.1) +
  geom_line(aes(group = Mark), alpha = 0.3) +
  facet_wrap(~Year, ncol = 1) +
  scale_y_log10() +
  xlab("Julian Day") +
  ylab("Bd Load (Copies/Swab)") +
  theme(legend.key.size = unit(0.8, "cm"))
}

#### Checking gaps between sampling occasions in all populations ----

## Note: Run after data_manip.R

gg.ppsp1 <- expand.grid(
  pop_spec   = unique(capt_history$pop_spec)
, capt_gaps  = seq(365)
) %>% left_join(
  .
, capt_history %>% 
  group_by(pop_spec, Year) %>% 
  summarize(capt_gaps = unique(capture_gap)) %>%
  mutate(gap_exists = 1) 
) %>% filter(!is.na(gap_exists)) %>% {
  ggplot(., aes(capt_gaps, pop_spec)) + geom_point(aes(colour = as.factor(Year)), size = 3, alpha = 0.7) +
    scale_color_brewer(palette = "Dark2", name = "Year", guide = 'none') +
    geom_hline(yintercept = 2.5, linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = 7.5, linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = 9.5, linetype = "dashed", size = 0.5) + 
    geom_hline(yintercept = 14.5, linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = 82, size = 0.75, colour = "dodgerblue4", linetype = "dashed") +
    geom_vline(xintercept = 135, size = 0.75, colour = "firebrick4") +
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
    xlab("Days Bewtween Sampling Occasions") +
    ylab("Population")
}

gg.ppsp2 <- expand.grid(
  pop_spec   = unique(capt_history$pop_spec)
, capt_gaps  = seq(365)
) %>% left_join(
  .
, capt_history %>% 
  group_by(pop_spec, Year) %>% 
  summarize(capt_gaps = unique(capture_gap)) %>%
  mutate(gap_exists = 1) 
) %>% filter(!is.na(gap_exists), capt_gaps < 82) %>% {
  ggplot(., aes(capt_gaps, pop_spec)) + geom_point(aes(colour = as.factor(Year)), size = 3, alpha = 0.7) +
    scale_color_brewer(palette = "Dark2", name = "Year") +
    geom_hline(yintercept = 2.5, linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = 7.5, linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = 9.5, linetype = "dashed", size = 0.5) + 
    geom_hline(yintercept = 14.5, linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = 9, size = 0.5, colour = "firebrick2") +
    xlab("Days Bewtween Sampling Occasions") +
    theme(
      axis.text.y = element_blank()
    , plot.margin = unit(c(0.2, 4, 0.2, 0.2), "cm")) + 
    ylab("")
}

gridExtra::grid.arrange(gg.ppsp1, gg.ppsp2, ncol = 2)

#### Exploring sampling schemes in specific populations ----

samp_days <- sampling %>% ungroup()

samp_days %<>% 
  rename(J_date = julian) %>% 
  dplyr::select(pop_spec, Species, CaptureDate, Year, Month, J_date) %>%
  mutate(J_date.deg     = (360 / 365) * J_date) %>% 
  mutate(y_coord = (Year - min(Year))) 

between_season_segment <- capt_history.phi %>% 
  filter(offseason == 1) %>% 
  group_by(pop_spec, capture_date) %>%
  slice(1) %>% as.data.frame() %>% 
  mutate(J_date = as.POSIXlt(capture_date)$yday) %>%
  dplyr::select(pop_spec, sampled, capture_date, Year, Month, J_date, capture_gap) %>%
  mutate(J_date_end = J_date + capture_gap) %>%
  mutate(
    J_date.deg     = (360 / 365) * J_date
  , J_date_end.deg = (360 / 365) * J_date_end
  ) %>% 
  mutate(
    y_coord = (Year - min(samp_days$Year))
  ) %>%
  mutate(
    y_end = ifelse(J_date_end.deg > 360, y_coord + 1, y_coord)
  ) 

between_season_segment.t <- between_season_segment %>% filter(y_coord != y_end)
between_season_segment.n <- between_season_segment %>% filter(y_coord == y_end)

if (nrow(between_season_segment.t) > 1) {

for (i in 1:nrow(between_season_segment.t)) {
  temp_entry <- between_season_segment.t[i, ]
  temp_entry <- rbind(temp_entry, temp_entry)
  temp_entry[2, ]$J_date_end.deg <- temp_entry[1, ]$J_date_end.deg - 360
  temp_entry[2, ]$J_date.deg     <- 0
  temp_entry[1, ]$J_date_end.deg <- 360
  temp_entry[2, ]$y_coord <- temp_entry[2, ]$y_end
  
  if (i == 1) {
    temp_entry.f <- temp_entry
  } else {
    temp_entry.f <- rbind(temp_entry.f, temp_entry)
  }
}

between_season_segment <- rbind(between_season_segment.n, temp_entry.f)

}

## Quick manual special consideration for ANBO.SonomaMountain which wasn't sampled in 2021 but was in 2022
between_season_segment.sm <- between_season_segment %>% filter(pop_spec == "ANBO.SonomaMountain")
between_season_segment    <- between_season_segment %>% filter(pop_spec != "ANBO.SonomaMountain")
between_season_segment.sm <- rbind(between_season_segment.sm, between_season_segment.sm[6, ])
between_season_segment.sm[6, ]$J_date_end.deg <- 360
between_season_segment.sm[6, ]$y_end <- between_season_segment.sm[6, ]$y_end + 1
between_season_segment.sm[7, ]$y_coord <- between_season_segment.sm[7, ]$y_coord + 1
between_season_segment.sm[7, ]$y_end <- between_season_segment.sm[7, ]$y_end + 1
between_season_segment.sm[7, ]$J_date_end.deg <- between_season_segment.sm[7, ]$J_date_end.deg - 360
between_season_segment <- rbind(between_season_segment, between_season_segment.sm)

within_season_segment <- capt_history.phi %>% 
  filter(phi_ones == 0 & offseason == 0) %>% 
  group_by(pop_spec, capture_date) %>%
  slice(1) %>% as.data.frame() %>% 
  mutate(J_date = as.POSIXlt(capture_date)$yday) %>%
  dplyr::select(pop_spec, sampled, capture_date, Year, Month, J_date, capture_gap) %>%
  mutate(J_date_end = J_date + capture_gap) %>%
  mutate(
    J_date.deg     = (360 / 365) * J_date
  , J_date_end.deg = (360 / 365) * J_date_end
  ) %>% 
  mutate(
    y_coord = (Year - min(samp_days$Year))
  ) %>%
  mutate(
    y_end = ifelse(J_date_end.deg > 360, y_coord + 1, y_coord)
  ) 

within_season_segment.t <- within_season_segment %>% filter(y_coord != y_end)
within_season_segment.n <- within_season_segment %>% filter(y_coord == y_end)

if (nrow(within_season_segment.t) > 1) {

for (i in 1:nrow(within_season_segment.t)) {
  temp_entry <- within_season_segment.t[i, ]
  temp_entry <- rbind(temp_entry, temp_entry)
  temp_entry[2, ]$J_date_end.deg <- temp_entry[1, ]$J_date_end.deg - 360
  temp_entry[2, ]$J_date.deg     <- 0
  temp_entry[1, ]$J_date_end.deg <- 360
  temp_entry[2, ]$y_coord <- temp_entry[2, ]$y_end
  
  if (i == 1) {
    temp_entry.f <- temp_entry
  } else {
    temp_entry.f <- rbind(temp_entry.f, temp_entry)
  }
}

within_season_segment <- rbind(within_season_segment.n, temp_entry.f)

}

closed_season_segment <- capt_history.phi %>% 
  filter(phi_ones == 1) %>% 
  group_by(pop_spec, capture_date) %>%
  slice(1) %>% as.data.frame() %>% 
  mutate(J_date = as.POSIXlt(capture_date)$yday) %>%
  dplyr::select(pop_spec, sampled, capture_date, Year, Month, J_date, capture_gap) %>%
  mutate(J_date_end = J_date + capture_gap) %>%
  mutate(
    J_date.deg     = (360 / 365) * J_date
  , J_date_end.deg = (360 / 365) * J_date_end
  ) %>% 
  mutate(
    y_coord = (Year - min(samp_days$Year)) 
  ) %>%
  mutate(
    y_end = ifelse(J_date_end.deg > 360, y_coord + 1, y_coord)
  ) 

closed_season_segment.t <- closed_season_segment %>% filter(y_coord != y_end)
closed_season_segment.n <- closed_season_segment %>% filter(y_coord == y_end)

if (nrow(closed_season_segment.t) > 1) {

for (i in 1:nrow(closed_season_segment.t)) {
  temp_entry <- closed_season_segment.t[i, ]
  temp_entry <- rbind(temp_entry, temp_entry)
  temp_entry[2, ]$J_date_end.deg <- temp_entry[1, ]$J_date_end.deg - 360
  temp_entry[2, ]$J_date.deg     <- 0
  temp_entry[1, ]$J_date_end.deg <- 360
  temp_entry[2, ]$y_coord <- temp_entry[2, ]$y_end
  
  if (i == 1) {
    temp_entry.f <- temp_entry
  } else {
    temp_entry.f <- rbind(temp_entry.f, temp_entry)
  }
}

closed_season_segment <- rbind(closed_season_segment.n, temp_entry.f)

}

which_pops_to_plot <- unique(data.all$pop_spec)

samp_days %>% 
  filter(pop_spec %in% which_pops_to_plot) %>% {
  ggplot(., aes(J_date.deg, y_coord)) + 
    geom_point(size = 1.0, alpha = 0.5) + 
    coord_polar() + 
    xlab("") + 
    ylab("") + 
    scale_y_continuous(
      lim    = c(0.5, 5.1)
     , breaks = unique((samp_days %>% 
         filter(pop_spec %in% which_pops_to_plot))$Year) %>% length() %>% seq()
    , labels = c(
        unique((samp_days %>% filter(pop_spec %in% which_pops_to_plot))$Year) %>%
        sort()
      )
      ) +
    scale_x_continuous(
      breaks = c(2, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
    , labels = c(
      "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"
    , "Sep", "Oct", "Nov", "Dec"
     )) +
    geom_segment(
      data = between_season_segment %>% filter(pop_spec %in% which_pops_to_plot)
    , aes(
      x    = J_date.deg
    , xend = J_date_end.deg
    , y    = y_coord
    , yend = y_end
      )
      , colour = "dodgerblue3" 
    ) +
    geom_segment(
      data = within_season_segment %>% filter(pop_spec %in% which_pops_to_plot)
    , aes(
      x    = J_date.deg
    , xend = J_date_end.deg
    , y    = y_coord
    , yend = y_end
      )
      , colour = "firebrick3" 
    ) +
    geom_segment(
      data = closed_season_segment %>% filter(pop_spec %in% which_pops_to_plot)
    , aes(
      x    = J_date.deg
    , xend = J_date_end.deg
    , y    = y_coord
    , yend = y_end
      )
      , colour = "springgreen2" 
    ) +
    facet_wrap(~pop_spec, nrow = 4) +
      theme(
      axis.text.x = element_text(size = 10)
    , axis.text.y = element_text(size = 11)
    , axis.ticks  = element_blank()
    , strip.text.x = element_text(size = 10)
    , plot.margin = unit(c(.1,.1,.1,.1), "cm")
    )
  }


#### Bd over time and temperature ----

## Broadest view of Bd loads across temperatures
capt_history.p %>% filter(swabbed == 1) %>% {
  ggplot(., aes(cumtemp_m, log_bd_load)) +
    geom_point(aes(colour = as.factor(Year))) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    facet_wrap(~Site, scales = "free") +
    theme(
      axis.text.x = element_text(size = 10)
    , axis.text.y = element_text(size = 10)
    ) 
}

## Broken down by a species
capt_history.p %>% filter(swabbed == 1) %>%
  filter(Species == "RANA") %>% {
  ggplot(., aes(cumtemp_m_s, log_bd_load)) +
    geom_jitter(aes(colour = as.factor(Year))) +
    geom_line(aes(group = Mark), alpha = 0.3) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    theme(
      axis.text.x = element_text(
        size = 10
      , angle = 300, hjust = 0)
    , axis.text.y = element_text(size = 10)
    ) +
    facet_grid(~Site, scales = "free") +
    xlab("Cumulative temperature") +
    ylab("Log Bd Load") +
  ggtitle("Rana Species")
  }

## Broken down by a species and year
capt_history.p %>% filter(swabbed == 1) %>%
  filter(Species == "RANA") %>% {
  ggplot(., aes(cumtemp_u_s, log_bd_load)) +
    geom_point(aes(colour = as.factor(Year))) +
    geom_line(aes(group = Mark), alpha = 0.5) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    theme(
      axis.text.x = element_text(
        size = 10
      , angle = 300, hjust = 0)
    , axis.text.y = element_text(size = 10)
    ) +
    facet_grid(Year~Site, scales = "free") +
    xlab("Cumulative temperature") +
    ylab("Log Bd Load") +
  ggtitle("Rana Species")
  }

## and also load vs proportion infected
capt_history.p %>% filter(swabbed == 1) %>%
  mutate(infected = ifelse(log_bd_load > 0, 1, 0)) %>%
  filter(Species == "RANA") %>% {
  ggplot(., aes(yday, log_bd_load)) +
    geom_jitter(aes(colour = as.factor(Year))) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    theme(
      axis.text.x = element_text(
        size = 10
      , angle = 300, hjust = 0)
    , axis.text.y = element_text(size = 10)
    ) +
    facet_grid(Year~Site, scales = "free") +
    xlab("Julian Day") +
    ylab("Log Bd Load") +
  ggtitle("Rana Species")
  }

## FIGURE S4

## and also proportion infected 
capt_history.bd_load %>% filter(swabbed == 1) %>%
  mutate(infected = ifelse(bd_load > 0, 1, 0)) %>%
  group_by(pop_spec, capture_date
    , Year, Site, Species, julian) %>%
  summarize(prop_inf = sum(infected) / n()) %>%
  filter(
    pop_spec %in% c(
      "RANA.LostHorse", "RANA.SanFrancisquito"
    , "NOVI.KettleMoraine", "NOVI.Springfield"
    , "AMCI.SMNWR_W", "PSMA.MatthewsPond"
    , "ANBO.JonesPond", "ANBO.Blackrock"
    )
  ) %>% mutate(
    pop_spec = plyr::mapvalues(pop_spec
      , from = c(
      "RANA.LostHorse", "RANA.SanFrancisquito"
    , "NOVI.KettleMoraine", "NOVI.Springfield"
    , "AMCI.SMNWR_W", "PSMA.MatthewsPond"
    , "ANBO.JonesPond", "ANBO.Blackrock"
    )
    , to = c(
      "RALU
Lost Horse"
    , "RADR
San Francisquito"
    , "NOVI
Kettle Moraine"
    , "NOVI
Springfield"
    , "AMCI
SMNWR W"
    , "PSMA
Matthews Pond"
    , "ANBO
Jones Pond"
    , "ANBO
Blackrock"
    )
    )
  ) %>% {
  ggplot(., aes(julian, prop_inf)) +
    geom_point(aes(colour = as.factor(Year)), size = 2) +
    geom_line() +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    theme(
      axis.text.x = element_text(
        size = 10
      , angle = 300, hjust = 0)
    , axis.text.y = element_text(size = 10)
    , strip.text.x = element_text(size = 11)
    ) +
  # facet_wrap(~pop_spec, scales = "free", nrow = 1) +
    facet_grid(Year~pop_spec, scales = "free") +
    xlab("Julian Day") +
    ylab("Proportion Infected")
  }

### FIGURE S2 and FIGURE S3

## Checking Rana loads vs dates and temps
RANAgg1 <- capt_history.bd_load %>% 
  filter(swabbed == 1) %>%
  filter(Species == "ANBO") %>% {
  ggplot(., aes(
      julian
    , bd_load)) +
    geom_point(aes(colour = as.factor(Year)), alpha = 0.4) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    scale_y_log10() +
    theme(
      axis.text.y = element_text(size = 12)
    , axis.text.x = element_text(size = 10)
    , axis.title.x = element_text(size = 12)
    , strip.text.x = element_text(size = 10)
    , plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
    ) +
    scale_x_continuous(lim = c(0, 365)) +
    facet_wrap(~Site, nrow = 1) +
    xlab("Julian Day") +
    ylab("Log Bd Load")
  }; RANAgg1

## and temperature of each of these populations 
RANAgg2 <- temp.all %>%
  filter(Site %in% 
#  c(
#    "DilmanMeadows", "FoxCreek", "JonesPond"
#  , "LostHorse", "SanFrancisquito", "SummitMeadow"
#  , "ThreeCreeks"
#  )
  c(
    "Blackrock", "JonesPond"
  , "SonomaMountain", "TwoMedicine"
  )
#  c(
#    "EmmaCarlin", "MudLake", "ScotiaBarrens"
#  , "SMNWR_W", "SPR"
#  )
#  c(
#   "SanFrancisquito", "ScotiaBarrens"
# )
    ) %>% filter(Year != 2017) %>% {
  ggplot(., aes(yday, tmax)) +
    geom_line(aes(colour = as.factor(Year)), alpha = 0.5) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    theme(
      axis.text.x = element_text(size = 10)
    , axis.text.y = element_text(size = 10)
    , strip.text.x = element_blank()
    , plot.margin = unit(c(0.2, 0.2, 0.2, 0.84), "cm")
    ) +
    facet_wrap(~Site, nrow = 1) +
    geom_hline(yintercept = 17, linetype = "dashed", colour = "firebrick3") +
    geom_hline(yintercept = 25, linetype = "dashed", colour = "firebrick3") +
    xlab("Julian Day") +
    ylab("Mean Daily Temperature")
}; RANAgg2

gridExtra::grid.arrange(RANAgg1, RANAgg2, ncol = 1)

## All pops over time
capt_history.p %>% filter(swabbed == 1) %>%
  filter(Year != 2017) %>% {
  ggplot(., aes(yday_s, log_bd_load)) +
    geom_line(aes(group = Mark), alpha = 0.5) +
    geom_point(aes(colour = as.factor(Year)), alpha = 0.5) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    facet_grid(Year~Site) +
    theme(
      axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
      ) +
    xlab("Julian Day") +
    ylab("Log Bd Load")
}

## Focusing just on Scotia Barrens for temp
SBgg.1 <- capt_history.p %>% filter(swabbed == 1) %>%
  filter(Year != 2017) %>%
  filter(Site == "ScotiaBarrens") %>% {
  ggplot(., aes(cumtemp, log_bd_load)) +
    geom_line(aes(group = Mark), alpha = 0.5) +
    geom_point(aes(colour = as.factor(Year)), alpha = 0.5) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    facet_grid(~Year) +
    theme(
      axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
    , plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
      ) +
    xlab("Cumulative Temperature") +
    ylab("Log Bd Load") +
    ggtitle("Scotia Barrens")
}

SBgg.2 <- capt_history.p %>% filter(swabbed == 1) %>%
  filter(Year != 2017) %>%
  filter(Site == "ScotiaBarrens") %>% {
  ggplot(., aes(yday, log_bd_load)) +
    geom_line(aes(group = Mark), alpha = 0.5) +
    geom_point(aes(colour = as.factor(Year)), alpha = 0.5) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    facet_grid(~Year) +
    theme(
      axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
    , plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
      ) +
    xlab("Julian Day") +
    ylab("Log Bd Load") 
  }

SBgg.3 <- capt_history.p %>% filter(swabbed == 1) %>%
  filter(Year != 2017) %>%
  filter(Site == "ScotiaBarrens") %>% {
  ggplot(., aes(days_in_opt_m, log_bd_load)) +
    geom_line(aes(group = Mark), alpha = 0.5) +
    geom_point(aes(colour = as.factor(Year)), alpha = 0.5) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    facet_grid(~Year) +
    theme(
      axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
      ) +
    xlab("Days within Bd's thermal optima (17-25)") +
    ylab("Log Bd Load") 
  }

SBgg.4 <- temp_data.all %>% 
  filter(
    Year != 2017
  , Site == "ScotiaBarrens") %>% {
  ggplot(., aes(yday, tmean)) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    geom_line() +
    facet_grid(~Year) +
    theme(
      axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
      ) +
    geom_hline(yintercept = 17, linetype = "dashed", colour = "firebrick3") +
    geom_hline(yintercept = 25, linetype = "dashed", colour = "firebrick3") +
    xlab("Julian Day") +
    ylab("Daily average temperature") 
}

gridExtra::grid.arrange(SBgg.1, SBgg.2, ncol = 1)

####
## Bd vs MeHg in 10 populations
####

## FIGURE S5

capt_history.bd_load %>% filter(pop_spec %in% 
    unique(data.all$pop_spec)[c(4, 5, 12, 13, 14, 15, 16, 17, 18, 19)]) %>%
  filter(!is.na(merc)) %>% droplevels() %>% {
  ggplot(., aes(merc, bd_load)) + 
    geom_point(aes(shape = Sex), alpha = 0.5) +
    geom_line(aes(group = Mark), alpha = 0.3) +
    scale_y_log10() + 
    scale_x_log10() +
    facet_wrap(~pop_spec, ncol = 3) +
    xlab("Methyl Mercury") +
    ylab("Bd Load") +
    theme(
      strip.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
    , axis.text.x = element_text(size = 12)
    )
}

####
## Some temp exploration
####

temp_data.all %>% 
  filter(
    Year != 2017
  , Site == "ScotiaBarrens") %>% {
  ggplot(., aes(yday, tmax)) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    geom_line() +
    facet_grid(~Year) +
    theme(
      axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
      ) +
    geom_hline(yintercept = 17, linetype = "dashed", colour = "firebrick3") +
    geom_hline(yintercept = 25, linetype = "dashed", colour = "firebrick3") +
    xlab("Julian Day") +
    ylab("Daily average temperature") 
}

capt_history.p %>% filter(swabbed == 1) %>%
  filter(Year != 2017) %>%
  filter(Site == "ScotiaBarrens") %>% {
  ggplot(., aes(days_in_opt_u_s, log_bd_load)) +
    geom_line(aes(group = Mark), alpha = 0.3) +
    geom_point(aes(colour = as.factor(Year)), alpha = 0.8) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    theme(
      axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
      ) +
    xlab("Days within Bd's thermal optima (17-25)") +
    ylab("Log Bd Load") 
  }



