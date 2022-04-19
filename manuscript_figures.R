######################################################
## A place to store code for all manuscript figures ##
######################################################

#### Supplemental Newt Quadratic Figure ----

Bd_Newts_AllSites   <- read.csv("data/cleaned_cmr_csv/SB_NOVI.csv")

## Stupid dates in R
if (length(grep("/", Bd_Newts_AllSites$CaptureDate[1])) > 0) {

date_convert <- apply(matrix(Bd_Newts_AllSites$CaptureDate), 1, FUN = function (x) {
  a <- strsplit(x, split = "/")[[1]]
  b <- a[3]
  b <- strsplit(b, "")[[1]][c(3, 4)] %>% paste(collapse = "")
  paste(c(a[c(1, 2)], b), collapse = "/")
})

Bd_Newts_AllSites$CaptureDate <- date_convert
Bd_Newts_AllSites      %<>% mutate(CaptureDate = as.Date(CaptureDate, "%m/%d/%y"))

} else {
  
Bd_Newts_AllSites      %<>% mutate(CaptureDate = as.Date(CaptureDate))
  
}

Bd_Newts_AllSites$Year   <- apply(matrix(as.character(Bd_Newts_AllSites$CaptureDate)), 1, FUN = function (x) strsplit(x, "-")[[1]][1]) %>% as.numeric()
Bd_Newts_AllSites$Month  <- apply(matrix(as.character(Bd_Newts_AllSites$CaptureDate)), 1, FUN = function (x) strsplit(x, "-")[[1]][2]) %>% as.numeric()
Bd_Newts_AllSites$Julian <- as.POSIXlt(Bd_Newts_AllSites$CaptureDate)$yday

ggplot(
  Bd_Newts_AllSites %>% filter(BdSample == "Y") %>% droplevels()
  , aes(Julian, TargetCopies.swab)) + 
  geom_line(aes(group = IndividualID, colour = as.factor(Year))) +
  facet_wrap(~Year) +
  scale_y_log10() +
  scale_colour_brewer(palette = "Dark2", name = "Year") +
  xlab("Julian Day") +
  ylab("Bd Load (Copies/Swab)") + 
  theme(legend.key.size = unit(0.8, "cm"))

#### Checking gaps between sampling occasions in all populations ----

## Note: Run after data_manip.R

expand.grid(
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
    scale_color_brewer(palette = "Dark2", name = "Year") +
    geom_hline(yintercept = 2.5, linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = 7.5, linetype = "dashed", size = 0.5) +
    geom_hline(yintercept = 9.5, linetype = "dashed", size = 0.5) + 
    geom_hline(yintercept = 14.5, linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = 82, size = 0.75, colour = "dodgerblue4", linetype = "dashed") +
    geom_vline(xintercept = 135, size = 0.75, colour = "firebrick4") +
    xlab("Days Bewtween Sampling Occasions") +
    ylab("Population")
}

expand.grid(
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
    ylab("Population")
}


#### Exploring sampling schemes in specific populations ----

samp_days <- capt_history %>% group_by(pop_spec, capture_date) %>% slice(1) %>% 
  mutate(J_date = as.POSIXlt(capture_date)$yday) %>%
  dplyr::select(pop_spec, sampled, capture_date, Year, Month, J_date)
all_days  <- expand.grid(
  pop_spec = unique(samp_days$pop_spec)
, Year     = unique(samp_days$Year)
, J_date   = seq(365)
)
samp_days <- left_join(all_days, samp_days) %>%
  mutate(
    sampled    = ifelse(is.na(sampled), 0, sampled)
  , J_date.deg = (360 / 365) * J_date
  , y_coord    = 1) %>% filter(sampled == 1) %>%
  mutate(y_coord = y_coord * as.numeric(as.factor(Year)))


between_season_segment <- capt_history.phi %>% filter(offseason == 1) %>% 
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
    y_coord = (Year - min(samp_days$Year)) + 1
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

within_season_segment <- capt_history.phi %>% filter(phi_ones == 0 & offseason == 0) %>% 
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
    y_coord = (Year - min(samp_days$Year)) + 1
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

closed_season_segment <- capt_history.phi %>% filter(phi_ones == 1) %>% 
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
    y_coord = (Year - min(samp_days$Year)) + 1
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

which_pops_to_plot <- unique(data.all$pop_spec)[c(15:21)]

samp_days %>% filter(pop_spec %in% which_pops_to_plot) %>% {
  ggplot(., aes(J_date.deg, y_coord)) + 
    geom_point(size = 1.5, alpha = 0.5) + 
    coord_polar() + 
    xlab("") + 
    ylab("") + 
    scale_y_continuous(
      lim    = c(0.5, 5.1)
    , breaks = unique((samp_days %>% filter(pop_spec %in% which_pops_to_plot))$Year) %>% length() %>% seq()
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
    facet_wrap(~pop_spec, nrow = 2) +
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
    theme(
      axis.text.x = element_text(size = 10)
    , axis.text.y = element_text(size = 13)
    , axis.ticks  = element_blank()
    , strip.text = element_text(size = 13)
    , plot.margin = unit(c(.1,.1,.1,.1), "cm")
    ) 
}


#### Bd over time and temperature ----

## Broadest view of Bd loads across temperatures
capt_history.p %>% filter(swabbed == 1) %>% {
  ggplot(., aes(cumtemp, log_bd_load)) +
    geom_point(aes(colour = as.factor(Year))) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    facet_wrap(~Site, scales = "free") +
    theme(
      axis.text.x = element_text(size = 10)
    , axis.text.y = element_text(size = 10)
    ) 
}

## Broken down by a species and year
capt_history.p %>% filter(swabbed == 1) %>%
  filter(Species == "RANA") %>% {
  ggplot(., aes(cumtemp, log_bd_load)) +
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

## Checking Rana loads vs dates and temps
RANAgg1 <- capt_history.p %>% filter(swabbed == 1) %>%
  filter(Species == "NOVI") %>% {
  ggplot(., aes(yday, log_bd_load)) +
    geom_point(aes(colour = as.factor(Year))) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    theme(
      axis.text.y = element_text(size = 12)
    , axis.text.x = element_blank()
    , axis.title.x = element_blank()
    ) +
    scale_x_continuous(lim = c(0, 365)) +
    facet_wrap(~Site, nrow = 1) +
    xlab("Julian Day") +
    ylab("Log Bd Load") +
  ggtitle("Eastern Newt")
  }

## and temperature of each of these populations 
RANAgg2 <- temp_data.all %>%
  filter(Site %in% 
#  c(
#    "DilmanMeadows", "FoxCreek", "JonesPond"
#  , "LostHorse", "SanFrancisquito", "SummitMeadow"
#  , "ThreeCreeks"
#  )
#  c(
#    "Blackrock-C", "Blackrock-H", "JonesPond"
#  , "SonomaMountain", "TwoMedicine"
#  )
  c(
    "EmmaCarlin", "MudLake", "ScotiaBarrens"
  , "SMNWR_W", "SPR"
  )
    ) %>% filter(Year != 2017) %>% {
  ggplot(., aes(yday, tmean)) +
    geom_line(aes(colour = as.factor(Year)), alpha = 0.5) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    theme(
      axis.text.x = element_text(size = 10)
    , axis.text.y = element_text(size = 10)
    , strip.text.x = element_blank()
    ) +
    facet_wrap(~Site, nrow = 1) +
    geom_hline(yintercept = 17, linetype = "dashed", colour = "firebrick3") +
    geom_hline(yintercept = 25, linetype = "dashed", colour = "firebrick3") +
    xlab("Julian Day") +
    ylab("Mean Daily Temperature") 
}

RANAgg3 <- capt_history.p %>% filter(swabbed == 1) %>%
  filter(Species == "RANA") %>% {
  ggplot(., aes(yday, cumtemp)) +
    geom_point(aes(colour = as.factor(Year))) +
    scale_colour_brewer(palette = "Dark2", name = "Year") +
    theme(
      axis.text.x = element_text(size = 12)
    , axis.text.y = element_text(size = 12)
    ) +
    facet_wrap(~Site, nrow = 1, scales = "free") +
    xlab("Julian Day") +
    ylab("Log Bd Load") +
  ggtitle("Rana Species")
  }

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
      ) +
    xlab("Julian Day") +
    ylab("Log Bd Load") 
  }

SBgg.3 <- capt_history.p %>% filter(swabbed == 1) %>%
  filter(Year != 2017) %>%
  filter(Site == "ScotiaBarrens") %>% {
  ggplot(., aes(days_in_opt, log_bd_load)) +
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
