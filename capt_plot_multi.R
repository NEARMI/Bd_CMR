print(capt_history %>%
    filter(pop_spec == which.dataset) %>% {
  ggplot(., aes(SecNumConsec, Mark, fill = as.factor(captured))) + 
    geom_tile() +
    geom_point(data = 
        capt_history %>% 
        filter(swabbed == 1, pop_spec == which.dataset)
      , aes(x = SecNumConsec, y = Mark, z = NULL), lwd = 0.7) +
    xlab("Sampling Event") +
    ylab("Individual") +
    scale_fill_manual(
        values = c("dodgerblue4", "firebrick4")
      , name   = "Detected?"
      , labels = c("No", "Yes")) +
    scale_x_continuous(breaks = seq(1, max(capt_history$SecNumConsec), by = 2)) +
    guides(alpha = FALSE) +
    geom_vline(
      data = capt_history.phi %>% filter(pop_spec == which.dataset) %>% filter(Mark == unique(Mark)[1]) %>% 
        dplyr::select(offseason, SecNumConsec) %>% filter(offseason == 1) %>% 
        mutate(SecNumConsec = SecNumConsec + 0.5)
      , aes(xintercept = SecNumConsec)) +
    geom_vline(
      data = capt_history.phi %>% filter(pop_spec == which.dataset) %>% 
        filter(Mark == unique(Mark)[1]) %>% 
        dplyr::select(time_gaps, SecNumConsec) %>% filter(time_gaps == 1) %>% 
        mutate(SecNumConsec = SecNumConsec + 0.5)
      , aes(xintercept = SecNumConsec), linetype = "dashed") +
    theme(
      axis.text.y = element_text(size = 6)
    , legend.text = element_text(size = 12)
    , legend.key.size = unit(.55, "cm")
    ) 
}
)

capt_history %>% 
  filter(pop_spec == which.dataset) %>%
  group_by(Mark) %>%
  summarize(
    num_capt = sum(captured)
  , num_swab = sum(swabbed)
  ) %>% reshape2::melt("Mark") %>% {
      ggplot(., aes(x = value)) +
      geom_histogram(bins = 10) +
      facet_wrap(~variable) +
      xlab("Number of Captures")
      }


