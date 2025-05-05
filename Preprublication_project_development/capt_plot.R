print(capt_history %>% {
  ggplot(., aes(SecNumConsec, Mark, fill = as.factor(captured))) + 
    geom_tile() +
    geom_point(data = 
        capt_history %>% 
        filter(swabbed == 1)
      , aes(x = SecNumConsec, y = Mark, z = NULL), lwd = 0.7) +
    xlab("Sampling Event") +
    ylab("Individual") +
    scale_fill_manual(
        values = c("dodgerblue4", "firebrick4")
      , name   = "Detected?"
      , labels = c("No", "Yes")) +
    scale_x_continuous(breaks = seq(1, max(capt_history$SecNumConsec), by = 2)) +
    guides(alpha = FALSE) +
    theme(
      axis.text.y = element_text(size = 6)
    , legend.text = element_text(size = 12)
    , legend.key.size = unit(.55, "cm")
    ) 
}
)

