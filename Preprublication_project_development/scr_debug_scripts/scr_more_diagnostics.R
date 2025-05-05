
out.pred.w.1 <- rbind(out.pred.bd, out.pred.s)
out.pred.b.1 <- rbind(out.pred.bd, out.pred.s)
out.pred.p.1 <- rbind(out.pred.bd, out.pred.s)
ind_order.1  <- ind_order

out.pred.w.2 <- rbind(out.pred.bd, out.pred.s)
out.pred.b.2 <- rbind(out.pred.bd, out.pred.s)
out.pred.p.2 <- rbind(out.pred.bd, out.pred.s)
ind_order.2  <- ind_order

out.pred.w.3 <- rbind(out.pred.bd, out.pred.s)
out.pred.b.3 <- rbind(out.pred.bd, out.pred.s)
out.pred.p.3 <- rbind(out.pred.bd, out.pred.s)
ind_order.3  <- ind_order

out.pred.w.1 %<>% mutate(pop = 1)
out.pred.b.1 %<>% mutate(pop = 1)
out.pred.p.1 %<>% mutate(pop = 1)
ind_order.1  %<>% mutate(pop = 1)

out.pred.w.2 %<>% mutate(pop = 2)
out.pred.b.2 %<>% mutate(pop = 2)
out.pred.p.2 %<>% mutate(pop = 2)
ind_order.2  %<>% mutate(pop = 2)

out.pred.w.3 %<>% mutate(pop = 3)
out.pred.b.3 %<>% mutate(pop = 3)
out.pred.p.3 %<>% mutate(pop = 3)
ind_order.3  %<>% mutate(pop = 3)

out.pred.w <- rbind(out.pred.w.1, out.pred.w.2, out.pred.w.3)
out.pred.b <- rbind(out.pred.b.1, out.pred.b.2, out.pred.b.3)
out.pred.p <- rbind(out.pred.p.1, out.pred.p.2, out.pred.p.3)
ind_order  <- rbind(ind_order.1, ind_order.2, ind_order.3)

out.pred.w %>%
  group_by(gap, variable, pop) %>%
  mutate(pop = as.factor(pop)) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% mutate(pop = plyr::mapvalues(pop
    , from = c("1", "2", "3")
    , to   = c(
      "Jones Lubrecht -- BUBO"
    , "Jones Lubrecht -- RALU"
    , "Lost Horse -- RALU"))
    , variable = plyr::mapvalues(variable
      , from = c("bd", "size")
      , to   = c("Bd Load (Copies/Swab)"
        , "Body Size (g) (scaled)"))) %>% {
  ggplot(., aes(gap, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = pop), alpha = 0.3) +
    geom_line(aes(colour = pop), lwd = 2) +
    facet_wrap(~variable, scales = "free") +
    theme(
      panel.spacing   = unit(2, "lines")
    , legend.key.size = unit(.55, "cm")
    , legend.title = element_text(size = 15)
    , legend.text  = element_text(size = 12)) +
    scale_fill_brewer(name = "Population", palette = "Dark2") +
    scale_colour_brewer(name = "Population", palette = "Dark2") +
    xlab("") +
    ylab("Probability animal remains in population")
  }

out.pred.b %>%
  group_by(gap, variable, pop) %>%
  mutate(pop = as.factor(pop)) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% mutate(pop = plyr::mapvalues(pop
    , from = c("1", "2", "3")
    , to   = c(
      "Jones Lubrecht -- BUBO"
    , "Jones Lubrecht -- RALU"
    , "Lost Horse -- RALU"))
    , variable = plyr::mapvalues(variable
      , from = c("bd", "size")
      , to   = c("Bd Load (Copies/Swab)"
        , "Body Size (g) (scaled)"))) %>% {
  ggplot(., aes(gap, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = pop), alpha = 0.3) +
    geom_line(aes(colour = pop), lwd = 2) +
    facet_wrap(~variable, scales = "free") +
    theme(
      panel.spacing   = unit(2, "lines")
    , legend.key.size = unit(.55, "cm")
    , legend.title = element_text(size = 15)
    , legend.text  = element_text(size = 12)) +
    scale_fill_brewer(name = "Population", palette = "Dark2") +
    scale_colour_brewer(name = "Population", palette = "Dark2") +
    xlab("") +
    ylab("Between Season Survival")
  }

out.pred.p %>%
  group_by(gap, variable, pop) %>%
  mutate(pop = as.factor(pop)) %>%
  summarize(
    lwr = quantile(value, 0.1)
  , mid = quantile(value, 0.5)
  , upr = quantile(value, 0.9)
  ) %>% mutate(pop = plyr::mapvalues(pop
    , from = c("1", "2", "3")
    , to   = c(
      "Jones Lubrecht -- BUBO"
    , "Jones Lubrecht -- RALU"
    , "Lost Horse -- RALU"))
    , variable = plyr::mapvalues(variable
      , from = c("bd", "size")
      , to   = c("Bd Load (Copies/Swab)"
        , "Body Size (g) (scaled)"))) %>% {
  ggplot(., aes(gap, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = pop), alpha = 0.3) +
    geom_line(aes(colour = pop), lwd = 2) +
    facet_wrap(~variable, scales = "free") +
    theme(
      panel.spacing   = unit(2, "lines")
    , legend.key.size = unit(.55, "cm")
    , legend.title = element_text(size = 15)
    , legend.text  = element_text(size = 12)) +
    scale_fill_brewer(name = "Population", palette = "Dark2") +
    scale_colour_brewer(name = "Population", palette = "Dark2") +
    xlab("") +
    ylab("Detection Probability")
        }

ind_order %>% mutate(pop = plyr::mapvalues(pop
    , from = c("1", "2", "3")
    , to   = c(
      "Jones Lubrecht -- BUBO"
    , "Jones Lubrecht -- RALU"
    , "Lost Horse -- RALU"))) %>% {
  ggplot(., aes(order_real, mid)) + 
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr)) +
    geom_point(aes(order_real, tot_bd), colour = "firebrick3") +
    xlab("Individual") + 
    ylab("Random Effect Deviate") +
    facet_wrap(~pop, scales = "free") +
    theme(axis.text.x = element_text(size = 8)) +
    geom_hline(yintercept = mean(ind_order$mid)) +
    geom_hline(yintercept = mean(ind_order$tot_bd, na.rm = T), colour = "firebrick3")
}

ind_order %>% mutate(pop = plyr::mapvalues(pop
    , from = c("1", "2", "3")
    , to   = c(
      "Jones Lubrecht -- BUBO"
    , "Jones Lubrecht -- RALU"
    , "Lost Horse -- RALU"))) %>% {
      ggplot(., aes(order_real, order_pred)) + 
  geom_point(size = 2) +
  facet_wrap(~pop, scales = "free") +
  scale_color_brewer(palette = "Dark2", name = "Individual type", labels = c("Middle", "Lowest Bd", "Highest Bd")) +
  xlab("Individual Ordered by Average of Measured Bd") +
  ylab("Predicted Order")
    }


