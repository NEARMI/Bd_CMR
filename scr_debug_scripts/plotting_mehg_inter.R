pred.vals <- expand.grid(
  bd   = seq(0, 14, by = 1)
, len  = seq(-3, 3, by = 0.5)
, mehg = seq(-3, 3, by = 0.5)
)

pred.est <- matrix(data = 0, nrow = nrow(pred.vals), ncol = dim(stan.fit.samples[[1]])[1])

for (j in 1:nrow(pred.est)) {
   pred.est[j, ] <- plogis(
    stan.fit.samples$beta_offseason[, 1] +
    stan.fit.samples$beta_offseason[, 2] * pred.vals[j, ]$bd   +
    stan.fit.samples$beta_offseason[, 3] * pred.vals[j, ]$len  +
    stan.fit.samples$beta_offseason[, 4] * pred.vals[j, ]$mehg +
    stan.fit.samples$beta_offseason[, 5] * pred.vals[j, ]$mehg * pred.vals[j, ]$bd +
    stan.fit.samples$beta_offseason_sex[, 1]
   )
}

pred.vals <- cbind(pred.vals, pred.est)

pred.vals %<>% pivot_longer(., c(-bd, -len, -mehg), names_to = "iter", values_to = "est")

pred.vals.gg <- pred.vals %>%
  group_by(bd, len, mehg) %>%
  summarize(
    lwr   = quantile(est, 0.025)
  , lwr_n = quantile(est, 0.200)
  , mid   = quantile(est, 0.500)
  , upr_n = quantile(est, 0.800)
  , upr   = quantile(est, 0.975)
  )

pred.vals.gg %>% filter(len == 0, mehg %in% c(-2, -1, 0, 1, 2)) %>%
  mutate(mehg = as.factor(mehg)) %>% {
  ggplot(., aes(bd, mid)) + 
    geom_ribbon(aes(ymin = lwr_n, ymax = upr_n, fill = mehg), alpha = 0.3) +
    geom_line(aes(colour = mehg), size = 1) + 
    scale_colour_discrete(name = "MeHg") +
    scale_fill_discrete(name = "MeHg") +
      facet_wrap(~mehg) +
    xlab("Bd Load") + ylab("Between Season Survival Probability") +
      ggtitle(this_pop)
}
