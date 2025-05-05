
c(
  which(capt_history.phi$phi_zeros == 1)
, which(capt_history.phi$phi_ones == 1 & capt_history.phi$phi_zeros == 0)
, which(capt_history.phi$offseason == 0 & capt_history.phi$phi_ones == 0 & capt_history.phi$phi_zeros == 0)
, which(capt_history.phi$offseason == 1 & capt_history.phi$phi_zeros == 0)
) %>% length()

which(
which(
  stan.fit.samples.l$phi %>% colMeans() == 1
)
  != 
which(
  stan.fit.samples.nl$phi %>% colMeans() == 1
)
)

data.frame(
  first_vals  = stan.fit.samples.l$phi %>% colMeans()
, second_vals = stan.fit.samples.nl$phi %>% colMeans()
) %>% filter(first_vals != 0, first_vals != 1) %>% {
  ggplot(., aes(first_vals, second_vals)) + geom_point() +
    geom_abline(intercept = 0, slope = 1)
}

data.frame(
  first_vals  = stan.fit.samples.l$X %>% colMeans()
, second_vals = stan.fit.samples.nl$X %>% colMeans()
) %>% filter(first_vals != 0, first_vals != 1) %>% {
  ggplot(., aes(first_vals, second_vals)) + geom_point() +
    geom_abline(intercept = 0, slope = 1)
}

data.frame(
  first_vals  = stan.fit.samples.l$chi %>% colMeans()
, second_vals = stan.fit.samples.nl$chi %>% colMeans()
) %>% filter(first_vals != 0, first_vals != 1) %>% {
  ggplot(., aes(first_vals, second_vals)) + geom_point() +
    geom_abline(intercept = 0, slope = 1)
}

data.frame(
  first_vals  = stan.fit.samples.nl$beta_offseason_sex %>% colMeans()
, second_vals = stan.fit.samples.l$beta_offseason_sex %>% colMeans()
) %>% filter(first_vals != 0, first_vals != 1) %>% {
  ggplot(., aes(first_vals, second_vals)) + geom_point() +
    geom_abline(intercept = 0, slope = 1)
}

data.frame(
  first_vals  = stan.fit.samples.nl$beta_offseason
, second_vals = stan.fit.samples.l$beta_offseason
) %>% reshape2::melt() %>% {
  ggplot(., aes(x = value)) + geom_histogram(bins = 50) +
    facet_wrap(~variable, ncol = 1)
}


