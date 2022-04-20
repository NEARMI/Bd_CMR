## Checking phi, p, and chi indexing for individuals only captured once

 ## !!!! This will be an important debugging script to show in the final repo so dont delete

## Mark
mark_check <- 33

capt_history %>% 
  filter(Mark == mark_check) %>% 
  dplyr::select(-Month, -Year, -month_year, -pop_spec, -Region, -State, -Site) %>%
  mutate(
    pred_phi = c(stan.fit.samples$phi[, which(capt_history.phi$Mark == mark_check)] %>% colMeans(), 1)
  , pred_p   = stan.fit.samples$p[, which(capt_history.p$Mark == mark_check)] %>% colMeans()
  , pred_chi = stan.fit.samples$chi[, which(capt_history.p$Mark == mark_check)] %>% colMeans()
  ) %>%
  as.data.frame() 

## Loop from first to last observation
first_p_one <- (capture_range$first[mark_check] + 1)
last_e      <- capture_range$final[mark_check]

## entries that are 1 ~ bernoulli( )
entry_one <- phi_first_index[mark_check] - 1 + first_p_one - 1
entry_two <- phi_first_index[mark_check] - 1 + last_e - 1

stan.fit.samples$phi[1:5, which(capt_history.phi$Mark == mark_check)]
stan.fit.samples$phi[1:5, entry_one]
stan.fit.samples$phi[1:5, entry_two]

## entries for y ~ bernoulli
entry_one <- p_first_index[mark_check] - 1 + first_p_one
entry_two <- p_first_index[mark_check] - 1 + last_e

stan.fit.samples$p[1:5, which(capt_history.p$Mark == mark_check)]
stan.fit.samples$p[1:5, entry_one]
stan.fit.samples$p[1:5, entry_two]

## entries of chi
stan.fit.samples$chi[, which(capt_history.p$Mark == mark_check)] %>% colMeans()

p_first_index[mark_check]:(p_first_index[mark_check] + 25 - 1)


p_first_index[mark_check] - 1 + last_e

stan.fit.samples$chi[1:5, which(capt_history.p$Mark == mark_check)]
stan.fit.samples$chi[1:5, p_first_index[mark_check] - 1 + last_e]

