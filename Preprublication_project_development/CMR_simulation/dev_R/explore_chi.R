capt_history.p %>% 
   dplyr::select(Mark, Month, Year, time_gaps, p_zeros, pred_p, pred_chi) %>% 
   filter(Mark == 4) %>% as.data.frame()


chi_sub <- stan.fit.samples$chi[1, p_first_index[1]:(p_first_index[1] + rep(colSums(n_occ), n_ind.per)[1] - 1)]

p_sub   <- stan.fit.samples$p[1, p_first_index[1]:(rep(colSums(n_occ), n_ind.per)[1])]
phi_sub <- stan.fit.samples$phi[1, phi_first_index[1]:(rep(colSums(n_occ) - 1, n_ind.per)[1])]

test_occ <- 29
check_chi     <- numeric(test_occ)
check_chi[test_occ] <- 1
for (t in 1:(rep(colSums(n_occ), n_ind.per)[1] - 1)) {
  t_curr <- test_occ - t
  t_next <- t_curr + 1
  check_chi[t_curr] <- (1 - phi_sub[t_curr]) + phi_sub[t_curr] * (1 - p_sub[t_next - 1]) * chi_sub[t_next]	
}

check_chi == chi_sub
