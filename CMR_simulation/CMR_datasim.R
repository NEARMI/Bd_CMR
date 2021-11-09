####
## Run the sim to create the data
####

for (pop_ind in 1:n_pop) {

expdat  <- bd.simulate(
  periods   = periods[pop_ind, ]
, times     = times[pop_ind, ]
, all_ind   = all_ind[pop_ind, ]
, bd_beta   = bd_beta[pop_ind, ]
, bd_sigma  = bd_sigma[pop_ind, ]
, bd_theta  = bd_theta[pop_ind, ]
, obs_noise = obs_noise[pop_ind, ]
, bd_noinf  = bd_noinf[pop_ind, ]
)

# expdat %>% filter(periods == 1) %>% {ggplot(., aes(times, log_bd_load)) + geom_line(aes(group = ind))}    

one_pop <- bd.sampling(
  expdat    = expdat
, all_ind   = all_ind[pop_ind, ]
, new_ind   = new_ind[[pop_ind]]
, times     = times[pop_ind, ]
, periods   = periods[pop_ind, ]
, when_samp = when_samp[pop_ind, ]
, samp      = samp[[pop_ind]]
, bd_perc   = bd_perc[[pop_ind]]
, inbetween = inbetween[[pop_ind]]
, between_season_duration = between_season_duration[pop_ind, ]
, bd_mort   = bd_mort[pop_ind, ]
, bd_detect = bd_detect[pop_ind, ]
, p_mort    = p_mort[pop_ind, ]
, background_mort = background_mort[pop_ind, ]
, pop_ind   = pop_ind
)

one_pop.long <- bd.stan_org(
  one_pop  = one_pop
, pop_ind  = pop_ind
, times    = times[pop_ind, ]
, periods  = periods[pop_ind, ]
, samp     = samp[[pop_ind]]
)

print(paste("Population", pop_ind, "Simulated", sep = " "))

### Put the pops together
if (pop_ind == 1) {
ind_occ_phi.all   <- one_pop.long$ind_occ_phi
ind_occ_p.all     <- one_pop.long$ind_occ_p
X_bd.m.all        <- one_pop.long$X_bd.m
capture_range.all <- one_pop.long$one_pop$capture_range
ind_occ_size.all  <- one_pop.long$one_pop$ind_occ_size
all_ind.all       <- one_pop.long$one_pop$all_ind
each_ind.all      <- one_pop.long$one_pop$all_ind
ind_occ.all       <- one_pop.long$one_pop$all_ind * sum(samp[[pop_ind]])
ind_occ_min1.all  <- one_pop.long$one_pop$all_ind * (sum(samp[[pop_ind]]) - 1)
pop_cov.bd.all    <- one_pop.long$pop_cov.bd
ind_in_pop.all    <- rep(pop_ind, length(unique(one_pop.long$ind_occ_p$ind)))
expdat.all        <- one_pop$expdat
} else  {
ind_occ_phi.all   <- rbind(ind_occ_phi.all, one_pop.long$ind_occ_phi)
ind_occ_p.all     <- rbind(ind_occ_p.all, one_pop.long$ind_occ_p)
X_bd.m.all        <- rbind(X_bd.m.all, one_pop.long$X_bd.m)
capture_range.all <- rbind(capture_range.all, one_pop.long$one_pop$capture_range)
ind_occ_size.all  <- c(ind_occ_size.all, one_pop.long$one_pop$ind_occ_size)
all_ind.all       <- all_ind.all + one_pop.long$one_pop$all_ind
each_ind.all      <- c(each_ind.all, one_pop.long$one_pop$all_ind)
ind_occ.all       <- ind_occ.all + one_pop.long$one_pop$all_ind * sum(samp[[pop_ind]])
ind_occ_min1.all  <- ind_occ_min1.all + one_pop.long$one_pop$all_ind * (sum(samp[[pop_ind]]) - 1)
pop_cov.bd.all    <- rbind(pop_cov.bd.all, one_pop.long$pop_cov.bd)
ind_in_pop.all    <- c(ind_in_pop.all, rep(pop_ind, length(unique(one_pop.long$ind_occ_p$ind))))
expdat.all        <- rbind(expdat.all, one_pop$expdat)
}

}
