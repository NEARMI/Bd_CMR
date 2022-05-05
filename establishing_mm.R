##############################################################
## Building the complete model matrices for the joint model ##
##############################################################

## NOTE: This is ``non-dynamic'' in the sense that it needs to be changed manually to match
 ## whatever the linear predictors look like in the stan model

## See bottom of script for current linear predictors (again, non-dynamic, code must be updated
 ## whenever these linear predictors change)

####
## Model matrices for the survival part of the model
####

## create the model matrices, separated for intercepts and slopes (which all use the same model matrix)
if (n_spec > 1) {
fe_mm_phi_int   <- model.matrix(~Species + Sex, capt_history.phi)[phi_off_index, ]
fe_mm_phi_slope <- model.matrix(~Species*log_bd_load, capt_history.phi)[phi_off_index, ((n_spec + 1):(n_spec*2))]
fe_mm_phi_slope <- ifelse(fe_mm_phi_slope != 0, 1, 0)  
} else {
fe_mm_phi_int   <- model.matrix(~Sex, capt_history.phi)[phi_off_index, ]
}

re_mm_phi <- model.matrix(~-1+pop_spec, capt_history.phi)[phi_off_index, ]
re_mm_phi <- ifelse(re_mm_phi != 0, 1, 0)

####
## Model matrices for the detection part of the model
####

## create the model matrices, separated for intercepts and slopes (which all use the same model matrix)
if (n_spec > 1) {
fe_mm_p_int   <- model.matrix(~Species + Sex, capt_history.p)[p_est_index, ]
} else {
fe_mm_p_int   <- model.matrix(~Sex, capt_history.p)[p_est_index, ]
}
fe_mm_p_slope <- model.matrix(~-1+drawdown_cont + veg_cont, capt_history.p)[p_est_index, ]

re_mm_p <- model.matrix(~-1+pop_spec+date_fac, capt_history.p)[p_est_index, ]
re_mm_p <- model.matrix(~-1+date_fac, capt_history.p)[p_est_index, ]
re_mm_p <- ifelse(re_mm_p != 0, 1, 0)

####
## Model matrices and other vectors for getting population size estimates
####

## Model matrix used to extract each unique species -by- sex intercept value
fe_mm_p_int.uni <- fe_mm_p_int %>% as.data.frame() %>% distinct()

## *** Need to come back through and make this part dynamic
if (n_spec > 1) {
colnames(fe_mm_p_int.uni) <- c(unique(capt_history.phi$Species) %>% as.character(), c("F", "U"))
} else {
colnames(fe_mm_p_int.uni) <- c("M", "F", "U")
}
spec_to_int <- matrix(
  data = seq(nrow(fe_mm_p_int.uni))
, nrow = n_spec
, ncol = n_sex
, byrow = T
)

## number of each sex captured each day 
n_capt_per_day_sex <- capt_history.p %>% group_by(date_fac, Sex) %>% summarize(num_capt = sum(captured)) %>%
    pivot_wider(., date_fac, values_from = num_capt, names_from = Sex) %>% ungroup() %>% dplyr::select(-date_fac) %>%
    as.matrix()

## Some pops dont have "U" for sex, but setting to zero won't affect their estimates
n_capt_per_day_sex[is.na(n_capt_per_day_sex)] <- 0

## the values of each of the continuous covariates used to estimate p on each sampling day (needed to get the
 ## daily detection probability which is needed to get that days population size)
fe_mm_p_slope_uni <- capt_history.p %>% ungroup() %>% group_by(date_fac) %>% slice(1) %>% ungroup() %>% dplyr::select(
  drawdown_cont, veg_cont
) %>% as.data.frame()

## the species being sampled on each day
spec_pop_se <- (capt_history.p %>% ungroup() %>% group_by(date_fac) %>% slice(1) %>% ungroup() %>% dplyr::select(Species) %>% as.data.frame() %>%
  mutate(Species = as.numeric(Species)))$Species


## The linear predictor strategy for phi
 ## int
# fe_mm.int * alpha + re_mm * offseason_pop 
 ## merc
# (fe_mm.slope * beta_merc) .* merc
 ## bd
# (fe_mm.slope * beta_bd + re.mm * offseason_pop_bd) .* bd
 ## len
# (fe_mm.slope * beta_len + re.mm * offseason_pop_len) .* len

## The linear predictor strategy for p
 ## int
# fe_mm.int * alpha_p + re_mm_pop * pop + re_mm_day * pop_day
 ## drawdown and veg
# fe_mm.slope * beta_p


