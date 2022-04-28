####
## Building the complete model matrices for the joint model
####

## NOTE: This is ``non-dynamic'' in the sense that it needs to be changed manually to match
 ## whatever the linear predictors look like in the stan model

## The linear predictor strategy for phi
 ## int
# fe_mm.int * alpha + re_mm * offseason_pop 
 ## merc
# (fe_mm.slope * beta_merc) .* merc
 ## bd
# (fe_mm.slope * beta_bd + re.mm * offseason_pop_bd) .* bd
 ## len
# (fe_mm.slope * beta_len + re.mm * offseason_pop_len) .* len

## create the model matrices, separated for intercepts and slopes (which all use the same model matrix)
fe_mm_phi_int   <- model.matrix(~Species + Sex, capt_history.phi)[phi_off_index, ]
fe_mm_phi_slope <- model.matrix(~Species*log_bd_load, capt_history.phi)[phi_off_index, 4:6]
fe_mm_phi_slope <- ifelse(fe_mm_phi_slope != 0, 1, 0)  

re_mm_phi <- model.matrix(~-1+pop_spec, capt_history.phi)[phi_off_index, ]
re_mm_phi <- ifelse(re_mm_phi != 0, 1, 0)


## The linear predictor strategy for p
 ## int
# fe_mm.int * alpha_p + re_mm_pop * pop + re_mm_day * pop_day
 ## drawdown and veg
# fe_mm.slope * beta_p


## create the model matrices, separated for intercepts and slopes (which all use the same model matrix)
fe_mm_p_int   <- model.matrix(~Species + Sex, capt_history.p)[p_est_index, ]
fe_mm_p_slope <- model.matrix(~-1+drawdown_cont + veg_cont, capt_history.p)[p_est_index, ]

re_mm_p <- model.matrix(~-1+pop_spec+date_fac, capt_history.p)[p_est_index, ]
re_mm_p <- model.matrix(~-1+date_fac, capt_history.p)[p_est_index, ]
re_mm_p <- ifelse(re_mm_p != 0, 1, 0)

