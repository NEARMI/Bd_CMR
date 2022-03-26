###################################
## Needed packages and functions ##
###################################

needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan", "gridExtra")
lapply(needed_packages, require, character.only = TRUE)
#source("../ggplot_theme.R")
set.seed(10002)
'%notin%' <- Negate('%in%')
