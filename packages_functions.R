###################################
## Needed packages and functions ##
###################################

needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan", "gridExtra")
lapply(needed_packages, require, character.only = TRUE)
#source("../ggplot_theme.R")
set.seed(10002)
'%notin%' <- Negate('%in%')

## Functions to scale MeHg

scale_newt <- function(x) {
 10^(0.2092 + 0.9882 * log10(x))
}

scale_frog <- function(x) {
 10^(0.3255 + 0.9511 * log10(x))
}
