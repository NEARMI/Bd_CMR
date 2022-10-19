###################################
## Needed packages and functions ##
###################################

needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan", "gridExtra")
lapply(needed_packages, require, character.only = TRUE)
#source("../ggplot_theme.R")
set.seed(10002)
'%notin%' <- Negate('%in%')

## Function to scale Bd from Ct to load

ct_to_load <- function(x) {
  ((10^((x-38.74861)/-3.366983))*(125/5))   ## Or maybe 140 and not 125 per Dan's code?
}

## Functions to scale MeHg

scale_newt <- function(x) {
 10^(0.2092 + 0.9882 * log10(x))
}

scale_frog <- function(x) {
 10^(0.3255 + 0.9511 * log10(x))
}
