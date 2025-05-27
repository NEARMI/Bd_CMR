###################################
## Needed packages and functions ##
###################################

needed_packages <- c(
  "magrittr", "dplyr", "tidyr", "readxl"
, "ggplot2", "rstan", "gridExtra", "RColorBrewer"
)

lapply(needed_packages, require, character.only = TRUE)

set.seed(10002)

## Some useful functions 

## For data wrangling
'%notin%'  <- Negate('%in%')

## Function to scale Bd from Ct to load
ct_to_load <- function(x) { ((10^((x-38.74861)/-3.366983))*(125/5)) }

## Functions to scale MeHg
scale_newt <- function(x) { 10^(0.2092 + 0.9882 * log10(x)) }
scale_frog <- function(x) { 10^(0.3255 + 0.9511 * log10(x)) }

## ggplot theme
theme_set(theme_bw())
suppressWarnings(
  theme_update(
    axis.text.x = element_text(size = 16)
    , axis.text.y = element_text(size = 16)
    , axis.title.x = element_text(size = 16)
    , axis.title.y = element_text(size = 16)
    , legend.title = element_text(size = 12)
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
    , strip.background = element_blank()
    , panel.margin = unit(0, "lines")
    , legend.key.size = unit(.55, "cm")
    , legend.key = element_rect(fill = "white")
    , panel.margin.y = unit(0.5, "lines")
    , panel.border = element_rect(colour = "black", fill = NA, size = 1)
    , strip.text.x = element_text(size = 16, colour = "black", face = "bold"))
)


