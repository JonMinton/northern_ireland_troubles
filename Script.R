
# Northern Ireland Mortality exploration 

# Jon Minton
# 7/6/2017


# Load packages -----------------------------------------------------------


rm(list = ls())

pacman::p_load(
  tidyverse,
  purrr,
  lattice,   latticeExtra,
  spatstat
)



# Harvest HMD (once) ------------------------------------------------------

# source("scripts/reharvest_hmd.R")


