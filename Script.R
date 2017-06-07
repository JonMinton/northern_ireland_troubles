
# Northern Ireland Mortality exploration 

# Jon Minton
# 7/6/2017


# Load packages -----------------------------------------------------------


rm(list = ls())

pacman::p_load(
  tidyverse,
  purrr, rgl,
  RColorBrewer,
  lattice,   latticeExtra,
  spatstat,
  MASS
)



# Harvest HMD (once) ------------------------------------------------------

# source("scripts/reharvest_hmd.R")

# Note: Need to combine East and West Germany for some analyses 


# Produce levelplot with effective qualitative colourscheme ---------------------------

dta <- read_csv("data/tidy/new_counts.csv")

dta %>% 
  filter(country_code == "GBR_NIR") %>% 
  mutate(mr = deaths / exposure) %>% 
  mutate(lmr = log(mr, base = 10)) -> dta_nir

qual_col <- colorRampPalette(rev(brewer.pal(12, "Paired")))(200)


# Ages up to 90, linear scale
dta_nir %>% 
  filter(age <= 90) %>% 
  levelplot(
    mr ~ year * age | sex, data = .,
    col.regions = qual_col,
    at = seq(from = 0, to = 0.40, by = 0.002),
    scales=list(alternating=3),
    aspect= "iso"
  )

# Ages up to 90, log scale
dta_nir %>% 
  filter(age <= 90) %>% 
  levelplot(
    lmr ~ year * age | sex, data = .,
    col.regions = qual_col,
    at = seq(from = -5.0, to = 0.1, by = 0.1),
    scales=list(alternating=3),
    aspect= "iso"
  )


# Ages between 15 and 45, linear scale

dta_nir %>% 
  filter(age >= 15 & age <= 45) %>% 
  levelplot(
    mr ~ year * age | sex, data = .,
    col.regions = qual_col,
    at = seq(from = 0.000, to = 0.012, by = 0.00010),
    scales = list(alternating = 3),
    aspect= "iso"
  )

# Ages between 15 and 45, log scale
dta_nir %>% 
  filter(age >= 15 & age <= 45) %>% 
  levelplot(
    lmr ~ year * age | sex, data = .,
    col.regions = qual_col,
    at = seq(from = -4.5, to = -1.5, by = 0.05),
    scales = list(alternating = 3),
    aspect= "iso"
  )


# Structural model --------------------------------------------------------


# Improvements 1922-1955
# War 1939-1945 
# age 35+ ageing effect
# Adulthood onset effect (18+)
# post 1972 effect 

dta_nir %>% 
  filter(age >=15 & age <= 45) %>% 
  filter(sex == "male") %>% 
  mutate(
    improvement_phase = ifelse(year < 1955, T, F),
    yrs_since_start_improvement = year - 1922
    ) %>% 
  mutate(
    ww2 = ifelse(year >= 1939 & year <= 1945, T, F) 
  ) %>% 
  mutate(
    age_since_35 = max(0, age - 35)
  ) %>% 
  mutate(
    is_adult = age >= 18
  ) %>% 
  mutate(
    after_1972 = year >= 1972,
    yrs_since_1972 = max(0, year - 1972)
  ) -> dta_male_wth_flags

# No interaction
glm.nb(
  deaths ~ age + year + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y

# Interaction
glm.nb(
  deaths ~ age * year + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_ay

# Model w improvement phase only

glm.nb(
  deaths ~ age * I(yrs_since_start_improvement * improvement_phase) + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_imprvphase


AIC(mdl_a_y, mdl_ay, mdl_imprvphase)
# No need to model improvement phase separately, just yrs 


# Now add ww2 effect 

glm.nb(
  deaths ~ age * year + ww2 + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_ay_ww2

AIC(mdl_a_y, mdl_ay, mdl_imprvphase, mdl_ay_ww2)


