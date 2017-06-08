
# Northern Ireland Mortality exploration 

# Jon Minton
# 7/6/2017


# To dos ------------------------------------------------------------------

# 1. Create function factory for leveplot defauls

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
  mutate(P = year - min(year)) %>% 
  mutate(
    improvement_phase = ifelse(year < 1955, T, F),
    yrs_since_start_improvement = year - 1922
    ) %>% 
  mutate(
    ww2 = ifelse(year >= 1939 & year <= 1945, T, F) 
  ) %>% 
  mutate(
    age_since_35 = age - 35,
    age_since_35 = ifelse(age_since_35 < 0, 0, age_since_35)
  ) %>% 
  mutate(
    is_adult = age >= 18
  ) %>% 
  mutate(
    flg_72_77 = year >= 1972 & year <= 1977, # 1972-1977
    flg_78_83 = year >= 1978 & year <= 1983,# 1978-1983
    flg_84_89 = year >= 1984 & year <= 1989,# 1984-1989
    flg_90_95 = year >= 1990 & year <= 1995,# 1990-1995
    flg_96_00 = year >= 1996 & year <= 2000# 1996-2000
  ) %>% 
  mutate(
    trbls_lim = ifelse(year < 1972, 0, 1),
    trbls_0_05pc = ifelse(year < 1972, 0, (1 - 0.0005) ^ (year - 1972)),
    trbls_0_10pc = ifelse(year < 1972, 0, (1 - 0.0010) ^ (year - 1972)),
    trbls_0_25pc = ifelse(year < 1972, 0, (1 - 0.0025) ^ (year - 1972)),
    trbls_0_5pc = ifelse(year < 1972, 0, (1 - 0.005) ^ (year - 1972)),
    trbls_1pc = ifelse(year < 1972, 0, (1 - 0.01) ^ (year - 1972)),
    trbls_1_5pc = ifelse(year < 1972, 0, (1 - 0.015) ^ (year - 1972)),
    trbls_2pc = ifelse(year < 1972, 0, (1 - 0.02) ^ (year - 1972)),
    trbls_5pc = ifelse(year < 1972, 0, (1 - 0.05) ^ (year - 1972)),
    trbls_10pc = ifelse(year < 1972, 0, (1 - 0.10) ^ (year - 1972))
  )  -> dta_male_wth_flags

dta_nir %>% 
  filter(age >=15 & age <= 45) %>% 
  filter(sex == "female") %>% 
  mutate(P = year - min(year)) %>% 
  mutate(
    improvement_phase = ifelse(year < 1955, T, F),
    yrs_since_start_improvement = year - 1922
  ) %>% 
  mutate(
    ww2 = ifelse(year >= 1939 & year <= 1945, T, F) 
  ) %>% 
  mutate(
    age_since_35 = age - 35,
    age_since_35 = ifelse(age_since_35 < 0, 0, age_since_35)
  ) %>% 
  mutate(
    is_adult = age >= 18
  ) %>% 
  mutate(
    flg_72_77 = year >= 1972 & year <= 1977, # 1972-1977
    flg_78_83 = year >= 1978 & year <= 1983,# 1978-1983
    flg_84_89 = year >= 1984 & year <= 1989,# 1984-1989
    flg_90_95 = year >= 1990 & year <= 1995,# 1990-1995
    flg_96_00 = year >= 1996 & year <= 2000# 1996-2000
  ) %>% 
  mutate(
    trbls_lim = ifelse(year < 1972, 0, 1),
    trbls_0_05pc = ifelse(year < 1972, 0, (1 - 0.0005) ^ (year - 1972)),
    trbls_0_10pc = ifelse(year < 1972, 0, (1 - 0.0010) ^ (year - 1972)),
    trbls_0_25pc = ifelse(year < 1972, 0, (1 - 0.0025) ^ (year - 1972)),
    trbls_0_5pc = ifelse(year < 1972, 0, (1 - 0.005) ^ (year - 1972)),
    trbls_1pc = ifelse(year < 1972, 0, (1 - 0.01) ^ (year - 1972)),
    trbls_1_5pc = ifelse(year < 1972, 0, (1 - 0.015) ^ (year - 1972)),
    trbls_2pc = ifelse(year < 1972, 0, (1 - 0.02) ^ (year - 1972)),
    trbls_5pc = ifelse(year < 1972, 0, (1 - 0.05) ^ (year - 1972)),
    trbls_10pc = ifelse(year < 1972, 0, (1 - 0.10) ^ (year - 1972))
  ) -> dta_female_wth_flags

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

# Model age 35+ separately
glm.nb(
  deaths ~ age * year + age_since_35 + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_ay_35plus

glm.nb(
  deaths ~ age + year + age_since_35 + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus

AIC(mdl_a_y, mdl_ay, mdl_a_y_35plus, mdl_ay_35plus)

# Adulthood effect
glm.nb(
  deaths ~ age + year + is_adult + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_adlt

glm.nb(
  deaths ~ age * year + is_adult + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_ay_adlt

glm.nb(
  deaths ~ age * year + is_adult + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_ay_adlt

glm.nb(
  deaths ~ age * year + age_since_35 + is_adult + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_ay_35plus_adlt

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt

AIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt)

BIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt)

# BIC is more conservative, and suggests against adding interaction effects 

# Best BIC so far is 
# mdl_a_y_35plus_adlt

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2

BIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt,
    mdl_a_y_35plus_adlt_ww2
    )

# Improvement phase
glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv

AIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt,
    mdl_a_y_35plus_adlt_ww2, mdl_a_y_35plus_adlt_ww2_imrv
)


BIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt,
    mdl_a_y_35plus_adlt_ww2, mdl_a_y_35plus_adlt_ww2_imrv
)

# Suggests a further improvement now using both AIC and BIC

# Now to add flags for 
# 1972-1977
# 1978-1983
# 1984-1989
# 1990-1995
# 1996-2000

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + flg_72_77 + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_flg01

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + flg_72_77 + flg_78_83 + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_flg02

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + flg_72_77 + flg_78_83 + flg_84_89 + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_flg03

AIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt,
    mdl_a_y_35plus_adlt_ww2, mdl_a_y_35plus_adlt_ww2_imrv,
    mdl_a_y_35plus_adlt_ww2_imrv_flg01,
    mdl_a_y_35plus_adlt_ww2_imrv_flg02,
    mdl_a_y_35plus_adlt_ww2_imrv_flg03
)


BIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt,
    mdl_a_y_35plus_adlt_ww2, mdl_a_y_35plus_adlt_ww2_imrv,
    mdl_a_y_35plus_adlt_ww2_imrv_flg01,
    mdl_a_y_35plus_adlt_ww2_imrv_flg02,
    mdl_a_y_35plus_adlt_ww2_imrv_flg03
)


summary(mdl_a_y_35plus_adlt_ww2_imrv_flg01)


# Now to look at the same model spec for females

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + flg_72_77 + offset(log(exposure)),
  data = dta_female_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_flg01_female

summary(mdl_a_y_35plus_adlt_ww2_imrv_flg01_female)


# Now to build some models that assume either a 1%, 2%, 5% or 10% decay in 
#effect after 1972

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_lim + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trblslim

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_0_05pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trbls0_05pc

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_0_10pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trbls0_10pc

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_0_25pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trbls0_25pc

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_0_5pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trbls0_5pc

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_1pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trbls1pc

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_1_5pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trbls1_5pc

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_2pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trbls2pc

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_5pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trbls5pc

glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_10pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trbls10pc

AIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt,
    mdl_a_y_35plus_adlt_ww2, mdl_a_y_35plus_adlt_ww2_imrv,
    mdl_a_y_35plus_adlt_ww2_imrv_flg01,
    mdl_a_y_35plus_adlt_ww2_imrv_flg02,
    mdl_a_y_35plus_adlt_ww2_imrv_flg03,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_05pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_10pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_25pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls1pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls1_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls2pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls10pc
)


BIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt,
    mdl_a_y_35plus_adlt_ww2, mdl_a_y_35plus_adlt_ww2_imrv,
    mdl_a_y_35plus_adlt_ww2_imrv_flg01,
    mdl_a_y_35plus_adlt_ww2_imrv_flg02,
    mdl_a_y_35plus_adlt_ww2_imrv_flg03,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_05pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_10pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_25pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls1pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls1_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls2pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls10pc
)

# Both suggest 1% has better fit than alternatives 

# 0.1% better still.

# This suggests to me the baseline model hasn't accounted enough for
# the collapse in improvment in male mortality after the improvement phase. 


# What if I re-do without a separate year coefficient


glm.nb(
  deaths ~ age + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_lim + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_35plus_adlt_ww2_imrv_trblslim

glm.nb(
  deaths ~ age + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_0_05pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_35plus_adlt_ww2_imrv_trbls0_05pc

glm.nb(
  deaths ~ age + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_0_10pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_35plus_adlt_ww2_imrv_trbls0_10pc

glm.nb(
  deaths ~ age + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_0_25pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_35plus_adlt_ww2_imrv_trbls0_25pc

glm.nb(
  deaths ~ age + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_0_5pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_35plus_adlt_ww2_imrv_trbls0_5pc

glm.nb(
  deaths ~ age + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_1pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_35plus_adlt_ww2_imrv_trbls1pc

glm.nb(
  deaths ~ age + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_1_5pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_35plus_adlt_ww2_imrv_trbls1_5pc

glm.nb(
  deaths ~ age + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_2pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_35plus_adlt_ww2_imrv_trbls2pc

glm.nb(
  deaths ~ age + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_5pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_35plus_adlt_ww2_imrv_trbls5pc

glm.nb(
  deaths ~ age + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_10pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_35plus_adlt_ww2_imrv_trbls10pc

AIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt,
    mdl_a_y_35plus_adlt_ww2, mdl_a_y_35plus_adlt_ww2_imrv,
    mdl_a_y_35plus_adlt_ww2_imrv_flg01,
    mdl_a_y_35plus_adlt_ww2_imrv_flg02,
    mdl_a_y_35plus_adlt_ww2_imrv_flg03,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_05pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_10pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_25pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls1pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls1_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls2pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls10pc,
    mdl_a_35plus_adlt_ww2_imrv_trblslim,
    mdl_a_35plus_adlt_ww2_imrv_trbls0_05pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls0_10pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls0_25pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls0_5pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls1pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls1_5pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls2pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls5pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls10pc
    
)


BIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt,
    mdl_a_y_35plus_adlt_ww2, mdl_a_y_35plus_adlt_ww2_imrv,
    mdl_a_y_35plus_adlt_ww2_imrv_flg01,
    mdl_a_y_35plus_adlt_ww2_imrv_flg02,
    mdl_a_y_35plus_adlt_ww2_imrv_flg03,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_05pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_10pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_25pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls1pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls1_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls2pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls10pc,
    mdl_a_35plus_adlt_ww2_imrv_trblslim,
    mdl_a_35plus_adlt_ww2_imrv_trbls0_05pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls0_10pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls0_25pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls0_5pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls1pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls1_5pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls2pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls5pc,
    mdl_a_35plus_adlt_ww2_imrv_trbls10pc
    
)


# No, removing year makes it much worse. 

# Alternative assumption is that there is both a step
# and also an exponential decay component to the additional mort

# Step + 0.5%
glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_lim + trbls_0_5pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_0_5pc

# Step + 1%
glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_lim + trbls_1pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_1pc

# Step + 2%
glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_lim + trbls_2pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_2pc

# Step + 5%
glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_lim + trbls_5pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_5pc

# Step + 10%
glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_lim + trbls_10pc + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_10pc

# What about the initial flag?

# Step + init flag
glm.nb(
  deaths ~ age + year + age_since_35 + is_adult + ww2 + I(yrs_since_start_improvement * improvement_phase) + 
    + trbls_lim + flg_72_77 + offset(log(exposure)),
  data = dta_male_wth_flags
) -> mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_flg01

AIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt,
    mdl_a_y_35plus_adlt_ww2, mdl_a_y_35plus_adlt_ww2_imrv,
    mdl_a_y_35plus_adlt_ww2_imrv_flg01,
    mdl_a_y_35plus_adlt_ww2_imrv_flg02,
    mdl_a_y_35plus_adlt_ww2_imrv_flg03,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_05pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_10pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_25pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls1pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls1_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls2pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls10pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_0_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_1pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_2pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_10pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_flg01
    
)


BIC(mdl_a_y, mdl_ay, mdl_a_y_adlt, mdl_ay_adlt, 
    mdl_a_y_35plus, mdl_ay_35plus, mdl_a_y_35plus_adlt, mdl_ay_35plus_adlt,
    mdl_a_y_35plus_adlt_ww2, mdl_a_y_35plus_adlt_ww2_imrv,
    mdl_a_y_35plus_adlt_ww2_imrv_flg01,
    mdl_a_y_35plus_adlt_ww2_imrv_flg02,
    mdl_a_y_35plus_adlt_ww2_imrv_flg03,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_05pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_10pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_25pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls0_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls1pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls1_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls2pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trbls10pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_0_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_1pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_2pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_5pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_10pc,
    mdl_a_y_35plus_adlt_ww2_imrv_trblslim_plus_flg01
)



