# Algeria data

# Jon Minton
# 27/7/2017


# To dos ------------------------------------------------------------------

# 1. Create function factory for leveplot defauls

# Load packages -----------------------------------------------------------


rm(list = ls())

pacman::p_load(
  MASS,
  readxl,
  tidyverse,
  purrr, rgl,
  RColorBrewer,
  lattice,   latticeExtra,
  spatstat
)


dta_male <- read_excel("data/algeria/Single_ages_mortality_surface_Algeria-V.-13.07.2017.xlsx", sheet = "Males")
dta_female <- read_excel("data/algeria/Single_ages_mortality_surface_Algeria-V.-13.07.2017.xlsx", sheet = "Females")

dta_male %>% 
  gather(key = "year", value = "risk", -qx) %>% 
  mutate(sex = "male") %>% 
  select(sex, age = qx, year, risk) -> dta_male

dta_female %>% 
  gather(key = "year", value = "risk", -qx) %>% 
  mutate(sex = "female") %>% 
  select(sex, age = qx, year, risk) -> dta_female

dta_tidy <- bind_rows(dta_male, dta_female) %>% 
  mutate(year = as.numeric(year))



qual_col <- colorRampPalette(rev(brewer.pal(12, "Paired")))(200)
rdbu_col <- colorRampPalette(rev(brewer.pal(5, "RdBu")))(200)


# Mortality risk
dta_tidy %>% 
  filter(age <= 90) %>% 
  levelplot(
    risk ~ year * age | sex, data = .,
    col.regions = qual_col,
    at = seq(from = 0, to = 0.40, by = 0.002),
    scales=list(alternating=3),
    aspect= "iso"
  )

# log10 mortality risk

dta_tidy %>% 
  filter(age <= 90) %>% 
  levelplot(
    log(risk, 10) ~ year * age | sex, data = .,
    col.regions = qual_col,
    at = seq(from = -0.5, to = -4.0, by = -0.05),
    scales=list(alternating=3),
    aspect= "iso"
  )



# 3d plot - identity scale - male

dta_tidy %>% 
  filter(age <= 90) %>% 
  filter(sex == "male") %>% 
  select(year, age, risk) %>% 
  spread(age, risk) -> tmp

yrs <- tmp$year
tmp$year <- NULL
ages <- names(tmp)
tmp <- as.matrix(tmp)

# Persp plot, males, all ages 
persp3d(z = tmp, col = "lightgray", 
        box = F, axes = F,
        aspect = c(1, 1, 91/93),
        xlab = "", ylab = "", zlab = ""
)



# 3d plot - log scale - men


dta_tidy %>% 
  filter(age <= 90) %>% 
  filter(sex == "male") %>% 
  mutate(risk = log(risk, 10)) %>% 
  select(year, age, risk) %>% 
  spread(age, risk) -> tmp

yrs <- tmp$year
tmp$year <- NULL
ages <- names(tmp)
tmp <- as.matrix(tmp)

# Persp plot, males, all ages 
persp3d(z = tmp, col = "lightgray", 
        box = F, axes = F,
        aspect = c(1, 1, 91/93),
        xlab = "", ylab = "", zlab = ""
)


dta_tidy %>% 
  filter(age <= 90) %>% 
  filter(sex == "female") %>% 
  mutate(risk = log(risk, 10)) %>% 
  select(year, age, risk) %>% 
  spread(age, risk) -> tmp

yrs <- tmp$year
tmp$year <- NULL
ages <- names(tmp)
tmp <- as.matrix(tmp)

# Persp plot, males, all ages 
persp3d(z = tmp, col = "lightgray", 
        box = F, axes = F,
        aspect = c(1, 1, 91/93),
        xlab = "", ylab = "", zlab = ""
)



dta_tidy %>% 
  filter(sex == "male") %>% 
  filter(age >= 15 & age <= 45) %>%
  mutate(lg_risk = log(risk, 10)) %>% 
  mutate(years_since_start = year - min(year)) -> dta_male_young_adult

dta_male_young_adult %>% 
    lm(lg_risk ~ as.factor(age) * (years_since_start), 
       data = .)  -> mod_lin_improvement

# predicted vs actual risk 

dta_male_young_adult %>% 
  modelr::add_predictions(mod_lin_improvement) %>% 
  dplyr::select(year, age, lg_risk, pred)  %>% 
  gather("type", "value", lg_risk:pred) %>% 
  levelplot(
    value ~ year * age | type, data = .,
    col.regions = qual_col,
    #    at = seq(from = -4.5, to = -1.5, by = 0.05),
    scales = list(alternating = 3),
    aspect= "iso"
  )

# And residuals

dta_male_young_adult %>% 
  modelr::add_predictions(mod_lin_improvement) %>% 
  dplyr::select(year, age, lg_risk, pred) %>% 
  mutate(residual = lg_risk - pred) %>% 
  levelplot(
    residual ~ year * age, data = .,
    col.regions = rdbu_col,
    at = seq(from = -0.4, to = 0.4, by = 0.05),
    scales = list(
      alternating = 3,
      x = list(at = seq(1980, 2020, by = 5), rot = 90)
    ),
    aspect= "iso"
  )



dta_tidy %>% 
  filter(sex == "female") %>% 
  filter(age >= 15 & age <= 45) %>%
  mutate(lg_risk = log(risk, 10)) %>% 
  mutate(years_since_start = year - min(year)) -> dta_male_young_adult

dta_male_young_adult %>% 
  lm(lg_risk ~ as.factor(age) * (years_since_start), 
     data = .)  -> mod_lin_improvement

# predicted vs actual risk 

dta_male_young_adult %>% 
  modelr::add_predictions(mod_lin_improvement) %>% 
  dplyr::select(year, age, lg_risk, pred)  %>% 
  gather("type", "value", lg_risk:pred) %>% 
  levelplot(
    value ~ year * age | type, data = .,
    col.regions = qual_col,
    #    at = seq(from = -4.5, to = -1.5, by = 0.05),
    scales = list(alternating = 3),
    aspect= "iso"
  )

# And residuals

dta_male_young_adult %>% 
  modelr::add_predictions(mod_lin_improvement) %>% 
  dplyr::select(year, age, lg_risk, pred) %>% 
  mutate(residual = lg_risk - pred) %>% 
  levelplot(
    residual ~ year * age, data = .,
    col.regions = rdbu_col,
    at = seq(from = -0.4, to = 0.4, by = 0.05),
    scales = list(
      alternating = 3,
      x = list(at = seq(1980, 2020, by = 5), rot = 90)
    ),
    aspect= "iso"
  )

