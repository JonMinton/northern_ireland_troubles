# Algeria data

# Jon Minton
# 27/7/2017


# To dos ------------------------------------------------------------------

# 1. Create function factory for leveplot defauls

# Load packages -----------------------------------------------------------


rm(list = ls())

# Package management ------------------------------------------------------


pacman::p_load(
  MASS,
  readxl,
  tidyverse,
  purrr, 
  rgl,  r2stl,
  RColorBrewer,
  lattice,   latticeExtra,
  spatstat
)


# Data Management ---------------------------------------------------------


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


# Levelplots of data ------------------------------------------------------


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



# 3D plots - rgl ----------------------------------------------------------


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


# stl files ---------------------------------------------------------------

convert_to_list_matrix <- function(tidy_dta){
  tidy_dta %>% 
    select(year, age, risk) %>% 
    spread(age, risk) -> tmp
  years <- tmp$year
  tmp$year <- NULL
  ages <- colnames(tmp) 
  tmp <- as.matrix(tmp)
  rownames(tmp) <- years
  colnames(tmp) <- ages
  out <- list(
    x = as.numeric(ages), 
    y = as.numeric(years), 
    z = tmp
    )
  
  out
}

dta_tidy %>% 
  filter(sex == "male") %>% 
  convert_to_list_matrix() -> out_male_all

dta_tidy %>% 
  filter(sex == "female") %>% 
  convert_to_list_matrix() -> out_female_all

dta_tidy %>% 
  filter(sex == "male") %>% 
  mutate(risk = log(risk, 10)) %>% 
  convert_to_list_matrix() -> log_male_all

dta_tidy %>% 
  filter(sex == "female") %>% 
  mutate(risk = log(risk, 10)) %>% 
  convert_to_list_matrix() -> log_female_all


r2stl(
    x = as.numeric(rownames(out_male_all$z)),
    y = as.numeric(colnames(out_male_all$z)),
    z = out_male_all$z,
    filename = "data/algeria/male.stl",
#    object.name = "tmp.stl",
    show.persp = T,
    z.expand = T
#    strict.stl = T,
#    min.height = 0.01
  )


r2stl(
  x = as.numeric(rownames(out_female_all$z)),
  y = as.numeric(colnames(out_female_all$z)),
  z = out_female_all$z,
  filename = "data/algeria/female.stl",
  #    object.name = "tmp.stl",
  show.persp = T,
  z.expand = T
  #    strict.stl = T,
  #    min.height = 0.01
)

r2stl(
  x = as.numeric(rownames(log_male_all$z)),
  y = as.numeric(colnames(log_male_all$z)),
  z = log_male_all$z,
  filename = "data/algeria/log_male.stl",
  #    object.name = "tmp.stl",
  show.persp = T,
  z.expand = T
  #    strict.stl = T,
  #    min.height = 0.01
)

r2stl(
  x = as.numeric(rownames(log_female_all$z)),
  y = as.numeric(colnames(log_female_all$z)),
  z = log_female_all$z,
  filename = "data/algeria/log_female.stl",
  #    object.name = "tmp.stl",
  show.persp = T,
  z.expand = T
  #    strict.stl = T,
  #    min.height = 0.01
)

# Modelling  --------------------------------------------------------------

k <- 0.09748423

dta_tidy %>% 
  filter(sex == "male") %>% 
  filter(age >= 15 & age <= 45) %>%
  mutate(lg_risk = log(risk, 10)) %>% 
  mutate(years_since_start = year - min(year)) %>% 
  mutate(
    conf_kpc = ifelse(year < 1994, 0, (1 - k) ^ (year - 1994))
  ) -> dta_male_young_adult

dta_male_young_adult %>% 
    lm(lg_risk ~ as.factor(age) * (years_since_start) + (year > 1982) + (year > 1997), 
       data = .)  -> mod_lin_improvement

dta_male_young_adult %>% 
  lm(lg_risk ~ as.factor(age) * years_since_start,
     data = .) -> mod_lin

dta_male_young_adult %>% 
  lm(lg_risk ~ as.factor(age) * years_since_start + (year > 1982) + (year > 1997),
     data = .) -> mod_lin_intercept_correction

dta_male_young_adult %>% 
  lm(lg_risk ~ as.factor(age) * (years_since_start + (year > 1982) + (year > 1997) ),
     data = .) -> mod_lin_age_spec_correction

dta_male_young_adult %>% 
  lm(lg_risk ~ as.factor(age) * (years_since_start + conf_kpc),
     data = .) -> mod_lin_conf

dta_male_young_adult %>% 
  lm(lg_risk ~ as.factor(age) * (years_since_start + conf_kpc) + (year > 1982) + (year > 1997),
     data = .) -> mod_lin_intercept_correction_conf

dta_male_young_adult %>% 
  lm(lg_risk ~ as.factor(age) * (years_since_start + conf_kpc + (year > 1982) + (year > 1997) ),
     data = .) -> mod_lin_age_spec_correction_conf


# F tests

anova(mod_lin, mod_lin_conf) # Conf +
anova(mod_lin_intercept_correction, mod_lin_intercept_correction_conf) # Conf +
anova(mod_lin_age_spec_correction, mod_lin_age_spec_correction_conf) # Conf +

AIC(mod_lin, mod_lin_conf,
    mod_lin_intercept_correction, mod_lin_intercept_correction_conf,
    mod_lin_age_spec_correction, mod_lin_age_spec_correction_conf)
# Best model: mod_lin_age_spec_correction_conf

BIC(mod_lin, mod_lin_conf,
    mod_lin_intercept_correction, mod_lin_intercept_correction_conf,
    mod_lin_age_spec_correction, mod_lin_age_spec_correction_conf)
# Best model: mod_lin_intercept_correction_conf



# predicted vs actual risk 

dta_male_young_adult %>% 
  modelr::add_predictions(mod_lin_intercept_correction_conf) %>% 
  mutate(int_pred = pred) %>% 
  modelr::add_predictions(mod_lin_age_spec_correction_conf) %>% 
  mutate(age_pred = pred) %>% 
  dplyr::select(year, age, lg_risk, int_pred, age_pred)  %>%
  gather("type", "value", lg_risk:age_pred) %>% 
  levelplot(
    value ~ year * age | type, data = .,
    col.regions = qual_col,
    #    at = seq(from = -4.5, to = -1.5, by = 0.05),
    scales = list(alternating = 3),
    aspect= "iso"
  )

# And residuals

dta_male_young_adult %>% 
  modelr::add_predictions(mod_lin_intercept_correction_conf) %>% 
  mutate(int_pred = pred) %>% 
  modelr::add_predictions(mod_lin_age_spec_correction_conf) %>% 
  mutate(age_pred = pred) %>% 
  dplyr::select(year, age, lg_risk, int_pred, age_pred)  %>%
  mutate(res_int = int_pred - lg_risk, res_age = age_pred - lg_risk) %>% 
  select(year, age, res_int, res_age) %>% 
  gather("type", "value", res_int:res_age) %>% 
  levelplot(
    value ~ year * age | type, data = .,
    col.regions = rdbu_col,
    at = seq(from = -0.4, to = 0.4, by = 0.05),
    scales = list(
      alternating = 3,
      x = list(at = seq(1980, 2020, by = 5), rot = 90)
    ),
    aspect= "iso"
  )



# Optimise decay rate 


k <- 0.09748423
optimise_decay_rate <- function(par){
  
  k <- par[["k"]]
  
  dta_tidy %>% 
    filter(sex == "male") %>% 
    filter(age >= 15 & age <= 45) %>%
    mutate(lg_risk = log(risk, 10)) %>% 
    mutate(years_since_start = year - min(year)) %>% 
    mutate(
      conf_kpc = ifelse(year < 1994, 0, (1 - k) ^ (year - 1994))
    ) -> dta_male_young_adult
  
  dta_male_young_adult %>% 
    lm(lg_risk ~ as.factor(age) * (years_since_start + conf_kpc + (year > 1982) + (year > 1997) ),
       data = .) -> mod_lin_age_spec_correction_conf
  
  
  AIC(mod_lin_age_spec_correction_conf)
}

opt_k <- optim(
  par = list(k = 0.0974), fn = optimise_decay_rate,
  lower = list(k = 0), upper = list(k = 2),
  method = "L-BFGS-B"
)

# Rate maximised at 0.1980003 

# Plot for many values of k

k_vals <- seq(0, 1, by = 0.01)
aics <- vector("numeric", length = length(k_vals))

for (i in seq_along(k_vals)){
  this_k <- k_vals[i]
  this_par <- list(k = this_k)
  aics[i] <- optimise_decay_rate(par = this_par)  
}

tmp_df <- 
  data_frame(k = k_vals, aic = aics)

tmp_df %>% ggplot(., aes(x = k, y = aic)) + 
  geom_line()

