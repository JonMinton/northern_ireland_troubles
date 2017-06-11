
# Northern Ireland Mortality exploration 

# Jon Minton
# 7/6/2017


# To dos ------------------------------------------------------------------

# 1. Create function factory for leveplot defauls

# Load packages -----------------------------------------------------------


rm(list = ls())

pacman::p_load(
  MASS,
  tidyverse,
  purrr, rgl,
  RColorBrewer,
  lattice,   latticeExtra,
  spatstat
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
rdbu_col <- colorRampPalette(rev(brewer.pal(5, "RdBu")))(200)

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

# Ages up to 90, log_10 scale
dta_nir %>% 
  filter(age <= 90) %>% 
  levelplot(
    lmr ~ year * age | sex, data = .,
    col.regions = qual_col,
    at = seq(from = -5.0, to = 0.0, by = 0.1),
    scales=list(alternating=3),
    aspect= "iso"
  ) -> fig_01

tiff("figures/fig01_log_both_genders_all_ages.tiff", width = 15, height = 10, units = "cm", res = 300)
print(fig_01)
dev.off()

# Blank figure

dta_nir %>% 
  filter(age <= 90) %>% 
  filter(sex == "male") %>%
  mutate(lmr = NA) %>% 
  levelplot(
    lmr ~ year * age, data = .,
    col.regions = qual_col,
    at = seq(from = -5.0, to = 0.0, by = 0.1),
    scales=list(alternating=3),
    colorkey = F,
    aspect= "iso"
  ) -> fig_02
print(fig_02)


svg("figures/fig02_log_both_genders_all_ages.svg")
print(fig_02)
dev.off()

# Figure 3: deaths between ages 18 and 40 

dta_nir %>% 
  filter(age >= 18 & age <= 40) %>% 
  group_by(year, sex) %>% 
  summarise(deaths = sum(deaths)) %>% 
  ggplot(aes(x = year, y = deaths, group = sex, colour = sex, shape = sex)) + 
  geom_point() + 
  theme_minimal() + 
  labs(x = "Year", y = "Deaths") + 
  geom_rect(
    data = data.frame(xmin = 1971, xmax = 1973, ymin = -Inf, ymax = Inf),
    aes(xmin = xmin, xmax = xmax, ymin=ymin, ymax = ymax), fill = "lightgrey", alpha = 0.5,
    inherit.aes = F
    ) + 
  theme(legend.position = c(0.8, 0.8))

ggsave("figures/fig03_deaths_trend.png", dpi = 300, width = 10, height = 10, units = "cm")

# RGL figure

dta_nir %>% 
  filter(age <= 90) %>% 
  filter(sex == "male") %>% 
  select(year, age, mr) %>% 
  spread(age, mr) -> tmp

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

## FEMALES 
dta_nir %>% 
  filter(age <= 90) %>% 
  filter(sex == "female") %>% 
  select(year, age, mr) %>% 
  spread(age, mr) -> tmp

yrs <- tmp$year
tmp$year <- NULL
ages <- names(tmp)
tmp <- as.matrix(tmp)

# Persp plot, females, all ages 
persp3d(z = tmp, col = "lightgray", 
        box = F, axes = F,
        aspect = c(1, 1, 91/93),
        xlab = "", ylab = "", zlab = ""
)


## Males - log scale, all ages
dta_nir %>% 
  filter(age <= 90) %>% 
  filter(sex == "male") %>% 
  mutate(mr = (deaths + 0.5) / (exposure + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  select(year, age, lmr) %>%
  mutate(lmr = ifelse(is.na(lmr), mean(lmr, na.rm = T), lmr)) %>% 
#  mutate(lmr = scales::rescale(lmr, to = c(0, 1))) %>% 
  spread(age, lmr) -> tmp

yrs <- tmp$year
tmp$year <- NULL
ages <- names(tmp)
tmp <- as.matrix(tmp)
tmp <- tmp[1:92, 1:91]

# Persp plot, males, all ages 
persp3d(z = tmp, col = "lightgray", 
        box = F, axes = F,
#        aspect = c(1, 1, 91/93),
        xlab = "", ylab = "", zlab = ""
)

dta_nir %>% 
  filter(age <= 90) %>% 
  filter(sex == "female") %>% 
  mutate(mr = (deaths + 0.5) / (exposure + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  select(year, age, lmr) %>%
  mutate(lmr = ifelse(is.na(lmr), mean(lmr, na.rm = T), lmr)) %>% 
  #  mutate(lmr = scales::rescale(lmr, to = c(0, 1))) %>% 
  spread(age, lmr) -> tmp

yrs <- tmp$year
tmp$year <- NULL
ages <- names(tmp)
tmp <- as.matrix(tmp)
tmp <- tmp[1:92, 1:91]

# Persp plot, males, all ages 
persp3d(z = tmp, col = "lightgray", 
        box = F, axes = F,
        #        aspect = c(1, 1, 91/93),
        xlab = "", ylab = "", zlab = ""
)




# Ages up to 90, log scale
dta_nir %>% 
  filter(age <= 90) %>% 
  levelplot(
    lmr ~ P * age | sex, data = .,
    col.regions = qual_col,
    at = seq(from = -5.0, to = 0.1, by = 0.1),
    scales=list(alternating=3),
    aspect= "iso"
  )


# Ages between 15 and 45, linear scale

dta_nir %>% 
  filter(age >= 15 & age <= 45) %>% 
  levelplot(
    mr ~ P * age | sex, data = .,
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
  ) -> fig04 

tiff("figures/fig04_lmr_zoomed.tiff", height = 15, width = 15, units = "cm", res = 300)
print(fig04)
dev.off()



dta_nir %>% 
  filter(age >= 15 & age <= 45) %>% 
  filter(sex == "male") %>% 
  select(year, age, mr) %>% 
  spread(age, mr) -> tmp

yrs <- tmp$year
tmp$year <- NULL
ages <- names(tmp)
tmp <- as.matrix(tmp)

# Persp plot, males, all ages 
persp3d(z = tmp, col = "lightgray", 
        box = F, axes = F,
        aspect = c(1, 1, 31/93),
        xlab = "", ylab = "", zlab = ""
)


# Linear Model Specs ------------------------------------------------------


# Log mort rate, based on age
# year pre 1938 flag 
# ww2 flag
# 45-55 flag
# after 55 flag

k <- 0.09748423

dta_nir %>% 
  filter(sex == "male") %>% 
  filter(age >= 15 & age <= 45) %>% 
  dplyr::select(year, age, lmr) %>% 
  mutate(yr_pre_38 = year <= 1938) %>%
  mutate(
    early_phase = ifelse(
      yr_pre_38,
      year - 1922,
      0
    )
  ) %>% 
  mutate(impr_phase = year >= 1939 & year <= 1955) %>% 
  mutate(
    yr_since_imp = ifelse(
      impr_phase, 
      year - 1939,
      0
    )
  ) %>% 
  mutate(
    post_55 = year > 1955,
    yr_since_imp2 = ifelse(
      post_55, 
      year - 1955,
      0
    )
  ) %>% 
  mutate(
    trbls_kpc = ifelse(year < 1972, 0, (1 - k) ^ (year - 1972))
  ) -> dta_nir_sub_lmr



dta_nir_sub_lmr %>% 
  lm(lmr ~ as.factor(age) * (yr_pre_38 + impr_phase + yr_since_imp), 
     data = .) -> mod_lmr_01


dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_01) %>% 
  dplyr::select(year, age, lmr, pred )  %>% 
  gather("type", "value", lmr:pred) %>% 
  levelplot(
    value ~ year * age | type, data = .,
    col.regions = qual_col,
    #    at = seq(from = -4.5, to = -1.5, by = 0.05),
    scales = list(alternating = 3),
    aspect= "iso"
  )

# Now to add Troubles w/ 9.76% decay 
dta_nir_sub_lmr %>% 
  lm(lmr ~ as.factor(age) * (yr_pre_38 + impr_phase + yr_since_imp + trbls_kpc), 
     data = .) -> mod_lmr_02

dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_02) %>% 
  dplyr::select(year, age, lmr, pred )  %>% 
  gather("type", "value", lmr:pred) %>% 
  levelplot(
    value ~ year * age | type, data = .,
    col.regions = qual_col,
    #    at = seq(from = -4.5, to = -1.5, by = 0.05),
    scales = list(alternating = 3),
    aspect= "iso"
  )

# Now to look at the residuals 

dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_02) %>% 
  dplyr::select(year, age, lmr, pred ) %>% 
  mutate(residuals = pred - lmr) %>% 
  levelplot(
    residuals ~ year * age, data = .,
    col.regions = rdbu_col,
    at = seq(from = -0.6, to = 0.6, by = 0.05),
    scales = list(
      alternating = 3,
      x = list(at = seq(1925, 2010, by = 5), rot = 90)
    ),
    aspect= "iso"
  )

# Looks like there's another smaller trend after 1955

dta_nir_sub_lmr %>% 
  lm(lmr ~ as.factor(age) * (yr_pre_38 + impr_phase + yr_since_imp + yr_since_imp2), 
     data = .) -> mod_lmr_03

dta_nir_sub_lmr %>% 
  lm(lmr ~ as.factor(age) * (yr_pre_38 + impr_phase + yr_since_imp + yr_since_imp2 + trbls_kpc), 
     data = .) -> mod_lmr_04

dta_nir_sub_lmr %>% 
  lm(lmr ~ as.factor(age) * (yr_pre_38 + early_phase + impr_phase + yr_since_imp + yr_since_imp2 + trbls_kpc), 
     data = .) -> mod_lmr_05


dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_05) %>% 
  dplyr::select(year, age, actual = lmr , predicted = pred)  %>% 
  gather("type", "value", actual:predicted) %>% 
  levelplot(
    value ~ year * age | type, data = .,
    col.regions = qual_col,
    at = seq(from = -4.5, to = -1.5, by = 0.05),
    scales = list(alternating = 3),
    aspect= "iso"
  ) -> fig05

tiff("figures/fig05_prediction_surface.tiff", height = 10, width = 20, units = "cm", res = 300)
print(fig05)
dev.off()


dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_05) %>% 
  dplyr::select(year, age, lmr, pred ) %>% 
  mutate(residuals = pred - lmr) %>% 
  levelplot(
    residuals ~ year * age, data = .,
    col.regions = rdbu_col,
    at = seq(from = -0.6, to = 0.6, by = 0.05),
    scales = list(
      alternating = 3,
      x = list(at = seq(1925, 2010, by = 5), rot = 90)
    ),
    aspect= "iso"
  ) -> fig06

tiff("figures/fig06_residuals_surface.tiff", height = 10, width = 20, units = "cm", res = 300)
print(fig06)
dev.off()




dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_04) %>% 
  dplyr::select(year, age, lmr, pred ) %>% 
  mutate(residuals = pred - lmr) %>% 
  levelplot(
    residuals ~ year * age, data = .,
    col.regions = rdbu_col,
    at = seq(from = -0.6, to = 0.6, by = 0.05),
    scales = list(
      alternating = 3,
      x = list(at = seq(1925, 2010, by = 5), rot = 90)
    ),
    aspect= "iso"
  )


dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_05) %>% 
  dplyr::select(year, age, lmr, pred ) %>% 
  mutate(residuals = pred - lmr) %>% 
  levelplot(
    residuals ~ year * age, data = .,
    col.regions = rdbu_col,
    at = seq(from = -0.6, to = 0.6, by = 0.05),
    scales = list(
      alternating = 3,
      x = list(at = seq(1925, 2010, by = 5), rot = 90)
    ),
    aspect= "iso"
  )

# So, model 4 has the BEST AIC, and also has the most random looking pattern of residual scatter. 

# Let's look at RMS for each 

get_rms <- function(res){
  out <- res^2
  out <- mean(out, na.rm = T)
  out <- out^0.5
  out
}

dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_01) %>% 
  mutate(residuals = pred - lmr) %>% .$residuals -> tmp1
rms_01 <- get_rms(tmp1)

dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_02) %>% 
  mutate(residuals = pred - lmr) %>% .$residuals -> tmp2
rms_02 <- get_rms(tmp2)

dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_03) %>% 
  mutate(residuals = pred - lmr) %>% .$residuals -> tmp3
rms_03 <- get_rms(tmp3)

dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_04) %>% 
  mutate(residuals = pred - lmr) %>% .$residuals -> tmp4
rms_04 <- get_rms(tmp4)

rms_04 / rms_01
# About a 17% fall in RMS error

dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_05) %>% 
  mutate(residuals = pred - lmr) %>% .$residuals -> tmp5
rms_05 <- get_rms(tmp5)

rms_05 / rms_01

# About a 19% fall in RMS error




dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_05, var = "pred_active") %>%
  dplyr::select(year, age, lmr, pred_active) -> preds_trbls

dta_nir_sub_lmr %>% 
  mutate(trbls_kpc = 0) -> newdat 

broom::augment(mod_lmr_05, newdata = newdat) %>% 
  dplyr::select(year, age, pred_notrbls = .fitted) %>% 
  tbl_df -> preds_notrbls

preds_trbls %>% 
  inner_join(preds_notrbls) %>% 
  gather("type", "value", lmr:pred_notrbls) -> preds_to_compare

preds_to_compare %>% 
  levelplot(
    value ~ year * age | type, data = .,
    col.regions = qual_col,
    #    at = seq(from = -4.5, to = -1.5, by = 0.05),
    scales = list(alternating = 3),
    aspect= "iso"
  )
  


# optimise decay rate -----------------------------------------------------




find_decay_rate <- function(k = 0.04, what = "AIC"){
  dta_nir %>% 
    filter(sex == "male") %>% 
    filter(age >= 15 & age <= 45) %>% 
    dplyr::select(year, age, lmr) %>% 
    mutate(yr_pre_38 = year <= 1938) %>%
    mutate(
      early_phase = ifelse(
        yr_pre_38,
        year - 1922,
        0
      )
    ) %>% 
    mutate(impr_phase = year >= 1939 & year <= 1955) %>% 
    mutate(
      yr_since_imp = ifelse(
        impr_phase, 
        year - 1939,
        0
      )
    ) %>% 
    mutate(
      post_55 = year > 1955,
      yr_since_imp2 = ifelse(
        post_55, 
        year - 1955,
        0
      )
    ) %>% 
    mutate(
      trbls_kpc = ifelse(year < 1972, 0, (1 - k) ^ (year - 1972))
    ) -> this_dta
  
  this_dta %>% 
    lm(lmr ~ as.factor(age) * (yr_pre_38 + early_phase + impr_phase + yr_since_imp + yr_since_imp2 + trbls_kpc), 
       data = .) -> this_model
  if (what == "AIC"){
    out <- AIC(this_model)
  } else if (what == "BIC"){
    out <- BIC(this_model)
  }
  
  return(out) 
}

dec_find <- optimize(find_decay_rate, c(0, 1))

# Best AIC at 0.09748423

dec_find <- optimize(find_decay_rate, c(0, 1), what = "BIC")

# Best BIC at 0.09748423


decay_aic <- data_frame(
  decay_rates = seq(0, 1, by = 0.005)
  ) %>% 
  mutate(
    aic = map_dbl(decay_rates, find_decay_rate)
  )

ggplot(decay_aic) +
  geom_line(aes(x = decay_rates, y = aic)) + 
  theme_minimal() + 
  labs(x = "Decay Rate", y = "Penalised model fit (AIC)")

ggsave("figures/fig07_fit_decay.tiff", height = 8, width = 8, units = "cm", dpi = 300)


# Half life 

log(0.5) / log (1 - 0.09748423)
# 6.76 years



# Estimate impact of counterfactual ---------------------------------------


dta_nir_sub_lmr %>% 
  modelr::add_predictions(mod_lmr_05, var = "pred_active") %>%
  dplyr::select(year, age, lmr, pred_active) -> preds_trbls

dta_nir_sub_lmr %>% 
  mutate(trbls_kpc = 0) -> newdat 

broom::augment(mod_lmr_05, newdata = newdat) %>% 
  dplyr::select(year, age, pred_notrbls = .fitted) %>% 
  tbl_df -> preds_notrbls

preds_trbls %>% 
  inner_join(preds_notrbls) %>% 
  gather("type", "value", lmr:pred_notrbls) -> preds_to_compare

dta_nir %>% 
  filter(sex == "male") %>% 
  select(year, age, deaths, exposure) %>% 
  right_join(preds_to_compare) %>% 
  mutate(pred_deaths = exposure * 10 ^ value) %>% 
  select(-value) %>% 
  spread(type, pred_deaths) %>% 
  select(-lmr) %>% 
  mutate(difference = pred_active - pred_notrbls) -> dif_estimates


dif_estimates %>% 
  filter(year >= 1970) %>% 
  levelplot(
    difference ~ year * age , data = .,
    col.regions = qual_col,
    #    at = seq(from = -4.5, to = -1.5, by = 0.05),
    scales = list(alternating = 3),
    aspect= "iso"
  ) -> fig_08

tiff("figures/fig_08_trbls_effect.tiff", height = 8, width = 10, units = "cm", res = 300)
print(fig_08)
dev.off()


# Table 1 about here 

dif_estimates %>% 
  mutate(age_grp = cut(age, breaks = seq(15, 45, by = 5), include.lowest = T)) %>% 
  filter(year >= 1970) %>% 
  group_by(year, age_grp) %>%
  summarise(excess_deaths = sum(difference)) %>%
  spread(age_grp, excess_deaths) %>% 
  mutate(total = `[15,20]` +`(20,25]` +`(25,30]`+ `(30,35]`+ `(35,40]`+ `(40,45]`) %>% 
  mutate_at(vars(`[15,20]`:total), funs(round))







# 3d plot 
dif_estimates %>% 
  select(year, age, difference) %>% 
  spread(age, difference) -> tmp

yrs <- tmp$year
tmp$year <- NULL
ages <- names(tmp)
tmp <- as.matrix(tmp)
dimnames(tmp) <- list(yrs, ages)

persp3d(z = tmp, col = "lightgray", 
        box = F, axes = F,
        #        aspect = c(1, 1, 91/93),
        xlab = "", ylab = "", zlab = ""
)


# So, how many extra deaths in total?

dif_estimates %>% 
  arrange(year) %>% 
  group_by(year) %>% 
  summarise(difference = sum(difference)) %>% 
  ggplot(.) + 
  geom_line(aes(x = year, y = difference))


dif_estimates %>% 
  arrange(year) %>% 
  group_by(year) %>% 
  summarise(difference = sum(difference)) %>% 
  ungroup() %>% 
  mutate(cum_difference = cumsum(difference)) %>% 
  ggplot(.) + 
  geom_line(aes(x = year, y = cum_difference))

dif_estimates %>% 
  arrange(year) %>% 
  group_by(year) %>% 
  summarise(difference = sum(difference)) %>% 
  ungroup() %>% 
  mutate(cum_difference = cumsum(difference)) %>% 
  arrange(rev(year))

# A total of 2776 additional deaths predicted using this model alone
# So around 50-60% of the attributed deaths appear explained by this 
# stylised model 



# Model 5 but for females

k <- 0.09748423

dta_nir %>% 
  filter(sex == "female") %>% 
  filter(age >= 15 & age <= 45) %>% 
  mutate(mr = (deaths + 0.5) / (exposure + 0.5)) %>% 
  mutate(lmr = log(mr, 10)) %>% 
  dplyr::select(year, age, lmr) %>% 
  mutate(yr_pre_38 = year <= 1938) %>%
  mutate(
    early_phase = ifelse(
      yr_pre_38,
      year - 1922,
      0
    )
  ) %>% 
  mutate(impr_phase = year >= 1939 & year <= 1955) %>% 
  mutate(
    yr_since_imp = ifelse(
      impr_phase, 
      year - 1939,
      0
    )
  ) %>% 
  mutate(
    post_55 = year > 1955,
    yr_since_imp2 = ifelse(
      post_55, 
      year - 1955,
      0
    )
  ) %>% 
  mutate(
    trbls_kpc = ifelse(year < 1972, 0, (1 - k) ^ (year - 1972))
  ) -> dta_nir_sub_lmr_female


dta_nir_sub_lmr_female %>% 
  lm(lmr ~ as.factor(age) * (yr_pre_38 + early_phase + impr_phase + yr_since_imp + yr_since_imp2 + trbls_kpc), 
     data = .) -> mod_lmr_05_female

dta_nir_sub_lmr_female %>% 
  modelr::add_predictions(mod_lmr_05_female, var = "pred_active") %>%
  dplyr::select(year, age, lmr, pred_active) -> preds_trbls_female

preds_trbls_female %>% 
  gather(key = "type", value = "value", lmr:pred_active) %>% 
  levelplot(
    value ~ year * age | type, data = .,
    col.regions = qual_col,
  #    at = seq(from = -4.5, to = -1.5, by = 0.05),
    scales = list(alternating = 3),
    aspect= "iso"
  )

dta_nir_sub_lmr_female %>% 
  lm(lmr ~ as.factor(age) * (yr_pre_38 + early_phase + impr_phase + yr_since_imp + yr_since_imp2), 
     data = .) -> mod_lmr_05_female_notrbls

AIC(mod_lmr_05_female, mod_lmr_05_female_notrbls)


coef(mod_lmr_05_female) -> trbls_effects_females
trbls_effects_females <- trbls_effects_females[188:217]
tmp <- data_frame(age = 16:45, female = trbls_effects_females)

# Troubles effects by age 
coef(mod_lmr_05) -> trbls_effects
trbls_effects <- trbls_effects[188:217]
tmp2 <- data_frame(age = 16:45, male = trbls_effects)

tmp <- full_join(tmp, tmp2) %>% 
  gather(sex, effect, - age)

tmp %>% ggplot() + 
  geom_point(aes(x = age, y = effect, group = sex, colour = sex, shape = sex)) +
  theme_minimal() + 
  labs(x = "Age", y = "Troubles Effect Coefficient") + 
  geom_vline(xintercept = 18, linetype= "dashed") +
  geom_hline(yintercept = 0)

ggsave("figures/fig09_troubles_coeff.tiff", height = 12, width = 12, units = "cm", dpi = 300)



