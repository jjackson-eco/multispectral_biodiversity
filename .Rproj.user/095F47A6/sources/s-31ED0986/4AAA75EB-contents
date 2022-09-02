##############################################################
##                                                          ##
##    RainDrop: Multispectral indicators of Biodiversity    ##
##                                                          ##
##  Biodiversity vs spectral diversity - Regression models  ## 
##                                                          ##
##              All Uncalibrated data July 2021             ##
##                                                          ##
##                        8th Feb 2022                      ##
##                                                          ##
##############################################################
# updated 2022-09-02
rm(list = ls())

# lib paths for uni computer
source("../lib_path_uni_computer.R")

library(tidyverse)
library(patchwork)
library(flextable)
library(ggridges)
library(ggdist)
library(cowplot)
library(brms)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Loading data ####

load("data/biodiversity_data_early.RData", verbose = TRUE)
load("data/multispec_data_uncalibrated.RData", verbose = TRUE)
load("../../RainDropRobotics/Data/raindrop_biodiversity_2016_2020.RData",
     verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Merging data ####

# Biomass summary
biomass_sum <- biomass_july2021 %>% 
  group_by(block, treatment, experiment) %>% 
  mutate(tot_biomass = sum(biomass)) %>% 
  ungroup() %>% 
  filter(is.na(functional_group) == FALSE) %>% 
  pivot_wider(names_from = functional_group, values_from = biomass) 

# converting 0's to NA's so they are ommited
biomass_sum[biomass_sum == 0] <- NA

# Merge for overall plot/height summary
multispectral_species_diversity <- multispec_data_uncalibrated %>% 
  left_join(x = ., y = percent_cover_july2021, by = c("block", "treatment", "experiment")) %>% 
  left_join(x = ., y = biomass_sum, by = c("block", "treatment", "experiment")) %>% 
  group_by(block, treatment, experiment, height, replicate) %>% 
  summarise(cv = mean(cv), sd = mean(sd), 
            skew = mean(skew), kurtosis = mean(kurtosis),
            richness = richness[1], simpsons = simpsons[1], 
            shannon = shannon[1], tot_biomass = tot_biomass[1],
            graminoids = Graminoids[1], legumes = Legumes[1], 
            forbs = Forbs[1]) %>% 
  ungroup() %>% 
  filter(height < 10) %>%  # remove 10m as it is only one image
  mutate(height = 
           case_when(
             height <= 3 ~ 2,
             height <= 5 & height > 3 ~ 4,
             height <= 7 & height > 5 ~ 6,
             height <= 9 & height > 7 ~ 8
           ), #adjust to more approximate heights to account for odd sample sizes
         tot_biomass = if_else(experiment == "DROUGHTNET", tot_biomass*4, tot_biomass*5)) #convert biomass to 1 m2 depending on experiment 

# Merge for full data across bands 
multispectral_species_diversity_allbands <- multispec_data_uncalibrated %>% 
  left_join(x = ., y = percent_cover_july2021, by = c("block", "treatment", "experiment")) %>% 
  left_join(x = ., y = biomass_sum, by = c("block", "treatment", "experiment"))  %>% 
  group_by(block, treatment, experiment, height, replicate, band) %>% 
  summarise(band_nm = band_nm[1], cv = mean(cv), sd = mean(sd), 
            skew = mean(skew), kurtosis = mean(kurtosis),
            richness = richness[1], simpsons = simpsons[1], 
            shannon = shannon[1], tot_biomass = tot_biomass[1],
            graminoids = Graminoids[1], legumes = Legumes[1], 
            forbs = Forbs[1]) %>% 
  ungroup() %>% 
  filter(height < 10) %>%  # remove 10m as it is only one image
  mutate(height = 
           case_when(
             height <= 3 ~ 2,
             height <= 5 & height > 3 ~ 4,
             height <= 7 & height > 5 ~ 6,
             height <= 9 & height > 7 ~ 8
           )) #adjust to more approximate heights to account for odd sample sizes

# CSV data, for future work
multispectral_species_diversity_full <- multispec_data_uncalibrated %>% 
  left_join(x = ., y = percent_cover_july2021, by = c("block", "treatment", "experiment")) %>% 
  left_join(x = ., y = biomass_sum, by = c("block", "treatment", "experiment")) %>% 
  dplyr::select(file, id, image_date = rec_date, experiment, block, 
                treatment, height, band, band_nm, replicate, cv, sd, skew, 
                kurtosis, date_cover = date.x, richness, simpsons, shannon,
                date_biomass = date.y, tot_biomass, graminoids = Graminoids,
                forbs = Forbs, legumes = Legumes, bryophytes = Bryophytes, 
                woody = Woody) %>% 
  mutate(height = 
           case_when(
             height <= 3 ~ 2,
             height <= 5 & height > 3 ~ 4,
             height <= 7 & height > 5 ~ 6,
             height <= 9 & height > 7 ~ 8
           )) #adjust to more approximate heights to account for odd sample sizes

write_csv(multispectral_species_diversity_full, 
          file = "data/multispec_biodiversity_uncalibrated.csv")

## Diversity across plots summary
summary(multispectral_species_diversity$richness)
summary(multispectral_species_diversity$tot_biomass)

## Summary statistics
multispectral_species_diversity %>% 
  mutate(plot = paste0(block,treatment,experiment)) %>% 
  group_by(plot) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = n))+ geom_histogram()


multispectral_species_diversity_full %>% 
  mutate(plot = paste0(block,treatment,experiment)) %>% 
  group_by(plot) %>% 
  summarise(experiment = experiment[1],n = n()) %>% 
  group_by(experiment) %>% 
  mutate(n_exp = n())

multispectral_species_diversity %>% 
  mutate(plot = paste0(block,treatment,experiment)) %>% 
  group_by(plot) %>% 
  summarise(experiment = experiment[1], treatment = treatment[1], n = n()) %>% 
  group_by(treatment) %>% 
  summarise(n = n())

## Summary of graminoids percent cover
percent_cover %>% 
  filter(species_level == "Yes" & month == "June") %>% 
  group_by(year, block, treatment) %>% 
  mutate(cum_perc = (percent_cover/sum(percent_cover))*100) %>% 
  ungroup() %>% 
  group_by(year, block, treatment, plant_type) %>% 
  summarise(tot_cum_perc = sum(cum_perc)) %>% 
  ungroup() %>% 
  group_by(plant_type) %>% 
  summarise(mn = mean(tot_cum_perc, na.rm = TRUE))

## Summary of graminoids biomass
biomass %>% 
  filter(harvest == "Mid") %>% 
  group_by(year, block, treatment) %>% 
  mutate(biomass_perc = (biomass_g/sum(biomass_g))*100) %>% 
  ungroup() %>% 
  group_by(group) %>% 
  summarise(mn = mean(biomass_perc))

## Temperature and precipitation
ECN <- read_csv("../RainDropRobotics/Data/ECN_MA_Automated_Weather_Station_data_2016_2020_compatible_with_ECN_Database_headers.csv") %>% 
  dplyr::select(year = 1, day = 2, hour = 3, temp = 7, rainfall = 10) %>% 
  filter(temp >= -274, rainfall >= 0) %>% 
  group_by(year, day) %>% 
  summarise(temp = mean(temp), rainfall = sum(rainfall))

summary(ECN$temp)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Model selection ####

## Data convert for bayesian models
mssp_data <- multispectral_species_diversity %>% 
  dplyr::select(experiment, block, treatment, height, replicate, 
                shannon, simpsons, tot_biomass, 
                skew, coefficient_variation = cv) %>% 
  mutate(plot = paste0(experiment, "_", block, "_", treatment),
         height_s = as.numeric(scale(height)),
         shannon = as.numeric(scale(shannon)),
         simpsons = as.numeric(scale(simpsons)),
         tot_biomass = as.numeric(scale(tot_biomass)),
         skew = as.numeric(scale(skew)),
         coefficient_variation = as.numeric(scale(coefficient_variation)))

## Priors
ms_priors_simple <- c(prior(normal(0, 1), class = Intercept),
                      prior(exponential(2), class = sd, group = "plot"))
ms_priors <- c(prior(normal(0, 1), class = Intercept),
               prior(normal(0,1), class = b),
               prior(exponential(2), class = sd, group = "plot"))

##______________________________________________________________________________
# 3a. CV and Shannon- core question

## Base model
set.seed(666)
cv_base <- brm(
  coefficient_variation ~ 1 + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors_simple,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV height
set.seed(666)
cv_height <- brm(
  coefficient_variation ~ 1 + height + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV shannon simple
set.seed(666)
cv_shannon_simple <- brm(
  coefficient_variation ~ 1 + shannon + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV shannon
set.seed(666)
cv_shannon <- brm(
  coefficient_variation ~ 1 + shannon + height + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV shannon full
set.seed(666)
cv_shannon_full <- brm(
  coefficient_variation ~ 1 + shannon*height + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Model comparisons
cv_base <- add_criterion(cv_base, criterion = c("loo","waic"))
cv_height <- add_criterion(cv_height, criterion = c("loo","waic"))
cv_shannon_simple <- add_criterion(cv_shannon_simple, criterion = c("loo","waic"))
cv_shannon <- add_criterion(cv_shannon, criterion = c("loo","waic"))
cv_shannon_full <- add_criterion(cv_shannon_full, criterion = c("loo","waic"))

as.data.frame(loo_compare(cv_base, cv_height, cv_shannon_simple, 
                          cv_shannon, cv_shannon_full, criterion = "loo"))

plot(cv_shannon)

##______________________________________________________________________________
# 3b. CV and Simpsons

## CV simpsons simple
set.seed(666)
cv_simpsons_simple <- brm(
  coefficient_variation ~ 1 + simpsons + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV simpsons
set.seed(666)
cv_simpsons <- brm(
  coefficient_variation ~ 1 + simpsons + height + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV simpsons full
set.seed(666)
cv_simpsons_full <- brm(
  coefficient_variation ~ 1 + simpsons*height + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Model comparisons
cv_simpsons_simple <- add_criterion(cv_simpsons_simple, criterion = c("loo","waic"))
cv_simpsons <- add_criterion(cv_simpsons, criterion = c("loo","waic"))
cv_simpsons_full <- add_criterion(cv_simpsons_full, criterion = c("loo","waic"))

as.data.frame(loo_compare(cv_base, cv_height, cv_simpsons_simple, 
                          cv_simpsons, cv_simpsons_full, criterion = "loo"))

plot(cv_simpsons)

##______________________________________________________________________________
# 3c. CV and Simpsons by band

## Data convert for bayesian models
mssp_bands_data <- multispectral_species_diversity_allbands %>% 
  dplyr::select(experiment,block, treatment, band, band_nm, height, 
                replicate, shannon, simpsons, tot_biomass, 
                skew, coefficient_variation = cv) %>% 
  mutate(plot = paste0(experiment, "_", block, "_", treatment),
         height_s = as.numeric(scale(height)),
         shannon = as.numeric(scale(shannon)),
         simpsons = as.numeric(scale(simpsons)),
         tot_biomass = as.numeric(scale(tot_biomass)),
         skew = as.numeric(scale(skew)),
         coefficient_variation = as.numeric(scale(coefficient_variation)))

# priors
ms_priors_simple_band <- c(prior(normal(0, 0.7), class = Intercept),
                           prior(exponential(4), class = sd, group = "plot"))

ms_priors_band <- c(prior(normal(0, 0.7), class = Intercept),
                    prior(normal(0,0.7), class = b),
                    prior(exponential(4), class = sd, group = "plot"))


## CV with all data base
set.seed(666)
cvb_base <- brm(
  coefficient_variation ~ 1 + (1|plot),
  data = mssp_bands_data, family = gaussian(), prior = ms_priors_simple_band,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV simpsons with all data
set.seed(666)
cvb_simpsons <- brm(
  coefficient_variation ~ 1 + simpsons + height + (1|plot),
  data = mssp_bands_data, family = gaussian(), prior = ms_priors_band,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV band
set.seed(666)
cvb_band <- brm(
  coefficient_variation ~ 1 + simpsons + height + band + (1|plot),
  data = mssp_bands_data, family = gaussian(), prior = ms_priors_band,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV band full
set.seed(666)
cvb_band_full <- brm(
  coefficient_variation ~ 1 + simpsons*band + height + (1|plot),
  data = mssp_bands_data, family = gaussian(), prior = ms_priors_band,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

cvb_base <- add_criterion(cvb_base, criterion = c("loo","waic"))
cvb_simpsons <- add_criterion(cvb_simpsons, criterion = c("loo","waic"))
cvb_band <- add_criterion(cvb_band, criterion = c("loo","waic"))
cvb_band_full<- add_criterion(cvb_band_full, criterion = c("loo","waic"))

as.data.frame(loo_compare(cvb_base, cvb_simpsons, cvb_band, cvb_band_full, criterion = "loo"))

plot(cvb_band_full)

cvb_band_full

##______________________________________________________________________________
# 3d. CV and Shannon by band

## CV shannon with all data
set.seed(666)
cvb_shannon <- brm(
  coefficient_variation ~ 1 + shannon + height + (1|plot),
  data = mssp_bands_data, family = gaussian(), prior = ms_priors_band,
  chains = 4, cores = 4, iter = 3000, warmup = 1500
)

## CV band
set.seed(666)
cvb_band_shannon <- brm(
  coefficient_variation ~ 1 + shannon + height + band + (1|plot),
  data = mssp_bands_data, family = gaussian(), prior = ms_priors_band,
  chains = 4, cores = 4, iter = 3000, warmup = 1500
)

## CV band full
set.seed(666)
cvb_band_full_shannon <- brm(
  coefficient_variation ~ 1 + shannon*band + height + (1|plot),
  data = mssp_bands_data, family = gaussian(), prior = ms_priors_band,
  chains = 4, cores = 4, iter = 3000, warmup = 1500
)

cvb_shannon <- add_criterion(cvb_shannon, criterion = c("loo","waic"))
cvb_band_shannon <- add_criterion(cvb_band_shannon, criterion = c("loo","waic"))
cvb_band_full_shannon <- add_criterion(cvb_band_full_shannon, criterion = c("loo","waic"))

as.data.frame(loo_compare(cvb_base, cvb_shannon, cvb_band_shannon, cvb_band_full_shannon, criterion = "loo"))

plot(cvb_band_full_shannon)

cvb_band_full_shannon

##______________________________________________________________________________
# 3e. Skewness and Biomass

## Base model
set.seed(666)
sk_base <- brm(
  skew ~ 1 + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors_simple,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Height only
set.seed(666)
sk_height <- brm(
  skew ~ 1 + height + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Biomass simple
set.seed(666)
sk_biomass_simple <- brm(
  skew ~ 1 + tot_biomass + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Biomass 
set.seed(666)
sk_biomass <- brm(
  skew ~ 1 + tot_biomass + height + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Biomass full
set.seed(666)
sk_biomass_full <- brm(
  skew ~ 1 + tot_biomass*height + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Model comparisons
sk_base <- add_criterion(sk_base, criterion = c("loo","waic"))
sk_height <- add_criterion(sk_height, criterion = c("loo","waic"))
sk_biomass_simple <- add_criterion(sk_biomass_simple, criterion = c("loo","waic"))
sk_biomass <- add_criterion(sk_biomass, criterion = c("loo","waic"))
sk_biomass_full <- add_criterion(sk_biomass_full, criterion = c("loo","waic"))

as.data.frame(loo_compare(sk_base, sk_height, sk_biomass_simple, 
                          sk_biomass, sk_biomass_full, criterion = "loo"))

plot(sk_biomass)

##______________________________________________________________________________
# 3e. Skewness and Shannon

## shannon simple
set.seed(666)
sk_shannon_simple <- brm(
  skew ~ 1 + shannon + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## shannon 
set.seed(666)
sk_shannon <- brm(
  skew ~ 1 + shannon + height + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## shannon full
set.seed(666)
sk_shannon_full <- brm(
  skew ~ 1 + shannon*height + (1|plot),
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Model comparisons
sk_shannon_simple <- add_criterion(sk_shannon_simple, criterion = c("loo","waic"))
sk_shannon <- add_criterion(sk_shannon, criterion = c("loo","waic"))
sk_shannon_full <- add_criterion(sk_shannon_full, criterion = c("loo","waic"))

as.data.frame(loo_compare(sk_base, sk_height, sk_shannon_simple, 
                          sk_shannon, sk_shannon_full, criterion = "loo"))

plot(sk_shannon)

##______________________________________________________________________________
# 3f. CV by treatments
mssp_treat_data <- mssp_data %>% 
  mutate(treatment = if_else(treatment %in% c("Drought", "Irrigated", "Control"),
                             treatment, "Ambient"))

treatment_plot <- ggplot(mssp_treat_data, aes(x = treatment, 
                                              y = coefficient_variation, 
                                              colour = treatment,
                                              fill = treatment)) +
  geom_violin(alpha = 0.1, show.legend = F) + 
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, show.legend = F) +
  labs(x = "Treatment", y = "Spectral coefficient of variation") +
  scale_colour_viridis_d(option = "A", begin = 0.1, end = 0.9, 
                         aesthetics = c("colour", "fill")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggsave(treatment_plot, filename = "output/plots/uncalibrated/treatment_plot.jpeg",
       width = 11, height = 11, units = "cm", dpi = 800)

set.seed(666)
cv_treatment_base <- brm(
  coefficient_variation ~ 1 + (1|plot),
  data = mssp_treat_data, family = gaussian(),prior = ms_priors_simple,
  chains = 4, cores = 4, iter = 2000, warmup = 1000,
  control = list(adapt_delta = 0.9)
)

set.seed(666)
cv_treatment <- brm(
  coefficient_variation ~ 1 + treatment + (1|plot),
  data = mssp_treat_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000,
  control = list(adapt_delta = 0.9)
)

cv_treatment_base <- add_criterion(cv_treatment_base, criterion = c("loo", "waic"))
cv_treatment <- add_criterion(cv_treatment, criterion = c("loo", "waic"))

as.data.frame(loo_compare(cv_treatment_base, cv_treatment, criterion = "loo"))

plot(cv_treatment)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Best Model plots ####

##______________________________________________________________________________
# 4a. CV and Shannon index
cv_shannon_raw_plot <- ggplot(multispectral_species_diversity, 
                              aes(x = shannon, y = cv, 
                                  colour = height, group = height)) +
  geom_point(size = 4) +
  facet_wrap(~ height, scales = "free", ncol = 4) +
  scale_colour_viridis_c(begin = 0.2, end = 0.8, guide = "none") +
  labs(x = "Shannon-Weiner biodiversity index", y = "Spectral coefficient of variation") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())

ggsave(cv_shannon_raw_plot, 
       filename = "output/plots/uncalibrated/cv_shannon_raw.jpeg",
       width = 28, height = 10, units = "cm", dpi = 600)

## Set up prediction data
pdat_cv <- expand_grid(height = seq(2,8,2), 
                       shannon = seq(-1.35,3.5, length.out = 100),
                       plot = unique(mssp_data$plot))

## Posterior predictions
post_cv <- posterior_predict(cv_shannon, newdata = pdat_cv)

## Adding posterior predictions to your data
pdat_cv <- pdat_cv %>% 
  mutate(fit = colMeans(post_cv),
         upr = apply(post_cv, 2, function(x) quantile(x, c(0.9))),
         lwr = apply(post_cv, 2, function(x) quantile(x, c(0.1))),
         height_plot = paste0(height, " m")) %>% 
  group_by(height, shannon) %>% 
  summarise(height_plot = height_plot[1],
            fit = mean(fit), upr = mean(upr),
            lwr = mean(lwr))

## Posterior prediction plot
cv_shannon_post <- mssp_data %>% 
  mutate(height_plot = paste0(height, " m")) %>% 
  ggplot(aes(x = shannon, y = coefficient_variation, colour = height)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_smooth(data = pdat_cv, aes(x = shannon, ymax = upr, y = fit,
                                  ymin = lwr, fill = height, group = height_plot),
              alpha = 0.2, stat = "identity") +
  scale_colour_viridis_c(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  facet_wrap(~ height_plot, ncol = 5) +
  labs(x = "Shannon-Weiner biodiversity index", y = "Spectral coefficient of variation", tag = "a)") +
  theme_bw(base_size = 20) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())

## Posterior plots
post_plot <- tibble(brms::posterior_samples(cv_shannon)) %>% 
  dplyr::select(b_shannon, b_height) %>%
  mutate(draws = rep(1:1000, times = 4)) %>% 
  pivot_longer(-draws) %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_vline(xintercept = 0) +
  geom_density(alpha = 0.9, size = 0.3) +
  scale_fill_manual(values = c("grey21", "grey90"),
                    name = "Parameter",
                    labels = c(expression(paste(beta[Height])),
                               expression(paste(beta[Shannon])))) +
  guides(fill = guide_legend(label.hjust = 0, keywidth = 1, keyheight = 1)) +
  labs(x = "Posterior estimate", y = "Density") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank())

## Big plot with cowplot grob
cv_shannon_plot <- cv_shannon_post %>% 
  ggdraw() + draw_plot({post_plot}, 
                       x = 0.77, y = 0.63, 
                       width = 0.22, height = 0.2)
ggsave(cv_shannon_plot, 
       filename = "output/plots/uncalibrated/cv_shannon_figure.jpeg",
       width = 32, height = 14, units = "cm", dpi = 600)

##______________________________________________________________________________
# 4b. CV and Simpsons index

## Set up prediction data
pdat_cv_simp <- expand_grid(height = seq(2,8,2), 
                            simpsons = seq(-1.35,4, length.out = 100),
                            plot = unique(mssp_data$plot))

## Posterior predictions
post_cv_simp <- posterior_predict(cv_simpsons, newdata = pdat_cv_simp)

## Adding posterior predictions to your data
pdat_cv_simp <- pdat_cv_simp %>% 
  mutate(fit = colMeans(post_cv_simp),
         upr = apply(post_cv_simp, 2, function(x) quantile(x, c(0.9))),
         lwr = apply(post_cv_simp, 2, function(x) quantile(x, c(0.1))),
         height_plot = paste0(height, " m")) %>% 
  group_by(height, simpsons) %>% 
  summarise(height_plot = height_plot[1],
            fit = mean(fit), upr = mean(upr),
            lwr = mean(lwr))

## Posterior prediction plot
cv_simpsons_post <- mssp_data %>% 
  mutate(height_plot = paste0(height, " m")) %>% 
  ggplot(aes(x = simpsons, y = coefficient_variation, colour = height)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_smooth(data = pdat_cv_simp, aes(x = simpsons, ymax = upr, y = fit,
                                       ymin = lwr, fill = height, group = height_plot),
              alpha = 0.2, stat = "identity") +
  scale_colour_viridis_c(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  facet_wrap(~ height_plot, ncol = 5) +
  labs(x = "Simpson's biodiversity index", y = "Spectral coefficient of variation", tag = "b)") +
  theme_bw(base_size = 20) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())

## Just 2m
cv_simpsons_post_2m <- mssp_data %>% 
  filter(height == 2) %>% 
  ggplot(aes(x = simpsons, y = coefficient_variation)) +
  geom_point(size = 4, alpha = 0.9, colour = viridis::viridis(10)[8]) +
  geom_smooth(data = filter(pdat_cv_simp, height ==2), 
              aes(x = simpsons, ymax = upr, y = fit, 
                  ymin = lwr),
              colour = viridis::viridis(10)[8],
              fill = viridis::viridis(10)[8],
              alpha = 0.2, stat = "identity") +
  labs(x = "Biodiversity Value", y = "Variation in Spectral Radiance") +
  theme_bw(base_size = 20) +
  # ggdark::dark_theme_bw(base_size = 20) +
  theme(panel.grid = element_blank())

ggsave(cv_simpsons_post_2m, filename = "output/plots/uncalibrated/cv_simpsons_2m_light.jpeg",
       width = 15, height = 15, units = "cm", dpi = 1000)

## Posterior plots
post_plot_cv_simp <- tibble(brms::posterior_samples(cv_simpsons)) %>% 
  dplyr::select(b_simpsons, b_height) %>%
  mutate(draws = rep(1:1000, times = 4)) %>% 
  pivot_longer(-draws) %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_vline(xintercept = 0)  +
  geom_density(alpha = 0.9, size = 0.3) +
  scale_fill_manual(values = c("grey21", "grey90"),
                    name = "Parameter",
                    labels = c(expression(paste(beta[Height])),
                               expression(paste(beta[Simpsons])))) +
  guides(fill = guide_legend(label.hjust = 0, keywidth = 1, keyheight = 1)) +
  labs(x = "Posterior estimate", y = "Density") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(), 
        strip.text = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank())

## Big plot with cowplot grob
cv_simpsons_plot <- cv_simpsons_post %>% 
  ggdraw() + draw_plot({post_plot_cv_simp}, 
                       x = 0.77, y = 0.61, 
                       width = 0.22, height = 0.21)
ggsave(cv_simpsons_plot, 
       filename = "output/plots/uncalibrated/cv_simpsons_figure.jpeg",
       width = 32, height = 14, units = "cm", dpi = 600)

ggsave(cv_shannon_plot / cv_simpsons_plot, 
       filename = "output/plots/uncalibrated/cv_figure.jpeg",
       width = 32, height = 28, units = "cm", dpi = 600)

##______________________________________________________________________________
# 4c. Plot-level variance

## Hypothesis tests - ICC of plot-level variance
hypothesis(cv_shannon, 
           hypothesis = "sd_plot__Intercept^2 / (sd_plot__Intercept^2 + sigma^2) = 0",
           class = NULL)

hypothesis(cv_simpsons, 
           hypothesis = "sd_plot__Intercept^2 / (sd_plot__Intercept^2 + sigma^2) = 0",
           class = NULL)

varplot_cv_shannon <- tibble(brms::posterior_samples(cv_shannon)) %>% 
  dplyr::select(plot = sd_plot__Intercept, sigma) %>% 
  mutate(sim = 1:n()) %>% 
  pivot_longer(-sim) %>% 
  ggplot(aes(y = name, x = value)) +
  stat_halfeye(fill = viridis::viridis(1, begin = 0.5)) +
  scale_x_continuous(breaks = seq(0, 1.25, by = 0.25), limits = c(0,1.25)) +
  scale_y_discrete(labels = c(expression(paste(sigma[quadrat])),
                              expression(paste(sigma)))) +
  labs(x = "Posterior estimate", y = NULL, tag = "a)") +
  theme_ridges(font_size = 15,
               center_axis_labels = TRUE, grid = T, line_size = 0.3) 

varplot_cv_simpsons <- tibble(brms::posterior_samples(cv_simpsons)) %>% 
  dplyr::select(plot = sd_plot__Intercept, sigma) %>% 
  mutate(sim = 1:n()) %>% 
  pivot_longer(-sim) %>% 
  ggplot(aes(y = name, x = value)) +
  stat_halfeye(fill = viridis::viridis(1, begin = 0.2)) +
  scale_x_continuous(breaks = seq(0, 1.25, by = 0.25), limits = c(0,1.25)) +
  scale_y_discrete(labels = c(expression(paste(sigma[quadrat])),
                              expression(paste(sigma)))) +
  labs(x = "Posterior estimate", y = NULL, tag = "b)") +
  theme_ridges(font_size = 15,
               center_axis_labels = TRUE, grid = T, line_size = 0.3) 

ggsave(varplot_cv_shannon + varplot_cv_simpsons, 
       filename = "output/plots/uncalibrated/cv_sigma_plots.jpeg",
       width = 20, height = 8, units = "cm", dpi = 600)


##______________________________________________________________________________
# 4d. Band-level Simpsons with CV

## Set up prediction data
p_cv_simp_band <- expand_grid(height = 2, 
                              simpsons = seq(-1.35,4, length.out = 50),
                              band = unique(mssp_bands_data$band),
                              plot = unique(mssp_bands_data$plot))

## Posterior predictions
post_cv_simp_band <- posterior_predict(cvb_band_full, newdata = p_cv_simp_band)

## Adding posterior predictions to your data
pdat_cv_simp_band <- p_cv_simp_band %>% 
  mutate(fit = colMeans(post_cv_simp_band),
         upr = apply(post_cv_simp_band, 2, function(x) quantile(x, c(0.9))),
         lwr = apply(post_cv_simp_band, 2, function(x) quantile(x, c(0.1))),
         height_plot = paste0(height, "m")) %>% 
  group_by(height, band, simpsons) %>% 
  summarise(height_plot = height_plot[1],
            fit = mean(fit), upr = mean(upr),
            lwr = mean(lwr)) %>% 
  ungroup() %>% 
  mutate(band = factor(band, levels = c("B", "G", "R", "RE", "NIR")))

# Plot
cv_simpsons_band_post <- mssp_bands_data %>% 
  mutate(height_plot = paste0(height, "m"),
         band = factor(band, levels = c("B", "G", "R", "RE", "NIR"))) %>% 
  ggplot(aes(x = simpsons, y = coefficient_variation, colour = band)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(data = pdat_cv_simp_band, aes(x = simpsons, ymax = upr, y = fit,
                                            ymin = lwr, fill = band, group = band),
              alpha = 0.2, stat = "identity") +
  scale_colour_manual(values = c("blue", "green", "red", "brown", "grey69"),
                      aesthetics = c("colour", "fill"), guide = "none") +
  facet_grid(~ band) +
  labs(x = "Simpson's biodiversity index", y = "Spectral coefficient of variation", 
       tag = "b)") +
  theme_bw(base_size = 22) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())

##______________________________________________________________________________
# 4e. Band-level Shannon with CV

## Set up prediction data
p_cv_shan_band <- expand_grid(height = 2, 
                              shannon = seq(-1.35,4, length.out = 50),
                              band = unique(mssp_bands_data$band),
                              plot = unique(mssp_bands_data$plot))

## Posterior predictions
post_cv_shan_band <- posterior_predict(cvb_band_full_shannon, newdata = p_cv_shan_band)

## Adding posterior predictions to your data
pdat_cv_shan_band <- p_cv_shan_band %>% 
  mutate(fit = colMeans(post_cv_shan_band),
         upr = apply(post_cv_shan_band, 2, function(x) quantile(x, c(0.9))),
         lwr = apply(post_cv_shan_band, 2, function(x) quantile(x, c(0.1))),
         height_plot = paste0(height, "m")) %>% 
  group_by(height, band, shannon) %>% 
  summarise(height_plot = height_plot[1],
            fit = mean(fit), upr = mean(upr),
            lwr = mean(lwr)) %>% 
  ungroup() %>% 
  mutate(band = factor(band, levels = c("B", "G", "R", "RE", "NIR")))

# Plot
cv_shannon_band_post <- mssp_bands_data %>% 
  mutate(height_plot = paste0(height, "m"),
         band = factor(band, levels = c("B", "G", "R", "RE", "NIR"))) %>% 
  ggplot(aes(x = shannon, y = coefficient_variation, colour = band)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(data = pdat_cv_shan_band, aes(x = shannon, ymax = upr, y = fit,
                                            ymin = lwr, fill = band, group = band),
              alpha = 0.2, stat = "identity") +
  scale_colour_manual(values = c("blue", "green", "red", "brown", "grey69"),
                      aesthetics = c("colour", "fill"), guide = "none") +
  facet_grid(~ band) +
  labs(x = "Shannon-Weiner biodiversity index", 
       y = "Spectral coefficient of variation", tag = "a)") +
  theme_bw(base_size = 22) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())

ggsave(cv_shannon_band_post/ cv_simpsons_band_post,
       filename = "output/plots/uncalibrated/cv_band_figure.jpeg",
       width = 30, height = 26, units = "cm", dpi = 1500)

cvb_band_full

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Model selection tables ####

##______________________________________________________________________________
# 5a. CV and Shannon index

as.data.frame(loo_compare(cv_base, cv_height, cv_shannon_simple, 
                          cv_shannon, cv_shannon_full, criterion = "loo")) %>% 
  mutate(predictor_terms = c("height", "shannon + height", "shannon + height + shannon:height",
                             "shannon", "base model")) %>% 
  dplyr::select(predictor_terms, elpd_loo, se_elpd_loo, elpd_diff, se_diff, looic) %>% 
  flextable(cwidth = 1.2) %>% 
  set_header_labels(predictor_terms = "Model predictors",
                    elpd_loo = "LOO elpd",
                    se_elpd_loo = "LOO elpd error",
                    elpd_diff = "elpd difference",
                    se_diff = "elpd error difference",
                    looic = "LOO information criterion") %>% 
  colformat_double(digits = 2) %>% 
  save_as_docx(path = "output/plots/uncalibrated/table1_cv_shannon.docx")

##______________________________________________________________________________
# 5b. CV and Simpsons index

as.data.frame(loo_compare(cv_base, cv_height, cv_simpsons_simple, 
                          cv_simpsons, cv_simpsons_full, criterion = "loo")) %>% 
  mutate(predictor_terms = c("simpsons + height + simpsons:height", "simpsons + height", "height",
                             "simpsons", "base model")) %>% 
  dplyr::select(predictor_terms, elpd_loo, se_elpd_loo, elpd_diff, se_diff, looic) %>% 
  flextable(cwidth = 1.2) %>% 
  set_header_labels(predictor_terms = "Model predictors",
                    elpd_loo = "LOO elpd",
                    se_elpd_loo = "LOO elpd error",
                    elpd_diff = "elpd difference",
                    se_diff = "elpd error difference",
                    looic = "LOO information criterion") %>% 
  colformat_double(digits = 2) %>% 
  save_as_docx(path = "output/plots/uncalibrated/table2_cv_simpsons.docx")

##______________________________________________________________________________
# 5c. CV Shannon band
as.data.frame(loo_compare(cvb_base, cvb_shannon, cvb_band_shannon, cvb_band_full_shannon, criterion = "loo")) %>% 
  mutate(predictor_terms = c("shannon + band + shannon:band + height", "shannon + height + band", 
                             "shannon + height", "base model")) %>% 
  dplyr::select(predictor_terms, elpd_loo, se_elpd_loo, elpd_diff, se_diff, looic) %>% 
  flextable(cwidth = 1.2) %>% 
  set_header_labels(predictor_terms = "Model predictors",
                    elpd_loo = "LOO elpd",
                    se_elpd_loo = "LOO elpd error",
                    elpd_diff = "elpd difference",
                    se_diff = "elpd error difference",
                    looic = "LOO information criterion") %>% 
  colformat_double(digits = 2) %>% 
  save_as_docx(path ="output/plots/uncalibrated/table3_cv_band_shannon.docx")

##______________________________________________________________________________
# 5d. CV Simpsons band

as.data.frame(loo_compare(cvb_base, cvb_simpsons, cvb_band, cvb_band_full, criterion = "loo")) %>% 
  mutate(predictor_terms = c("simpsons + band + simpsons:band + height", "simpsons + height + band", 
                             "simpsons + height", "base model")) %>% 
  dplyr::select(predictor_terms, elpd_loo, se_elpd_loo, elpd_diff, se_diff, looic) %>% 
  flextable(cwidth = 1.2) %>% 
  set_header_labels(predictor_terms = "Model predictors",
                    elpd_loo = "LOO elpd",
                    se_elpd_loo = "LOO elpd error",
                    elpd_diff = "elpd difference",
                    se_diff = "elpd error difference",
                    looic = "LOO information criterion") %>% 
  colformat_double(digits = 2) %>% 
  save_as_docx(path ="output/plots/uncalibrated/table4_cv_band_simpsons.docx")

##______________________________________________________________________________
# 5e. Skewness and Biomass

as.data.frame(loo_compare(sk_base, sk_height, sk_biomass_simple, 
                          sk_biomass, sk_biomass_full, criterion = "loo")) %>% 
  mutate(predictor_terms = c("height", "biomass + height", "biomass + height + biomass:height", 
                             "base model", "biomass")) %>% 
  dplyr::select(predictor_terms, elpd_loo, se_elpd_loo, elpd_diff, se_diff, looic) %>% 
  flextable(cwidth = 1.2) %>% 
  set_header_labels(predictor_terms = "Model predictors",
                    elpd_loo = "LOO elpd",
                    se_elpd_loo = "LOO elpd error",
                    elpd_diff = "elpd difference",
                    se_diff = "elpd error difference",
                    looic = "LOO information criterion") %>% 
  colformat_double(digits = 2) %>% 
  save_as_docx(path = "output/plots/uncalibrated/table5_sk_biomass.docx")
