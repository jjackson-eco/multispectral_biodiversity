#########################################################
##                                                     ##
##  RainDrop: Multispectral indicators of Biodiversity ##
##                                                     ##
##         Biodiversity vs spectral diversity          ## 
##                                                     ##
##               DRAGNet uncalibrated                  ##
##                                                     ##
##                  Oct 19th 2021                      ##
##                                                     ##
#########################################################
# updated 2022-09-02
rm(list = ls())

# lib paths for uni computer
source("../lib_path_uni_computer.R")

library(tidyverse)
library(patchwork)
library(brms)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Loading data ####

load("data/biodiversity_data_dragnet_early.RData", verbose = TRUE)
load("data/multispec_data_dragnet_uncalibrated.RData", verbose = TRUE)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Merging data ####

# Biomass summary
biomass_sum <- biomass_dragnet %>% 
  group_by(block, treatment) %>% 
  mutate(tot_biomass = sum(biomass)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = functional_group, values_from = biomass)

# Merge for overall plot/height summary
multispectral_species_diversity <- multispec_data_dragnet_uncalibrated %>% 
  left_join(x = ., y = percent_cover_dragnet, by = c("block", "treatment")) %>% 
  left_join(x = ., y = biomass_sum, by = c("block", "treatment")) %>% 
  group_by(block, treatment, height, replicate) %>% 
  summarise(cv = mean(cv), sd = mean(sd), 
            skew = mean(skew), kurtosis = mean(kurtosis),
            richness = richness[1], simpsons = simpsons[1], 
            shannon = shannon[1], tot_biomass = tot_biomass[1],
            graminoids = Graminoids[1], legumes = Legumes[1], 
            forbs = Forbs[1]) %>% 
  ungroup() %>% 
  filter(height < 10) # remove 10m as it is only one image

# Merge for full data across bands 
multispectral_species_diversity_allbands <- multispec_data_dragnet_uncalibrated %>% 
  left_join(x = ., y = percent_cover_dragnet, by = c("block", "treatment")) %>% 
  left_join(x = ., y = biomass_sum, by = c("block", "treatment")) %>% 
  group_by(block, treatment, height, replicate, band) %>% 
  summarise(band_nm = band_nm[1], cv = mean(cv), sd = mean(sd), 
            skew = mean(skew), kurtosis = mean(kurtosis),
            richness = richness[1], simpsons = simpsons[1], 
            shannon = shannon[1], tot_biomass = tot_biomass[1],
            graminoids = Graminoids[1], legumes = Legumes[1], 
            forbs = Forbs[1]) %>% 
  ungroup() %>% 
  filter(height < 10) # remove 10m as it is only one image

# # CSV data, for future work
# multispectral_species_diversity_full <- multispec_data_dragnet_uncalibrated %>% 
#   left_join(x = ., y = percent_cover, by = c("block", "treatment")) %>% 
#   left_join(x = ., y = biomass_sum, by = c("block", "treatment")) %>% 
#   dplyr::select(file, id, image_date = rec_date, experiment, block, 
#                 treatment, height, band, band_nm, replicate, cv, sd, skew, 
#                 kurtosis, date_cover = date.x, richness, simpsons, shannon,
#                 date_biomass = date.y, tot_biomass, graminoids = Graminoids,
#                 forbs = Forbs, legumes = Legumes, bryophytes = Bryophytes, 
#                 woody = Woody)
# 
# write_csv(multispectral_species_diversity_full, 
#           file = "data/multispec_biodiversity_dragnet_uncalibrated.csv")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Plotting with all moments and proxys ####

## 3a. summary data
ms_sp_long <- multispectral_species_diversity %>% 
  dplyr::select(block, treatment, height, cv, sd, skew, kurtosis,
                richness, simpsons, shannon, tot_biomass) %>% 
  pivot_longer(cols = c("cv", "sd", "skew", "kurtosis"), 
                        names_to = "moment", values_to = "moment_value") %>% 
  pivot_longer(cols = c("richness", "simpsons", "shannon", "tot_biomass"),
               names_to = "biodiversity_indicator", 
               values_to = "biodiversity_value")


all_plot <- ggplot(ms_sp_long, aes(x = biodiversity_value, 
                       y = moment_value, colour = height)) +
  geom_point() +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  facet_wrap(~ moment + biodiversity_indicator, scales = "free") +
  scale_colour_viridis_c() +
  labs(x = "Biodiversity indicator value", y = "Multispectral moment value",
       colour = "Height (m)") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())


ggsave(all_plot, filename = "output/plots/dragnet_uncalibrated/all_moments_indicator_plot.jpeg",
       height = 25, width = 30, units = "cm", dpi = 600)

## 3b. separate bands - this is a bit of a mess
ms_sp_bands_long <- multispectral_species_diversity_allbands %>%
  dplyr::select(block, treatment, height, band, cv, sd, skew, kurtosis,
                richness, simpsons, shannon, tot_biomass) %>%
  pivot_longer(cols = c("cv", "sd", "skew", "kurtosis"),
               names_to = "moment", values_to = "moment_value") %>%
  pivot_longer(cols = c("richness", "simpsons", "shannon", "tot_biomass"),
               names_to = "biodiversity_indicator",
               values_to = "biodiversity_value")
# 
# all_bands_plot <- ggplot(ms_sp_bands_long, aes(x = biodiversity_value, 
#                                    y = moment_value, colour = band)) +
#   geom_point() +
#   geom_smooth(aes(group = interaction(height, band)), method = "lm", se = F) +
#   facet_wrap(band ~ moment + biodiversity_indicator, scales = "free") +
#   scale_colour_viridis_d() +
#   labs(x = "Biodiversity indicator value", y = "Multispectral moment value",
#        colour = "Height (m)") +
#   theme_bw(base_size = 13) +
#   theme(panel.grid = element_blank(),
#         strip.background = element_blank())
# 
# ggsave(all_bands_plot, filename = "output/plots/dragnet_uncalibrated/all_moments_indicator_plot_all_bands.jpeg",
#        height = 40, width = 40, units = "cm", dpi = 600)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Explore relationship with separate bands ####

# Going in more depth and looking at separate bands as a plot
cv_shannon_bands_long <- ms_sp_bands_long %>% 
  filter(moment == "cv" & 
           biodiversity_indicator %in% c("shannon", "richness") == T) %>% 
  mutate(band = factor(band, levels = c("B", "G", "R", "RE", "NIR")))

cv_shannon_bands_plot <- cv_shannon_bands_long %>% 
  ggplot(aes(x = biodiversity_value, y = moment_value, 
             colour = band)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(aes(group = interaction(height, band)), method = "lm", se = F) +
  facet_wrap( ~ band + biodiversity_indicator, scales = "free",
             ncol = 5) +
  scale_colour_manual(values = c("blue", "green", "red", "darkred", "grey")) +
  labs(x = "Biodiversity indicator value", 
       y = "Mean coefficient of variation") +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

ggsave(cv_shannon_bands_plot, filename = "output/plots/dragnet_uncalibrated/cv_sp_band_plot.jpeg",
       height = 18, width = 30, units = "cm", dpi = 600)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Model selection ####

## Data convert for bayesian models
mssp_data <- multispectral_species_diversity %>% 
  dplyr::select(block, treatment, height, replicate, shannon, tot_biomass, 
                skew, coefficient_variation = cv) %>% 
  mutate(height_s = as.numeric(scale(height)),
         shannon = as.numeric(scale(shannon)),
         tot_biomass = as.numeric(scale(tot_biomass)),
         skew = as.numeric(scale(skew)),
         coefficient_variation = as.numeric(scale(coefficient_variation)))

## Priors
ms_priors_simple <- c(prior(normal(0, 1), class = Intercept))
ms_priors <- c(prior(normal(0, 1), class = Intercept),
               prior(normal(0,1), class = b))

##______________________________________________________________________________
# 5a. CV and Shannon- core question

## Base model
set.seed(666)
cv_base <- brm(
  coefficient_variation ~ 1,
  data = mssp_data, family = gaussian(),prior = ms_priors_simple,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV height
set.seed(666)
cv_height <- brm(
  coefficient_variation ~ 1 + height,
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV shannon simple
set.seed(666)
cv_shannon_simple <- brm(
  coefficient_variation ~ 1 + shannon,
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV shannon
set.seed(666)
cv_shannon <- brm(
  coefficient_variation ~ 1 + shannon + height,
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV shannon full
set.seed(666)
cv_shannon_full <- brm(
  coefficient_variation ~ 1 + shannon*height,
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
# 5b. CV and Shannon by band

## Data convert for bayesian models
mssp_bands_data <- multispectral_species_diversity_allbands %>% 
  dplyr::select(block, treatment, band, band_nm, height, 
                replicate, shannon,tot_biomass, 
                skew, coefficient_variation = cv) %>% 
  mutate(height_s = as.numeric(scale(height)),
         shannon = as.numeric(scale(shannon)),
         tot_biomass = as.numeric(scale(tot_biomass)),
         skew = as.numeric(scale(skew)),
         coefficient_variation = as.numeric(scale(coefficient_variation)))

## CV with all data base
set.seed(666)
cvb_base <- brm(
  coefficient_variation ~ 1,
  data = mssp_bands_data, family = gaussian(), prior = ms_priors_simple,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV shannon with all data
set.seed(666)
cvb_shannon <- brm(
  coefficient_variation ~ 1 + shannon + height,
  data = mssp_bands_data, family = gaussian(), prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV band
set.seed(666)
cvb_band <- brm(
  coefficient_variation ~ 1 + shannon + height + band,
  data = mssp_bands_data, family = gaussian(), prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV band full
set.seed(666)
cvb_band_full <- brm(
  coefficient_variation ~ 1 + shannon*band + height,
  data = mssp_bands_data, family = gaussian(), prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

cvb_base <- add_criterion(cvb_base, criterion = c("loo","waic"))
cvb_shannon <- add_criterion(cvb_shannon, criterion = c("loo","waic"))
cvb_band <- add_criterion(cvb_band, criterion = c("loo","waic"))
cvb_band_full<- add_criterion(cvb_band_full, criterion = c("loo","waic"))

as.data.frame(loo_compare(cvb_base, cvb_shannon, cvb_band, cvb_band_full, criterion = "loo"))

plot(cvb_band_full)


##______________________________________________________________________________
# 5c. Skewness and Biomass

## Base model
set.seed(666)
sk_base <- brm(
  skew ~ 1,
  data = mssp_data, family = gaussian(),prior = ms_priors_simple,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Height only
set.seed(666)
sk_height <- brm(
  skew ~ 1 + height,
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Biomass simple
set.seed(666)
sk_biomass_simple <- brm(
  skew ~ 1 + tot_biomass,
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Biomass 
set.seed(666)
sk_biomass <- brm(
  skew ~ 1 + tot_biomass + height,
  data = mssp_data, family = gaussian(),prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Biomass full
set.seed(666)
sk_biomass_full <- brm(
  skew ~ 1 + tot_biomass*height,
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
# 5d. Skewness and Biomass by band

## Skewness Biomass base
set.seed(666)
skb_base <- brm(
  skew ~ 1,
  data = mssp_bands_data, family = gaussian(), prior = ms_priors_simple,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## Skewness Biomass with all data
set.seed(666)
skb_biomass <- brm(
  skew ~ 1 + tot_biomass + height,
  data = mssp_bands_data, family = gaussian(), prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV band
set.seed(666)
skb_band <- brm(
  skew ~ 1 + tot_biomass + height + band,
  data = mssp_bands_data, family = gaussian(), prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

## CV band full
set.seed(666)
skb_band_full <- brm(
  skew ~ 1 + tot_biomass*band + height,
  data = mssp_bands_data, family = gaussian(), prior = ms_priors,
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

skb_base <- add_criterion(skb_base, criterion = c("loo","waic"))
skb_biomass <- add_criterion(skb_biomass, criterion = c("loo","waic"))
skb_band <- add_criterion(skb_band, criterion = c("loo","waic"))
skb_band_full<- add_criterion(skb_band_full, criterion = c("loo","waic"))

as.data.frame(loo_compare(skb_base, skb_biomass, skb_band, skb_band_full, criterion = "loo"))

plot(skb_band)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Best Model plots ####

## Raw shannon - cv plots
cv_shannon_raw_plot <- ggplot(multispectral_species_diversity, 
       aes(x = shannon, y = cv, 
           colour = height, group = height)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~ height, scales = "free", ncol = 4) +
  scale_colour_viridis_c(begin = 0.2, end = 0.8, guide = "none") +
  labs(x = "Shannon biodiversity index", y = "Spectral coefficient of variation") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())

ggsave(cv_shannon_raw_plot, 
       filename = "output/plots/dragnet_uncalibrated/cv_shannon_raw.jpeg",
       width = 28, height = 10, units = "cm", dpi = 600)


pdat_cv <- expand_grid(height = seq(2,8,2), 
                       shannon = seq(-1.35,2.35, length.out = 100))

post_cv <- posterior_predict(cv_shannon, newdata = pdat_cv)

pdat_cv <- pdat_cv %>% 
  mutate(fit = colSums(post_cv))



cv_shannon_plot <- mssp_data %>% 
  ggplot(aes(x = shannon, y = coefficient_variation, colour = height)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_smooth(stat = "identity", aes(y = fit, group = height), se = F) +
  scale_colour_viridis_c(guide = "none") +
  facet_wrap(~height, ncol = 5) +
  labs(x = "Shannon-Weiner biodiviersity index", 
       y = "Spectral diversity (average coefficient of variation)") +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank())

ggsave(cv_shannon_plot, filename = "output/plots/dragnet_uncalibrated/cv_shannon_plot.jpeg",
       height = 15, width = 25, units = "cm", dpi = 600)
