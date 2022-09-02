#########################################################
##                                                     ##
##  RainDrop: Multispectral indicators of Biodiversity ##
##                                                     ##
##         Biodiversity vs spectral diversity          ## 
##                                                     ##
##            Supplementary exploration                ##
##                                                     ##
##                  8th Feb 2022                       ##
##                                                     ##
#########################################################
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
library(rnaturalearth)
library(viridis)

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

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Plotting with all moments and proxys ####

## 3a. CV 
rawplot_a <- ggplot(multispectral_species_diversity,
                              aes(x = richness, y = cv, 
                                  colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Species richness", y = "Spectral coefficient of variation", tag = "a)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

rawplot_b <- ggplot(multispectral_species_diversity,
                    aes(x = shannon, y = cv, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Shannon-Weiner index", 
       y = "Spectral coefficient of variation", tag = "b)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

rawplot_c <- ggplot(multispectral_species_diversity,
                    aes(x = simpsons, y = cv, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Simpson's index", 
       y = "Spectral coefficient of variation", tag = "c)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

rawplot_d <- ggplot(multispectral_species_diversity,
                    aes(x = tot_biomass, y = cv, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Total above ground biomass (g)", 
       y = "Spectral coefficient of variation", tag = "d)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

## 3b. SD
rawplot_e <- ggplot(multispectral_species_diversity,
                    aes(x = richness, y = sd, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Species richness", y = "Spectral standard deviation", tag = "e)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

rawplot_f <- ggplot(multispectral_species_diversity,
                    aes(x = shannon, y = sd, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Shannon-Weiner index", 
       y = "Spectral standard deviation", tag = "f)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

rawplot_g <- ggplot(multispectral_species_diversity,
                    aes(x = simpsons, y = sd, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Simpson's index", 
       y = "Spectral standard deviation", tag = "g)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

rawplot_h <- ggplot(multispectral_species_diversity,
                    aes(x = tot_biomass, y = sd, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Total above ground biomass (g)", 
       y = "Spectral standard deviation", tag = "h)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

## 3c. Skewness
rawplot_i <- ggplot(multispectral_species_diversity,
                    aes(x = richness, y = skew, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Species richness", y = "Spectral skewness", tag = "i)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

rawplot_j <- ggplot(multispectral_species_diversity,
                    aes(x = shannon, y = skew, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Shannon-Weiner index", 
       y = "Spectral skewness", tag = "j)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

rawplot_k <- ggplot(multispectral_species_diversity,
                    aes(x = simpsons, y = skew, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Simpson's index", 
       y = "Spectral skewness", tag = "k)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

rawplot_l <- ggplot(multispectral_species_diversity,
                    aes(x = tot_biomass, y = skew, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Total above ground biomass (g)", 
       y = "Spectral skewness", tag = "l)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

## 3d. Kurtosis
rawplot_m <- ggplot(multispectral_species_diversity,
                    aes(x = richness, y = kurtosis, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Species richness", y = "Spectral kurtosis", tag = "m)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

rawplot_n <- ggplot(multispectral_species_diversity,
                    aes(x = shannon, y = kurtosis, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Shannon-Weiner index", 
       y = "Spectral kurtosis", tag = "n)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

rawplot_o <- ggplot(multispectral_species_diversity,
                    aes(x = simpsons, y = kurtosis, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Simpson's index", 
       y = "Spectral kurtosis", tag = "o)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

rawplot_p <- ggplot(multispectral_species_diversity,
                    aes(x = tot_biomass, y = kurtosis, 
                        colour = as.factor(height))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = height), method = "lm", se = F) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  labs(x = "Total above ground biomass (g)", 
       y = "Spectral kurtosis", tag = "p)") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

## 3e. Height plot
rawplot_height <- tibble(height = seq(2,8, by = 2),
       height_lab = paste0(seq(2,8, by = 2), " m"),
       height_mid = c(1.5, 3.5, 5.5, 7.5),
       x  = seq(1,7, by = 2),
       x_end =  seq(2,8, by = 2)) %>% 
  ggplot(aes(colour = as.factor(height))) +
  geom_segment(aes(y = 1, yend = 1, x = x, xend = x_end),
               size = 5) +
  geom_text(aes(x = height_mid, y = 1.02, label = height_lab),
            colour = "black", size = 6) +
  scale_colour_viridis_d(guide = "none", begin = 0.2, end = 0.8, 
                         aesthetics = c("colour", "fill")) +
  coord_cartesian(ylim = c(0.96,1.04)) +
  theme_void()

## overall plot

layout <- "
ABCD
EFGH
IJKL
MNOP
QQQQ
"

rawplot <- rawplot_a + rawplot_b + rawplot_c + rawplot_d +
  rawplot_e + rawplot_f + rawplot_g + rawplot_h +
  rawplot_i + rawplot_j + rawplot_k + rawplot_l +
  rawplot_m + rawplot_n + rawplot_o + rawplot_p +
  rawplot_height + plot_layout(design = layout, heights = c(2,2,2,2,1))

ggsave(rawplot, 
       filename = "output/plots/uncalibrated/all_moments_indicator_plot.jpeg",
       width = 40, height = 40, units = "cm", dpi = 1000)

## 3e. separate bands - this is a bit of a mess
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

## 3f. Looking at the odd CV images
multispec_data_uncalibrated %>% 
  group_by(block, treatment, experiment, height, replicate) %>% 
  summarise(mn_cv = mean(cv), file = file[1]) %>% 
  filter(mn_cv < 20) 

## 3g. Treatment biodiversity comparison plots
treatment_summary_data <- multispectral_species_diversity %>% 
  filter(treatment == "Ambient" | experiment == "DRAGNET") %>% 
  group_by(block, treatment, experiment) %>% 
  summarise(richness = richness[1],
            simpsons = simpsons[1], 
            shannon = shannon[1], 
            tot_biomass = tot_biomass[1]) 

amb_rich <- ggplot(treatment_summary_data,
                   aes(x = experiment, y = richness, colour = experiment, 
                       fill = experiment)) +
  geom_violin(alpha = 0.2, colour = "white", show.legend = F) +
  geom_jitter(size = 4, width = 0.1, show.legend = F) +
  stat_summary(geom = "point", fun = mean, shape = 3, size = 6,
               show.legend = F) +
  scale_x_discrete(labels = c("DRAGNet", "DroughtNet")) +
  labs(x = "Ambient control experiment", y = "Species richness") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

amb_shan <- ggplot(treatment_summary_data,
                   aes(x = experiment, y = shannon, colour = experiment, 
                       fill = experiment)) +
  geom_violin(alpha = 0.2, colour = "white", show.legend = F) +
  geom_jitter(size = 4, width = 0.1, show.legend = F) +
  stat_summary(geom = "point", fun = mean, shape = 3, size = 6,
               show.legend = F) +
  scale_x_discrete(labels = c("DRAGNet", "DroughtNet")) +
  labs(x = "Ambient control experiment", y = "Shannon-Weiner biodiversity index") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

amb_simp <- ggplot(treatment_summary_data,
       aes(x = experiment, y = simpsons, colour = experiment, 
           fill = experiment)) +
  geom_violin(alpha = 0.2, colour = "white", show.legend = F) +
  geom_jitter(size = 4, width = 0.1, show.legend = F) +
  stat_summary(geom = "point", fun = mean, shape = 3, size = 6,
               show.legend = F) +
  scale_x_discrete(labels = c("DRAGNet", "DroughtNet")) +
  labs(x = "Ambient control experiment", y = "Simpson's biodiversity index") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

amb_biom <- ggplot(treatment_summary_data,
                   aes(x = experiment, y = tot_biomass, colour = experiment, 
                       fill = experiment)) +
  geom_violin(alpha = 0.2, colour = "white", show.legend = F) +
  geom_jitter(size = 4, width = 0.1, show.legend = F) +
  stat_summary(geom = "point", fun = mean, shape = 3, size = 6,
               show.legend = F) +
  scale_x_discrete(labels = c("DRAGNet", "DroughtNet")) +
  labs(x = "Ambient control experiment", y = expression(paste("Total biomass (g/", m^2, ")"))) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

ggsave((amb_rich + amb_shan) /
         (amb_simp + amb_biom), filename = "output/plots/uncalibrated/treatment_ambient_plot.jpeg",
       height = 30, width = 30, units = "cm", dpi = 600)

# stats
t.test(richness ~ experiment, data = treatment_summary_data)
t.test(shannon ~ experiment, data = treatment_summary_data)
t.test(log(simpsons) ~ experiment, data = treatment_summary_data)
t.test(log(tot_biomass) ~ experiment, data = treatment_summary_data)

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

ggsave(cv_shannon_bands_plot, filename = "output/plots/uncalibrated/cv_sp_band_plot.jpeg",
       height = 18, width = 30, units = "cm", dpi = 600)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Map of Wytham Woods ####
countries <- ne_countries(scale = 50,returnclass = "sf")

wytham_woods_map <- ggplot(countries) +
  geom_sf(size = 0.1) +
  geom_point(aes(x = -1.333083, y = 51.771333),
             colour = viridis(10)[6], size = 5) +
  annotate("text",x = -1, y = 51.771333, label = "Wytham\nWoods", hjust = 0, size = 12) +
  coord_sf(xlim = c(-8, 8), ylim = c(55, 50)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

ggsave(wytham_woods_map, filename = "output/plots/uncalibrated/wytham_woods_map.jpeg",
       width = 20, height = 15, units = "cm", dpi = 1000)

