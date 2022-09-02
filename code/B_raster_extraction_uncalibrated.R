#########################################################
##                                                     ##
##  RainDrop: Multispectral indicators of Biodiversity ##
##                                                     ##
##            TIFF raster data extraction              ## 
##                                                     ##
##                Uncalibrated images                  ##
##                                                     ##
##                  Oct 14th 2021                      ##
##                                                     ##
#########################################################
# updated 2022-09-02
rm(list = ls())

# lib paths for uni computer
source("../lib_path_uni_computer.R")

library(tidyverse)
library(patchwork)
library(raster)
library(rasterVis)
library(psych)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Setting up iteration info and meta-data ####

# all image files incl RGB and reference images
image_files_drag <- list.files("../cropped_image_data/un_corrected_20211014/")
image_files_drought <- list.files("../cropped_image_data/un_corrected_20211124/")

image_files <- c(image_files_drag, image_files_drought)

# meta-data extraction
image_meta_data <- bind_rows(lapply(image_files, function(x){
  
  crr_image = gsub(".tif", "", x)
  
  # parsing out the data from file name
  id = unlist(strsplit(crr_image, split = "_"))[1]
  rec_date = as.Date(unlist(strsplit(crr_image, split = "_"))[2], 
                     format = "%Y%m%d")
  experiment = unlist(strsplit(crr_image, split = "_"))[3]
  source = ifelse(experiment == "DRAGNET", 
                  "un_corrected_20211014/", "un_corrected_20211124/")
  block = unlist(strsplit(crr_image, split = "_"))[4]
  treatment = unlist(strsplit(crr_image, split = "_"))[5]
  height = as.numeric(unlist(strsplit(crr_image, split = "_"))[6])
  band = unlist(strsplit(crr_image, split = "_"))[7]
  replicate = as.numeric(unlist(strsplit(crr_image, split = "_"))[8])
  
  return(data.frame(source = source, file = x, id = id, rec_date = rec_date,
                    experiment = experiment, block = block, 
                    treatment = treatment, height = height, band = band, 
                    replicate = replicate))
  
}))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Extracting distribution summaries from raster data ####

# pulling only multispectral images and no reference images
image_files_ms <- image_files[-(grep("REF|RGB", image_files))]

# New columns to be added to the image metadata
ms_image_data <- image_meta_data %>% 
  # adding data for the treatments
  mutate(treatment = case_when(
    treatment == "N" & experiment == "DRAGNET" ~ "NPK",
    treatment == "D" & experiment == "DRAGNET" ~ "Disturbance",
    treatment == "ND" & experiment == "DRAGNET" ~ "NPK + Disturbance",
    treatment == "NC" & experiment == "DRAGNET" ~ "NPK Cessation",
    treatment == "D" & experiment == "DROUGHTNET" ~ "Drought",
    treatment == "A" & experiment == "DROUGHTNET" ~ "Ambient",
    treatment == "C" & experiment == "DROUGHTNET" ~ "Control",
    treatment == "I" & experiment == "DROUGHTNET" ~ "Irrigated")) %>% 
  filter(file %in% image_files_ms) %>% 
  mutate(cv = 0, sd = 0, skew = 0, kurtosis = 0, pixel_dim = 0)

# Iterate through and compute summary stats for the rasters
for(i in 1:nrow(ms_image_data)){
  
  # current data/file
  crr_dat = ms_image_data[i,]
  crr_img = crr_dat$file
  crr_source = crr_dat$source
  
  # load the current cropped image as a raster
  c_raster = raster(paste0("../cropped_image_data/", crr_source, crr_img))
  
  # compute the summary statistics
  ms_image_data[i,"cv"] <- cellStats(c_raster, cv)
  ms_image_data[i,"sd"] <- cellStats(c_raster, sd)
  ms_image_data[i,"skew"] <- cellStats(c_raster, skew)
  ms_image_data[i,"kurtosis"] <- cellStats(c_raster, kurtosi)
  ms_image_data[i,"pixel_dim"] <- dim(c_raster)[1]
  
  cat('\r', "Your job is ", round(i/nrow(ms_image_data), 3)*100, " % Complete              ", sep = "")
  
}

multispec_data_uncalibrated <- ms_image_data %>% 
  # Adding information on the wavelength
  mutate(band_nm = case_when(
    band == "B" ~ 450,
    band == "G" ~ 560,
    band == "R" ~ 650,
    band == "RE" ~ 730,
    band == "NIR" ~ 840
  ))

save(multispec_data_uncalibrated, 
     file = "data/multispec_data_uncalibrated.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Simple exploration of the image data ####

#__________________________________________
### 3a. Overall distributions
cv_plot <- ggplot(multispec_data_uncalibrated, aes(x = cv)) +
  geom_histogram(fill = "darkgrey", colour = "black", 
                 size = 0.08,bins = 50) +
  labs(x = "Coefficient of Variation", y = "Frequency") +
  theme_bw(base_size = 12)

sd_plot <- ggplot(multispec_data_uncalibrated, aes(x = sd)) +
  geom_histogram(fill = "darkgrey", colour = "black", 
                 size = 0.08,bins = 50) +
  labs(x = "Standard deviation", y = "Frequency") +
  theme_bw(base_size = 12)

sk_plot <- ggplot(multispec_data_uncalibrated, aes(x = skew)) +
  geom_histogram(fill = "darkgrey", colour = "black", 
                 size = 0.08,bins = 50) +
  labs(x = "Skewness", y = "Frequency") +
  theme_bw(base_size = 12)

kt_plot <- ggplot(multispec_data_uncalibrated, aes(x = kurtosis)) +
  geom_histogram(fill = "darkgrey", colour = "black", 
                 size = 0.08,bins = 50) +
  labs(x = "Kurtosis", y = "Frequency") +
  theme_bw(base_size = 12)

(cv_plot + sd_plot) / (sk_plot + kt_plot)

ggsave((cv_plot + sd_plot) / (sk_plot + kt_plot),
       filename = "output/plots/uncalibrated/distribution_summary.jpeg",
       width = 16, height = 10, units = "cm", dpi = 600)

## finding the weird Kurtosis value - makes complete sense as it is very dark
multispec_data_uncalibrated %>% 
  filter(kurtosis > 50) %>% 
  pull(file)

#__________________________________________
### 3b. Relationships with height

cv_height <- ggplot(multispec_data_uncalibrated, aes(x = height, y = cv)) +
  geom_jitter(size = 1, alpha = 0.4,width = 0.2) +
  geom_smooth(se = F, colour = "black", size = 0.5, linetype = "dashed") +
  labs(x = "Height", y = "Coefficient of Variation") +
  theme_bw(base_size = 12)

sd_height <- ggplot(multispec_data_uncalibrated, aes(x = height, y = sd)) +
  geom_jitter(size = 1, alpha = 0.4,width = 0.2) +
  geom_smooth(se = F, colour = "black", size = 0.5, linetype = "dashed") +
  labs(x = "Height", y = "Standard Deviation") +
  theme_bw(base_size = 12)

sk_height <- ggplot(multispec_data_uncalibrated, aes(x = height, y = skew)) +
  geom_jitter(size = 1, alpha = 0.4,width = 0.2) +
  geom_smooth(se = F, colour = "black", size = 0.5, linetype = "dashed") +
  labs(x = "Height", y = "Skewness") +
  theme_bw(base_size = 12)


(cv_height + sd_height) / (sk_height + plot_spacer())

ggsave((cv_height + sd_height) / (sk_height + plot_spacer()),
       filename = "output/plots/uncalibrated/heigh_distribution.jpeg",
       width = 16, height = 10, units = "cm", dpi = 600)

#__________________________________________
### 3c. Relationships with Band

cv_band <- ggplot(multispec_data_uncalibrated,
       aes(x = band_nm, y = cv, colour = band)) +
  geom_jitter(width = 15, alpha = 0.4) +
  scale_colour_manual(values = c("blue", "green", "grey", "red", "darkred"),
                      guide = "none") +
  labs(x = "Wavelength (nm)", y = "Coefficient of Variation") +
  theme_bw(base_size = 12)

sd_band <- ggplot(multispec_data_uncalibrated,
                  aes(x = band_nm, y = sd, colour = band)) +
  geom_jitter(width = 15, alpha = 0.4) +
  scale_colour_manual(values = c("blue", "green", "grey", "red", "darkred"),
                      guide = "none") +
  labs(x = "Wavelength (nm)", y = "Standard Deviation") +
  theme_bw(base_size = 12)

sk_band <- ggplot(multispec_data_uncalibrated,
                  aes(x = band_nm, y = skew, colour = band)) +
  geom_jitter(width = 15, alpha = 0.4) +
  scale_colour_manual(values = c("blue", "green", "grey", "red", "darkred"),
                      guide = "none") +
  labs(x = "Wavelength (nm)", y = "Skewness") +
  theme_bw(base_size = 12)

kt_band <- ggplot(multispec_data_uncalibrated,
                  aes(x = band_nm, y = kurtosis, colour = band)) +
  geom_jitter(width = 15, alpha = 0.4) +
  scale_colour_manual(values = c("blue", "green", "grey", "red", "darkred"),
                      guide = "none") +
  coord_cartesian(ylim = c(0,50)) +
  labs(x = "Wavelength (nm)", y = "Kurtosis") +
  theme_bw(base_size = 12)

(cv_band + sd_band) / (sk_band + kt_band)

ggsave((cv_band + sd_band) / (sk_band + kt_band),
       filename = "output/plots/uncalibrated/distribution_summary_band.jpeg",
       width = 16, height = 10, units = "cm", dpi = 600)

