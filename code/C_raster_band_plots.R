#########################################################
##                                                     ##
##  RainDrop: Multispectral indicators of Biodiversity ##
##                                                     ##
##            TIFF raster data extraction              ## 
##                                                     ##
##               DRAGNet uncalibrated                  ##
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

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load rasters ####

b1586 <- raster("../cropped_image_data/un_corrected_20211014/1586MS_20210629_DRAGNET_D_ND_2_B_1.tif")
g1587 <- raster("../cropped_image_data/un_corrected_20211014/1587MS_20210629_DRAGNET_D_ND_2_G_1.tif")
r1588 <- raster("../cropped_image_data/un_corrected_20211014/1588MS_20210629_DRAGNET_D_ND_2_R_1.tif")
re1589 <- raster("../cropped_image_data/un_corrected_20211014/1589MS_20210629_DRAGNET_D_ND_2_RE_1.tif")
nir1590 <- raster("../cropped_image_data/un_corrected_20211014/1590MS_20210629_DRAGNET_D_ND_2_NIR_1.tif")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Convert to data frames ####

b_df <- as.data.frame(as(b1586, "SpatialPixelsDataFrame")) %>% 
  dplyr::rename(value = X1586MS_20210629_DRAGNET_D_ND_2_B_1)

g_df <- as.data.frame(as(g1587, "SpatialPixelsDataFrame")) %>% 
  dplyr::rename(value = X1587MS_20210629_DRAGNET_D_ND_2_G_1)

r_df <- as.data.frame(as(r1588, "SpatialPixelsDataFrame")) %>% 
  dplyr::rename(value = X1588MS_20210629_DRAGNET_D_ND_2_R_1)

re_df <- as.data.frame(as(re1589, "SpatialPixelsDataFrame")) %>% 
  dplyr::rename(value = X1589MS_20210629_DRAGNET_D_ND_2_RE_1)

nir_df <- as.data.frame(as(nir1590, "SpatialPixelsDataFrame")) %>% 
  dplyr::rename(value = X1590MS_20210629_DRAGNET_D_ND_2_NIR_1)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. DN value plots ####

b_plot <- ggplot(b_df) +
  geom_tile(aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "black", high = "blue") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_colorbar("DN value", barwidth = 0.5, barheight = 5)) +
  labs(title = "Blue - 450nm") +
  theme(legend.title = element_text(size = 8),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

g_plot <- ggplot(g_df) +
  geom_tile(aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "black", high = "green") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_colorbar("DN value", barwidth = 0.5, barheight = 5)) +
  labs(title = "Green - 560nm") +
  theme(legend.title = element_text(size = 8),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

r_plot <- ggplot(r_df) +
  geom_tile(aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "black", high = "red") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_colorbar("DN value", barwidth = 0.5, barheight = 5)) +
  labs(title = "Red - 650nm") +
  theme(legend.title = element_text(size = 8),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

re_plot <- ggplot(re_df) +
  geom_tile(aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "black", high = "darkred") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_colorbar("DN value", barwidth = 0.5, barheight = 5)) +
  labs(title = "Red edge - 730nm") +
  theme(legend.title = element_text(size = 8),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

nir_plot <- ggplot(nir_df) +
  geom_tile(aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "black", high = "grey") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_colorbar("DN value", barwidth = 0.5, barheight = 5)) +
  labs(title = "Near Infrared - 840nm") +
  theme(legend.title = element_text(size = 8),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Distributions ####

b_d <- ggplot(b_df, aes(x = value)) +
  geom_density(fill = "blue", alpha = 0.4) +
  labs(x = "Digital Number (DN) value", y = "Density") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 5))

g_d <- ggplot(g_df, aes(x = value)) +
  geom_density(fill = "green", alpha = 0.4) +
  labs(x = "Digital Number (DN) value", y = "Density") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 5))

r_d <- ggplot(r_df, aes(x = value)) +
  geom_density(fill = "red", alpha = 0.4) +
  labs(x = "Digital Number (DN) value", y = "Density") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 5))

re_d <- ggplot(b_df, aes(x = value)) +
  geom_density(fill = "darkred", alpha = 0.4) +
  labs(x = "Digital Number (DN) value", y = "Density") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 5))

nir_d <- ggplot(b_df, aes(x = value)) +
  geom_density(fill = "grey", alpha = 0.4) +
  labs(x = "Digital Number (DN) value", y = "Density") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 5))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Large Plot ####

full_plot <- b_plot + g_plot + r_plot + re_plot + nir_plot +
  b_d + g_d + r_d + re_d + nir_d +
  plot_layout(nrow = 2, heights = c(3,1))

ggsave(full_plot, filename = "output/plots/dragnet_uncalibrated/band_plot.jpeg",
       width = 52, height = 12, units = "cm", dpi = 600)


