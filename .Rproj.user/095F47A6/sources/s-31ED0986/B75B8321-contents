#########################################################
##                                                     ##
##  RainDrop: Multispectral indicators of Biodiversity ##
##                                                     ##
##            Percent cover raw data processing        ## 
##                                                     ##
##                   Uncalibrated                      ##
##                                                     ##
##                   Nov 30th 2021                     ##
##                                                     ##
#########################################################
# updated 2022-09-02
rm(list = ls())

library(tidyverse)
library(patchwork)

# lib paths for uni computer
source("../lib_path_uni_computer.R")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Loading data ####

# DRAGNET
dragnet_cover_raw <- read_csv("../../DRAGNet/Raw data/DRAGNET_percent_cover_2021.csv")
dragnet_biomass_raw <- read_csv("../../DRAGNet/Raw data/DRAGNET_dry_biomass_2021.csv")
weight_calibration <- read_csv("../../DRAGNet/Raw data/weight_calibrations_2021.csv")

tin_weight <- weight_calibration[weight_calibration$bag_type == "Tin",]$weight

# DROUGHTNET
load("../../DROUGHTNet/percent_cover_2021.RData")
droughtnet_cover_raw <- percentcover_2021
rm(percentcover_2021)

droughtnet_biomass_raw <- read_csv("../../DROUGHTNet/droughtnet_biomass_july2021.csv") %>% 
  mutate(`Functional Group` = case_when(
    `Functional Group` == "Grass" ~ "Graminoids",
    `Functional Group` == "Forb" ~ "Forbs",
    `Functional Group` == "Legume" ~ "Legumes",
    `Functional Group` == "Moss" ~ "Bryophytes"),
         Treatment = if_else(Treatment == "Procedural", "Control", Treatment))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Percent Cover ####

# dragnet
percent_cover_dragnet <- dragnet_cover_raw %>% 
  # adjustment to column values to align with multispectral data
  mutate(group = if_else(group == "Forbe" | group == "Forbs", "Forb", group),
         date = as.Date(date, format = "%d/%m/%Y"),
         season = if_else(date < "2021-09-01", "Early", "Late"),
         treatment = case_when(
           treatment == "NPK+" ~ "NPK",
           treatment == "Disturbance" ~ "Disturbance",
           treatment == "NPK disturbance" ~ "NPK + Disturbance",
           treatment == "NPK cessation" ~ "NPK Cessation")) %>% 
  # filter to only include the early season and at species level
  filter(season == "Early" &
           is.na(species) == FALSE) %>% 
  group_by(block, treatment) %>% 
  # diversity indices for each group - here cover is a proportion
  summarise(experiment = "DRAGNET",
            date = date[1], time = time[1],
            richness = n(),
            simpsons = sum(cover^2),
            shannon = -sum(cover * log(cover))) %>% 
  ungroup()

# droughtnet
percent_cover_droughtnet <- droughtnet_cover_raw %>% 
  mutate(cover = percent_cover/100) %>% 
  filter(species_level == 1) %>% 
  group_by(block, treatment) %>% 
  summarise(experiment = "DROUGHTNET", 
            date = date[1], time = NA,
            richness = n(),
            simpsons = sum(cover^2),
            shannon = -sum(cover * log(cover))) %>% 
  ungroup()

# combine
percent_cover_july2021 <- bind_rows(percent_cover_dragnet, 
                                    percent_cover_droughtnet)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Biomass ####

# dragnet
biomass_dragnet <- dragnet_biomass_raw %>% 
  dplyr::select(date, block, treatment, clip, bag_type, 
                functional_group, biomass_uc = weight) %>% 
  # Changing dates and adding collection season
  mutate(date = as.Date(date, format = "%d/%m/%Y"),
         season = if_else(date < "2021-09-01", "Early", "Late")) %>% 
  filter(season == "Early") %>% 
  # adding the weight calibrations from the bags
  left_join(x = ., y = weight_calibration, by = "bag_type") %>% 
  mutate(biomass = biomass_uc - weight,
         biomass = if_else(biomass < 0, 0, biomass)) %>% 
  # Summarising biomass over clips
  group_by(block, treatment, functional_group) %>% 
  summarise(experiment = "DRAGNET",
            date = date[1], biomass = sum(biomass)) %>% 
  ungroup() %>% 
  dplyr::select(5,1,2,4,3,6)

# droughtnet - only has early growing season
biomass_droughtnet <- droughtnet_biomass_raw %>% 
  dplyr::rename_all(.funs = tolower) %>% # SO USEFUL
  dplyr::select(block, treatment, functional_group = `functional group`,
                date, biomass = `dry biomass (g)`) %>% 
  mutate(date = as.Date(date, format = "%d/%m/%Y"),
         biomass = if_else(is.na(biomass) == TRUE, 0, biomass),
         experiment = "DROUGHTNET") %>% 
  dplyr::select(4,1,2,6,3,5)

# combine
biomass_july2021 <- bind_rows(biomass_dragnet, biomass_droughtnet)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Save data ####

save(percent_cover_july2021, biomass_july2021, file = "data/biodiversity_data_early.RData")
save(percent_cover_dragnet, biomass_dragnet, file = "data/biodiversity_data_dragnet_early.RData")


