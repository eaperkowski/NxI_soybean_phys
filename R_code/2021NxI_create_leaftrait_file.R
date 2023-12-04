#####################################################################
# Libraries
#####################################################################
library(LeafArea)
library(dplyr)
library(stringr)
library(readr)
library(reshape)
library(plantecophys)
library(tidyverse)

#####################################################################
# Load and clean biomass, CN, fluorescence data
#####################################################################
## Read in biomass file
biomass <- read.csv("../data/2021NxI_biomass_TLA.csv")

## List files
file.list <- list.files("../costech_results",
                        recursive = TRUE,
                        pattern = "\\.csv$",
                        full.names = TRUE)

file.list <- setNames(file.list, stringr::str_extract(basename(file.list), 
                                                      '.*(?=\\.csv)'))

## Read files, merge to data frame, separate by organ type
concat.plates <- lapply(file.list, read.csv) %>%
  reshape::merge_all() %>%
  filter(sample.type == "unknown" & id != "QC") %>%
  separate(id, into = c("rep", "n.trt", "inoc", "organ")) %>%
  unite(col = "id", rep:inoc) %>%
  dplyr::select(id, organ, n.weight.percent, c.weight.percent)


## Create different objects for each organ, rename cols
focal <- concat.plates %>%
  filter(organ == "focal") %>%
  dplyr::select(id, focal.n = n.weight.percent, focal.c = c.weight.percent)

leaves <- concat.plates %>%
  filter(organ == "tl") %>%
  dplyr::select(id, leaf.n = n.weight.percent, leaf.c = c.weight.percent)

stems <- concat.plates %>%
  filter(organ == "ts") %>%
  dplyr::select(id, stem.n = n.weight.percent, stem.c = c.weight.percent)

roots <- concat.plates %>%
  filter(organ == "tr") %>%
  dplyr::select(id, root.n = n.weight.percent, root.c = c.weight.percent)

nods <- concat.plates %>%
  filter(organ == "nod") %>%
  dplyr::select(id, nodule.n = n.weight.percent, nodule.c = c.weight.percent)

## Merge organs back into data frame
leaf.cn <- biomass %>%
  full_join(focal) %>%
  full_join(leaves) %>%
  full_join(stems) %>%
  full_join(roots) %>%
  full_join(nods) %>%
  filter(!row_number() %in% 18) %>%
  group_by(id) %>%
  mutate(total.leaf.n = (focal.biomass + leaf.biomass) * (leaf.n/100),
         total.stem.n = stem.biomass * (stem.n/100),
         total.root.n = ifelse(is.na(root.biomass), NA,
                               root.biomass * (root.n/100)),
         total.root.c = ifelse(is.na(root.biomass), NA,
                               root.biomass * (root.c/100)),
         total.nod.n = ifelse(nodule.biomass == 0, 0,
                              ifelse(is.na(root.biomass), NA,
                                     nodule.biomass * (nodule.n/100))),
         total.nod.c = ifelse(nodule.biomass == 0, 0,
                              ifelse(is.na(root.biomass), NA,
                                     nodule.biomass * (nodule.c/100)))) %>%
  mutate(wp.total.n = ifelse(is.na(root.biomass) | is.na(stem.biomass) | 
                               is.na(leaf.biomass) | is.na(leaf.n) | is.na(stem.n) | 
                               is.na(root.n), NA,
                             sum(total.leaf.n, total.root.n, total.stem.n, 
                                 total.nod.n, na.rm = TRUE)),
         bg.total.c = ifelse(is.na(root.biomass) | is.na(root.c), NA, 
                             sum(total.root.c, total.nod.c, na.rm = TRUE)),
         n.cost = bg.total.c / wp.total.n,
         total.biomass = ifelse(is.na(root.biomass), 
                                NA, sum(focal.biomass, leaf.biomass, 
                                        stem.biomass, root.biomass, 
                                        nodule.biomass, na.rm = TRUE)),
         bvr = total.biomass / 6) %>%
  separate(col = "id", 
           sep = "(_*)[_]_*",
           into = c("rep", "n.trt", "inoc"),
           remove = FALSE) %>%
  mutate(rep = gsub("r", "", rep),
         rep = str_pad(rep, width = 2, side = "left", pad = "0"),
         id = tolower(id)) %>%
  arrange(rep)

length(which(!is.na(leaf.cn$n.cost))) ## Should read 54 if done correctly

#####################################################################
# Determine leaf areas
#####################################################################
## Path for leaf area and image J software
imagej.localaddress <- "/Applications/ImageJ.app"
imagepath <- "../leaf_area/"

## Run 'LeafArea' function
leaf.area <- run.ij(path.imagej = imagej.localaddress,
                    set.directory = imagepath,
                    distance.pixel = 117.9034,
                    known.distance = 1,
                    set.memory = 30)
head(leaf.area)

## Separate id into rep, n.trt, and inoc and rename leaf area column
## Then, remove R from rep name a add a leading zero for all single 
## digit reps. Also add block and other leaf traits (sla, narea)
leaf.traits <- leaf.area %>%
  dplyr::rename(id = sample,
                focal.area = total.leaf.area) %>%
  mutate(id = tolower(id)) %>%
  full_join(leaf.cn) %>%
  group_by(id) %>%
  separate(col = "id", 
           sep = "(_*)[_]_*",
           into = c("rep", "n.trt", "inoc"),
           remove = FALSE) %>%
  mutate(rep = gsub("r", "", rep),
         block = ifelse(rep == 1 | rep == 2 | rep == 3 | rep == 4, 1, 
                        ifelse(rep == 5 | rep == 6 | rep == 7 | rep == 8,
                               2, 
                               ifelse(rep == 9 | rep == 10 | rep == 11 | rep == 12,
                                      3, ifelse(rep == 13 | rep == 14 | rep == 15 | rep == 16, 4, 
                                                NA)))),
         rep = str_pad(rep, width = 2, side = "left", pad = "0"),
         id = tolower(id),
         total.leaf.area = focal.area + total.leaf.area,
         n.trt = ifelse(n.trt == "hn", 630, 70),
         leaf.biomass = leaf.biomass + focal.biomass,
         nod.root.biomass = nodule.biomass / root.biomass) %>%
  arrange(rep) %>%
  select(id, n.trt, inoc, rep, block, leaf.biomass:nod.root.biomass, total.biomass,
         total.leaf.area, leaf.n:n.cost, bvr)

## Check data frame
head(leaf.traits)


## Write .csv file for leaf trait data
write.csv(leaf.traits, "../data/2021NxI_trait_data.csv", row.names = FALSE)

