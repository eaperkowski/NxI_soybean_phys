#####################################################################
# Libraries
#####################################################################
library(LeafArea)
library(stringr)
library(dplyr)
library(tidyr)
library(car)
library(emmeans)

#####################################################################
# Quantify leaf area
#####################################################################
imagej.localaddress <- "/Applications/ImageJ.app"
imagepath <- "../leaf_area/"

leaf.area <- run.ij(path.imagej = imagej.localaddress,
                    set.directory = imagepath,
                    distance.pixel = 117.9034,
                    known.distance = 1,
                    low.size = 0.00001,
                    set.memory = 30)

leaf.area <- leaf.area %>%
  tidyr::separate(col = "sample", 
           sep = "(_*)[_]_*",
           into = c("rep", "n.trt", "inoc")) %>%
  dplyr::mutate(rep = gsub("R", "", rep))

gsub("R", "", leaf.area$sample)

leaf.area$block[leaf.area$rep == "R1" |
                  leaf.area$rep == "R2" |
                  leaf.area$rep == "R3"]


