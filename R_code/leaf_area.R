#####################################################################
# File summary information
#####################################################################
# Initializes leaf trait database by formatting appropriate column
# names, quantifying leaf area for both site visit trips, and 
# creates column for leaf area trait data

#####################################################################
# Libraries
#####################################################################
library(LeafArea)
library(stringr)
library(dplyr)
library(tidyr)

#####################################################################
# Quantify leaf area
#####################################################################
imagej.localaddress <- "/Applications/ImageJ.app"
imagepath <- "/Users/eaperkowski/git/joseph_greenhouse_phys_2021/leaf_area/"


leaf.area <- run.ij(path.imagej = imagej.localaddress,
                           set.directory = imagepath,
                           distance.pixel = 123,
                           known.distance = 1,
                           low.size = 0.00001,
                           set.memory = 30)

leaf.area <- leaf.area %>%
  tidyr::separate(col = "sample", 
           sep = "(_*)[_]_*",
           into = c("rep", "n.trt", "inoc")) %>%


