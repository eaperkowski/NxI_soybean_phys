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
# Load temp standardization function for Vcmax/Jmax and leaf demand
# function from Dong et al. (2022)
#####################################################################
source("/Users/eaperkowski/git/r_functions/temp_standardize.R")
source("/Users/eaperkowski/git/r_functions/leafDemand.R")

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
         rep = str_pad(rep, width = 2, side = "left", pad = "0"),
         id = tolower(id),
         total.leaf.area = focal.area + total.leaf.area) %>%
  arrange(rep) %>%
  dplyr::mutate(sla = (focal.area / focal.biomass), # SLA is in cm^2 g^-1
                narea = (focal.n/100) / sla * 10000,
                leaf.cn = focal.c / focal.n)

## Check data frame
head(leaf.traits)

#####################################################################
# Determine respiration values - to later be merged to A/Ci files
# to fit curves with explicit respiration
#####################################################################
## Load files into large list of data frames
file.list <- list.files(path = "../licor_data_cleaned/resp/",
                        recursive = TRUE,
                        pattern = "\\.csv$",
                        full.names = TRUE)
file.list <- setNames(file.list, file.list)
df.resp <- lapply(file.list, read.csv)

## Merge data frames into single data frame, summarize respiration for each
## ID (12 measurements per ID), change negative values to absolute numbers
resp.merged <- df.resp %>%
  merge_all() %>%
  group_by(id) %>%
  dplyr::select(id, A, TleafEB) %>%
  mutate(resp = abs(A)) %>%
  summarize(rd = mean(resp, na.rm = TRUE),
            tleaf = mean(TleafEB, na.rm = TRUE),
            rd25 = temp_standardize(rd,
                                    estimate.type = "Rd",
                                    pft = "C3H",
                                    standard.to = 25,
                                    tLeaf = tleaf,
                                    tGrow = 30)) %>%
  data.frame()
resp.merged

#####################################################################
# Load A/Ci curves, put in central data frame, add respiration means,
# then run fitacis function
#####################################################################
file.list <- list.files(path = "../licor_data_cleaned/aci",
                        recursive = TRUE,
                        pattern = "\\.csv$",
                        full.names = TRUE)
file.list <- setNames(file.list, file.list)
df.aci <- lapply(file.list, read.csv)

## Subset A/Ci curve to only columns necessary for 'fitaci'
## function. Add respiration mean values to each ID. Also,
## add keep row column and insert "no" when 400 ppm CO2 measurement
## is recorded after a 2000 ppm CO2 measurement. This is to later remove
## these values in the curve estimation
aci.merged <- df.aci %>%
  merge_all() %>%
  group_by(id) %>%
  dplyr::select(id, machine, A, Ci, Ca, gsw, 
         CO2_s,	CO2_r,	H2O_s,	H2O_r,
         Qin, VPDleaf, Flow,	Tair,	TleafEB) %>%
  arrange(id) %>%
  left_join(resp.merged, by = "id") %>%
  separate(col = "id",
           sep = "(_*)[_]_*",
           into = c("rep", "n.trt", "inoc"),
           remove = FALSE) %>%
  group_by(id) %>%
  mutate(rd.curvefit = temp_standardize(rd,
                                        estimate.type = "Rd",
                                        pft = "C3H",
                                        standard.to = mean(TleafEB),
                                        tLeaf = tleaf,
                                        tGrow = 30),
         rep = gsub("r", "", rep),
         rep = str_pad(rep, width = 2, side = "left", pad = "0"),
         keep.row = ifelse(lag(CO2_r > 1501, n = 1L),"no","yes"),
         keep.row = tidyr::replace_na(keep.row, "yes")) %>%
  data.frame()
aci.merged

## Remove rows based on A/Ci fits, and also include all points measured
## at 0 ppm CO2. Workshop w/ Licor noted that 0ppm CO2 turns off mixing fan. 
#aci.merged$keep.row[aci.merged$A < -1.5] <- "no"
aci.merged$keep.row[c(31, 282, 340, 431, 444, 452, 487, 488, 489, 
                      506, 620, 621, 692, 693, 694, 705, 938)] <- "no"

#####################################################################
# Rep 1 cluster
#####################################################################
aci.fits <- aci.merged %>% filter(keep.row == "yes") %>%
  fitacis(group = "id",
          varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE, useRd = TRUE)

photo.params <- coef(aci.fits) %>%
  select(id:Jmax)


#####################################################################
## Extract coefficients and separate id into rep, n.trt, and inoc. 
## Also formate rep number and add block number
#####################################################################
aci.coef <- photo.params %>%
  separate(col = "id",
           sep = "(_*)[_]_*",
           into = c("rep", "n.trt", "inoc"),
           remove = FALSE) %>%
  mutate(rep = gsub("r", "", rep),
         rep = str_pad(rep, width = 2, side = "left", pad = "0")) %>%
  arrange(rep) %>%
  mutate(block = rep(1:4, each = 16)) 


#####################################################################
# Import HOBO data, determine mean temp over experiment. Will be added
# as additional column to aci.coef based on block number and used
# for tGrow call in temp_standardize fxn
#####################################################################
file.list <- list.files(path = "../hobo_temp",
                        recursive = TRUE,
                        pattern = "\\.csv$",
                        full.names = TRUE)
file.list <- setNames(file.list, file.list)
df.tgrow <- lapply(file.list, read.csv, strip.white = TRUE)

aci.coef <- df.tgrow %>%
  merge_all() %>%
  group_by(block) %>%
  dplyr::summarize(tGrow = mean(temp, na.rm = TRUE)) %>%
  slice(-5) %>%
  right_join(aci.coef, by = "block")

#####################################################################
# Extract A400, Ci:Ca, gsw values from each ID
#####################################################################
a.gs <- aci.merged %>%
  group_by(id) %>%
  filter(CO2_r > 350 & CO2_r < 425) %>%
  filter(row_number() == 1) %>%
  dplyr::select(id, rep, n.trt, inoc, machine, A, Ci, Ca, gsw) %>%
  mutate(iwue = A / gsw,
         ci.ca = Ci / Ca,
         n.trt = tolower(n.trt),
         inoc = tolower(inoc)) %>%
  data.frame()
a.gs

#####################################################################
# Merge respiration data, A/Ci coefficiencts, net photosynthesis
# data files to leaf.traits data frame
#####################################################################
aci.coef <- aci.coef %>%
  full_join(resp.merged) %>%
  full_join(a.gs) %>%
  full_join(leaf.traits) %>%
  data.frame()


#####################################################################
# Add machine to aci.coef file
#####################################################################
aci.coef <- aci.merged %>%
  group_by(id) %>%
  summarize(leaf.temp = mean(TleafEB, na.rm = TRUE),
            machine = unique(machine)) %>%
  left_join(aci.coef) %>%
  data.frame()

## Check aci.coef data frame. Should have id, n.trt, inoc, machine, 
## model coefficients, and block as designated columns
head(aci.coef)

#####################################################################
# Standardize Vcmax and Jmax to 25 deg C
#####################################################################
aci.coef <- aci.coef %>%
  group_by(id) %>%
  dplyr::rename(Tleaf = leaf.temp) %>%
  mutate(vcmax25 = temp_standardize(estimate = Vcmax,
                                    estimate.type = "Vcmax",
                                    standard.to = 25,
                                    tLeaf = Tleaf,
                                    tGrow = tGrow),
         jmax25 = temp_standardize(estimate = Jmax,
                                   estimate.type = "Jmax",
                                   standard.to = 25,
                                   tLeaf = Tleaf,
                                   tGrow = tGrow),
         rd25.vcmax25 = rd25 / vcmax25,
         jmax25.vcmax25 = jmax25 / vcmax25,
         vcmax.gs = Vcmax / gsw,
         narea.gs = narea / gsw,
         pnue = A / narea) %>%
  dplyr::select(id, rep, n.trt, inoc, block, machine, anet = A, vcmax = Vcmax, 
                vcmax25, jmax = Jmax, jmax25, rd, rd25, rd25.vcmax25, 
                jmax25.vcmax25, gsw, ci.ca, pnue, iwue, vcmax.gs, narea.gs, 
                sla, focal.area, focal.biomass, leaf.n, leaf.cn, narea, Tleaf,
                everything(), -tleaf) %>%
  dplyr::rename_all(tolower) %>%
  data.frame()

## Write .csv file for leaf trait data
write.csv(aci.coef, "../data/2021NxI_trait_data.csv", row.names = FALSE)

