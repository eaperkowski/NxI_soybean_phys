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

fluorescence <- read.csv("../data/2021NxI_fluorescence.csv")

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
## r1_hn_ni
r1_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r1_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r1_hn_ni)


photo.params <- data.frame(id = "r1_hn_ni", 
                           vcmax = r1_hn_ni$pars[[1]],
                           vcmax_se = r1_hn_ni$pars[[4]],
                           jmax = r1_hn_ni$pars[[2]],
                           jmax_se = r1_hn_ni$pars[[5]])

## r1_ln_ni
r1_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r1_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r1_ln_ni)


photo.params <- photo.params %>%
  add_row(id = "r1_ln_ni", 
          vcmax = r1_ln_ni$pars[[1]],
          vcmax_se = r1_ln_ni$pars[[4]],
          jmax = r1_ln_ni$pars[[2]],
          jmax_se = r1_ln_ni$pars[[5]])


## r1_ln_yi
r1_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r1_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r1_ln_yi)


photo.params <- photo.params %>%
  add_row(id = "r1_ln_yi", 
          vcmax = r1_ln_yi$pars[[1]],
          vcmax_se = r1_ln_yi$pars[[4]],
          jmax = r1_ln_yi$pars[[2]],
          jmax_se = r1_ln_yi$pars[[5]])

## r1_hn_yi
r1_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r1_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r1_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r1_hn_yi", 
          vcmax = r1_hn_yi$pars[[1]],
          vcmax_se = r1_hn_yi$pars[[4]],
          jmax = r1_hn_yi$pars[[2]],
          jmax_se = r1_hn_yi$pars[[5]])

#####################################################################
# Rep 2 cluster
#####################################################################
## r2_hn_ni
r2_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r2_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r2_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r2_hn_ni", 
          vcmax = r2_hn_ni$pars[[1]],
          vcmax_se = r2_hn_ni$pars[[4]],
          jmax = r2_hn_ni$pars[[2]],
          jmax_se = r2_hn_ni$pars[[5]])

## r2_ln_ni
r2_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r2_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r2_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r2_ln_ni", 
          vcmax = r2_ln_ni$pars[[1]],
          vcmax_se = r2_ln_ni$pars[[4]],
          jmax = r2_ln_ni$pars[[2]],
          jmax_se = r2_ln_ni$pars[[5]])


## r2_ln_yi
r2_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r2_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r2_ln_yi)


photo.params <- photo.params %>%
  add_row(id = "r2_ln_yi", 
          vcmax = r2_ln_yi$pars[[1]],
          vcmax_se = r2_ln_yi$pars[[4]],
          jmax = r2_ln_yi$pars[[2]],
          jmax_se = r2_ln_yi$pars[[5]])

## r2_hn_yi
r2_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r2_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE, citransition = 200)
plot(r2_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r2_hn_yi", 
          vcmax = r2_hn_yi$pars[[1]],
          vcmax_se = r2_hn_yi$pars[[4]],
          jmax = r2_hn_yi$pars[[2]],
          jmax_se = r2_hn_yi$pars[[5]])

#####################################################################
# Rep 3 cluster
#####################################################################
## r3_hn_ni
r3_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r3_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE, citransition = 200)
plot(r3_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r3_hn_ni", 
          vcmax = r3_hn_ni$pars[[1]],
          vcmax_se = r3_hn_ni$pars[[4]],
          jmax = r3_hn_ni$pars[[2]],
          jmax_se = r3_hn_ni$pars[[5]])

## r3_ln_ni
r3_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r3_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r3_ln_ni)



photo.params <- photo.params %>%
  add_row(id = "r3_ln_ni", 
          vcmax = r3_ln_ni$pars[[1]],
          vcmax_se = r3_ln_ni$pars[[4]],
          jmax = r3_ln_ni$pars[[2]],
          jmax_se = r3_ln_ni$pars[[5]])


## r3_ln_yi
r3_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r3_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE, citransition = 300)
plot(r3_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r3_ln_yi", 
          vcmax = r3_ln_yi$pars[[1]],
          vcmax_se = r3_ln_yi$pars[[4]],
          jmax = r3_ln_yi$pars[[2]],
          jmax_se = r3_ln_yi$pars[[5]])

## r2_hn_yi
r3_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r3_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE, citransition = 200)
plot(r3_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r3_hn_yi", 
          vcmax = r3_hn_yi$pars[[1]],
          vcmax_se = r3_hn_yi$pars[[4]],
          jmax = r3_hn_yi$pars[[2]],
          jmax_se = r3_hn_yi$pars[[5]])

#####################################################################
# Rep 4 cluster
#####################################################################
## r4_hn_ni
r4_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r4_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE, citransition = 300)
plot(r4_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r4_hn_ni", 
          vcmax = r4_hn_ni$pars[[1]],
          vcmax_se = r4_hn_ni$pars[[4]],
          jmax = r4_hn_ni$pars[[2]],
          jmax_se = r4_hn_ni$pars[[5]])

## r4_ln_ni
r4_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r4_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r4_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r4_ln_ni", 
          vcmax = r4_ln_ni$pars[[1]],
          vcmax_se = r4_ln_ni$pars[[4]],
          jmax = r4_ln_ni$pars[[2]],
          jmax_se = r4_ln_ni$pars[[5]])

## r4_ln_yi
r4_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r4_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r4_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r4_ln_yi", 
          vcmax = r4_ln_yi$pars[[1]],
          vcmax_se = r4_ln_yi$pars[[4]],
          jmax = r4_ln_yi$pars[[2]],
          jmax_se = r4_ln_yi$pars[[5]])

## r4_hn_yi
r4_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r4_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r4_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r4_hn_yi", 
          vcmax = r4_hn_yi$pars[[1]],
          vcmax_se = r4_hn_yi$pars[[4]],
          jmax = r4_hn_yi$pars[[2]],
          jmax_se = r4_hn_yi$pars[[5]])

#####################################################################
# Rep 5 cluster
#####################################################################
## r5_hn_ni
r5_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r5_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r5_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r5_hn_ni", 
          vcmax = r5_hn_ni$pars[[1]],
          vcmax_se = r5_hn_ni$pars[[4]],
          jmax = r5_hn_ni$pars[[2]],
          jmax_se = r5_hn_ni$pars[[5]])

## r5_ln_ni
r5_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r5_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r5_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r5_ln_ni", 
          vcmax = r5_ln_ni$pars[[1]],
          vcmax_se = r5_ln_ni$pars[[4]],
          jmax = r5_ln_ni$pars[[2]],
          jmax_se = r5_ln_ni$pars[[5]])

## r5_ln_yi
r5_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r5_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r5_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r5_ln_yi", 
          vcmax = r5_ln_yi$pars[[1]],
          vcmax_se = r5_ln_yi$pars[[4]],
          jmax = r5_ln_yi$pars[[2]],
          jmax_se = r5_ln_yi$pars[[5]])

## r5_hn_yi
r5_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r5_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r5_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r5_hn_yi", 
          vcmax = r5_hn_yi$pars[[1]],
          vcmax_se = r5_hn_yi$pars[[4]],
          jmax = r5_hn_yi$pars[[2]],
          jmax_se = r5_hn_yi$pars[[5]])

#####################################################################
# Rep 6 cluster
#####################################################################
## r6_hn_ni
r6_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r6_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r6_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r6_hn_ni", 
          vcmax = r6_hn_ni$pars[[1]],
          vcmax_se = r6_hn_ni$pars[[4]],
          jmax = r6_hn_ni$pars[[2]],
          jmax_se = r6_hn_ni$pars[[5]])

## r6_ln_ni
r6_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r6_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r6_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r6_ln_ni", 
          vcmax = r6_ln_ni$pars[[1]],
          vcmax_se = r6_ln_ni$pars[[4]],
          jmax = r6_ln_ni$pars[[2]],
          jmax_se = r6_ln_ni$pars[[5]])

## r6_ln_yi
r6_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r6_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r6_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r6_ln_yi", 
          vcmax = r6_ln_yi$pars[[1]],
          vcmax_se = r6_ln_yi$pars[[4]],
          jmax = r6_ln_yi$pars[[2]],
          jmax_se = r6_ln_yi$pars[[5]])

## r6_hn_yi
r6_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r6_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r6_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r6_hn_yi", 
          vcmax = r6_hn_yi$pars[[1]],
          vcmax_se = r6_hn_yi$pars[[4]],
          jmax = r6_hn_yi$pars[[2]],
          jmax_se = r6_hn_yi$pars[[5]])

#####################################################################
# Rep 7 cluster
#####################################################################
## r7_hn_ni
r7_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r7_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r7_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r7_hn_ni", 
          vcmax = r7_hn_ni$pars[[1]],
          vcmax_se = r7_hn_ni$pars[[4]],
          jmax = r7_hn_ni$pars[[2]],
          jmax_se = r7_hn_ni$pars[[5]])

## r7_ln_ni
r7_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r7_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r7_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r7_ln_ni", 
          vcmax = r7_ln_ni$pars[[1]],
          vcmax_se = r7_ln_ni$pars[[4]],
          jmax = r7_ln_ni$pars[[2]],
          jmax_se = r7_ln_ni$pars[[5]])

## r7_ln_yi
r7_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r7_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r7_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r7_ln_yi", 
          vcmax = r7_ln_yi$pars[[1]],
          vcmax_se = r7_ln_yi$pars[[4]],
          jmax = r7_ln_yi$pars[[2]],
          jmax_se = r7_ln_yi$pars[[5]])

## r7_hn_yi
r7_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r7_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r7_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r7_hn_yi", 
          vcmax = r7_hn_yi$pars[[1]],
          vcmax_se = r7_hn_yi$pars[[4]],
          jmax = r7_hn_yi$pars[[2]],
          jmax_se = r7_hn_yi$pars[[5]])

#####################################################################
# Rep 8 cluster
#####################################################################
## r8_hn_ni
r8_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r8_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r8_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r8_hn_ni", 
          vcmax = r8_hn_ni$pars[[1]],
          vcmax_se = r8_hn_ni$pars[[4]],
          jmax = r8_hn_ni$pars[[2]],
          jmax_se = r8_hn_ni$pars[[5]])

## r6_ln_ni
r8_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r8_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r8_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r8_ln_ni", 
          vcmax = r8_ln_ni$pars[[1]],
          vcmax_se = r8_ln_ni$pars[[4]],
          jmax = r8_ln_ni$pars[[2]],
          jmax_se = r8_ln_ni$pars[[5]])

## r8_ln_yi
r8_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r8_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r8_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r8_ln_yi", 
          vcmax = r8_ln_yi$pars[[1]],
          vcmax_se = r8_ln_yi$pars[[4]],
          jmax = r8_ln_yi$pars[[2]],
          jmax_se = r8_ln_yi$pars[[5]])

## r8_hn_yi
r8_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r8_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r8_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r8_hn_yi", 
          vcmax = r8_hn_yi$pars[[1]],
          vcmax_se = r8_hn_yi$pars[[4]],
          jmax = r8_hn_yi$pars[[2]],
          jmax_se = r8_hn_yi$pars[[5]])

#####################################################################
# Rep 9 cluster
#####################################################################
## r9_hn_ni
r9_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r9_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r9_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r9_hn_ni", 
          vcmax = r9_hn_ni$pars[[1]],
          vcmax_se = r9_hn_ni$pars[[4]],
          jmax = r9_hn_ni$pars[[2]],
          jmax_se = r9_hn_ni$pars[[5]])

## r9_ln_ni
r9_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r9_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r9_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r9_ln_ni", 
          vcmax = r9_ln_ni$pars[[1]],
          vcmax_se = r9_ln_ni$pars[[4]],
          jmax = r9_ln_ni$pars[[2]],
          jmax_se = r9_ln_ni$pars[[5]])

## r9_ln_yi
r9_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r9_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r9_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r9_ln_yi", 
          vcmax = r9_ln_yi$pars[[1]],
          vcmax_se = r9_ln_yi$pars[[4]],
          jmax = r9_ln_yi$pars[[2]],
          jmax_se = r9_ln_yi$pars[[5]])

## r9_hn_yi
r9_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r9_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r9_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r9_hn_yi", 
          vcmax = r9_hn_yi$pars[[1]],
          vcmax_se = r9_hn_yi$pars[[4]],
          jmax = r9_hn_yi$pars[[2]],
          jmax_se = r9_hn_yi$pars[[5]])

#####################################################################
# Rep 10 cluster
#####################################################################
## r10_hn_ni
r10_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r10_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r10_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r10_hn_ni", 
          vcmax = r10_hn_ni$pars[[1]],
          vcmax_se = r10_hn_ni$pars[[4]],
          jmax = r10_hn_ni$pars[[2]],
          jmax_se = r10_hn_ni$pars[[5]])

## r10_ln_ni
r10_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r10_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r10_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r10_ln_ni", 
          vcmax = r10_ln_ni$pars[[1]],
          vcmax_se = r10_ln_ni$pars[[4]],
          jmax = r10_ln_ni$pars[[2]],
          jmax_se = r10_ln_ni$pars[[5]])

## r10_ln_yi
r10_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r10_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r10_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r10_ln_yi", 
          vcmax = r10_ln_yi$pars[[1]],
          vcmax_se = r10_ln_yi$pars[[4]],
          jmax = r10_ln_yi$pars[[2]],
          jmax_se = r10_ln_yi$pars[[5]])

## r10_hn_yi
r10_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r10_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r10_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r10_hn_yi", 
          vcmax = r10_hn_yi$pars[[1]],
          vcmax_se = r10_hn_yi$pars[[4]],
          jmax = r10_hn_yi$pars[[2]],
          jmax_se = r10_hn_yi$pars[[5]])

#####################################################################
# Rep 11 cluster
#####################################################################
## r11_hn_ni
r11_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r11_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r11_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r11_hn_ni", 
          vcmax = r11_hn_ni$pars[[1]],
          vcmax_se = r11_hn_ni$pars[[4]],
          jmax = r11_hn_ni$pars[[2]],
          jmax_se = r11_hn_ni$pars[[5]])

## r11_ln_ni
r11_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r11_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r11_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r11_ln_ni", 
          vcmax = r11_ln_ni$pars[[1]],
          vcmax_se = r11_ln_ni$pars[[4]],
          jmax = r11_ln_ni$pars[[2]],
          jmax_se = r11_ln_ni$pars[[5]])

## r11_ln_yi
r11_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r11_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r11_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r11_ln_yi", 
          vcmax = r11_ln_yi$pars[[1]],
          vcmax_se = r11_ln_yi$pars[[4]],
          jmax = r11_ln_yi$pars[[2]],
          jmax_se = r11_ln_yi$pars[[5]])

## r11_hn_yi
r11_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r11_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r11_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r11_hn_yi", 
          vcmax = r11_hn_yi$pars[[1]],
          vcmax_se = r11_hn_yi$pars[[4]],
          jmax = r11_hn_yi$pars[[2]],
          jmax_se = r11_hn_yi$pars[[5]])

#####################################################################
# Rep 12 cluster
#####################################################################
## r12_hn_ni
r12_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r12_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r12_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r12_hn_ni", 
          vcmax = r12_hn_ni$pars[[1]],
          vcmax_se = r12_hn_ni$pars[[4]],
          jmax = r12_hn_ni$pars[[2]],
          jmax_se = r12_hn_ni$pars[[5]])

## r12_ln_ni
r12_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r12_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r12_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r12_ln_ni", 
          vcmax = r12_ln_ni$pars[[1]],
          vcmax_se = r12_ln_ni$pars[[4]],
          jmax = r12_ln_ni$pars[[2]],
          jmax_se = r12_ln_ni$pars[[5]])

## r12_ln_yi
r12_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r12_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r12_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r12_ln_yi", 
          vcmax = r12_ln_yi$pars[[1]],
          vcmax_se = r12_ln_yi$pars[[4]],
          jmax = r12_ln_yi$pars[[2]],
          jmax_se = r12_ln_yi$pars[[5]])

## r12_hn_yi
r12_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r12_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r12_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r12_hn_yi", 
          vcmax = r12_hn_yi$pars[[1]],
          vcmax_se = r12_hn_yi$pars[[4]],
          jmax = r12_hn_yi$pars[[2]],
          jmax_se = r12_hn_yi$pars[[5]])

#####################################################################
# Rep 13 cluster
#####################################################################
## r13_hn_ni
r13_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r13_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r13_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r13_hn_ni", 
          vcmax = r13_hn_ni$pars[[1]],
          vcmax_se = r13_hn_ni$pars[[4]],
          jmax = r13_hn_ni$pars[[2]],
          jmax_se = r13_hn_ni$pars[[5]])

## r13_ln_ni
r13_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r13_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r13_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r13_ln_ni", 
          vcmax = r13_ln_ni$pars[[1]],
          vcmax_se = r13_ln_ni$pars[[4]],
          jmax = r13_ln_ni$pars[[2]],
          jmax_se = r13_ln_ni$pars[[5]])

## r13_ln_yi
r13_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r13_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r13_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r13_ln_yi", 
          vcmax = r13_ln_yi$pars[[1]],
          vcmax_se = r13_ln_yi$pars[[4]],
          jmax = r13_ln_yi$pars[[2]],
          jmax_se = r13_ln_yi$pars[[5]])

## r13_hn_yi
r13_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r13_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r13_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r13_hn_yi", 
          vcmax = r13_hn_yi$pars[[1]],
          vcmax_se = r13_hn_yi$pars[[4]],
          jmax = r13_hn_yi$pars[[2]],
          jmax_se = r13_hn_yi$pars[[5]])

#####################################################################
# Rep 14 cluster
#####################################################################
## r14_hn_ni
r14_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r14_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r14_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r14_hn_ni", 
          vcmax = r14_hn_ni$pars[[1]],
          vcmax_se = r14_hn_ni$pars[[4]],
          jmax = r14_hn_ni$pars[[2]],
          jmax_se = r14_hn_ni$pars[[5]])

## r14_ln_ni
r14_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r14_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r14_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r14_ln_ni", 
          vcmax = r14_ln_ni$pars[[1]],
          vcmax_se = r14_ln_ni$pars[[4]],
          jmax = r14_ln_ni$pars[[2]],
          jmax_se = r14_ln_ni$pars[[5]])

## r14_ln_yi
r14_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r14_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r14_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r14_ln_yi", 
          vcmax = r14_ln_yi$pars[[1]],
          vcmax_se = r14_ln_yi$pars[[4]],
          jmax = r14_ln_yi$pars[[2]],
          jmax_se = r14_ln_yi$pars[[5]])

## r14_hn_yi
r14_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r14_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r14_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r14_hn_yi", 
          vcmax = r14_hn_yi$pars[[1]],
          vcmax_se = r14_hn_yi$pars[[4]],
          jmax = r14_hn_yi$pars[[2]],
          jmax_se = r14_hn_yi$pars[[5]])

#####################################################################
# Rep 15 cluster
#####################################################################
## r15_hn_ni
r15_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r15_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r15_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r15_hn_ni", 
          vcmax = r15_hn_ni$pars[[1]],
          vcmax_se = r15_hn_ni$pars[[4]],
          jmax = r15_hn_ni$pars[[2]],
          jmax_se = r15_hn_ni$pars[[5]])

## r15_ln_ni
r15_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r15_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r15_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r15_ln_ni", 
          vcmax = r15_ln_ni$pars[[1]],
          vcmax_se = r15_ln_ni$pars[[4]],
          jmax = r15_ln_ni$pars[[2]],
          jmax_se = r15_ln_ni$pars[[5]])

## r15_ln_yi
r15_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r15_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r15_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r15_ln_yi", 
          vcmax = r15_ln_yi$pars[[1]],
          vcmax_se = r15_ln_yi$pars[[4]],
          jmax = r15_ln_yi$pars[[2]],
          jmax_se = r15_ln_yi$pars[[5]])

## r15_hn_yi
r15_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r15_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r15_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r15_hn_yi", 
          vcmax = r15_hn_yi$pars[[1]],
          vcmax_se = r15_hn_yi$pars[[4]],
          jmax = r15_hn_yi$pars[[2]],
          jmax_se = r15_hn_yi$pars[[5]])

#####################################################################
# Rep 16 cluster
#####################################################################
## r16_hn_ni
r16_hn_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r16_hn_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r16_hn_ni)

photo.params <- photo.params %>%
  add_row(id = "r16_hn_ni", 
          vcmax = r16_hn_ni$pars[[1]],
          vcmax_se = r16_hn_ni$pars[[4]],
          jmax = r16_hn_ni$pars[[2]],
          jmax_se = r16_hn_ni$pars[[5]])

## r16_ln_ni
r16_ln_ni <- aci.merged %>% filter(keep.row == "yes" & id == "r16_ln_ni") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r16_ln_ni)

photo.params <- photo.params %>%
  add_row(id = "r16_ln_ni", 
          vcmax = r16_ln_ni$pars[[1]],
          vcmax_se = r16_ln_ni$pars[[4]],
          jmax = r16_ln_ni$pars[[2]],
          jmax_se = r16_ln_ni$pars[[5]])

## r16_ln_yi
r16_ln_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r16_ln_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r16_ln_yi)

photo.params <- photo.params %>%
  add_row(id = "r16_ln_yi", 
          vcmax = r16_ln_yi$pars[[1]],
          vcmax_se = r16_ln_yi$pars[[4]],
          jmax = r16_ln_yi$pars[[2]],
          jmax_se = r16_ln_yi$pars[[5]])

## r16_hn_yi
r16_hn_yi <- aci.merged %>% filter(keep.row == "yes" & id == "r16_hn_yi") %>%
  fitaci(varnames = list(ALEAF = "A",
                         Tleaf = "TleafEB",
                         Ci = "Ci",
                         PPFD = "Qin", Rd = "rd.curvefit"),
         fitTPU = FALSE, Tcorrect = FALSE)
plot(r16_hn_yi)

photo.params <- photo.params %>%
  add_row(id = "r16_hn_yi", 
          vcmax = r16_hn_yi$pars[[1]],
          vcmax_se = r16_hn_yi$pars[[4]],
          jmax = r16_hn_yi$pars[[2]],
          jmax_se = r16_hn_yi$pars[[5]])

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
  mutate(vcmax25 = temp_standardize(estimate = vcmax,
                                    estimate.type = "Vcmax",
                                    standard.to = 25,
                                    tLeaf = Tleaf,
                                    tGrow = tGrow),
         jmax25 = temp_standardize(estimate = jmax,
                                   estimate.type = "Jmax",
                                   standard.to = 25,
                                   tLeaf = Tleaf,
                                   tGrow = tGrow),
         rd25.vcmax25 = rd25 / vcmax25,
         jmax25.vcmax25 = jmax25 / vcmax25,
         vcmax.gs = vcmax / gsw,
         narea.gs = narea / gsw,
         pnue = A / narea) %>%
  dplyr::select(id, rep, n.trt, inoc, block, machine, anet = A, vcmax, vcmax25, 
                jmax, jmax25, rd, rd25, rd25.vcmax25, jmax25.vcmax25, gsw, 
                ci.ca, pnue, iwue, vcmax.gs, narea.gs, sla, focal.area, 
                focal.biomass, leaf.n, leaf.cn, narea, Tleaf,
                everything(), -tleaf, -vcmax_se, -jmax_se) %>%
  dplyr::rename_all(tolower) %>%
  data.frame()

## Write .csv file for leaf trait data
write.csv(aci.coef, "../data/2021NxI_trait_data.csv", row.names = FALSE)

