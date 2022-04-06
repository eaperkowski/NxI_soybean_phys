#####################################################################
# Libraries
#####################################################################
library(LeafArea)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(reshape)
library(plantecophys)

#####################################################################
# Load temp standardization function for Vcmax/Jmax
#####################################################################
source("/Users/eaperkowski/git/r_functions/temp_standardize.R")

#####################################################################
# Load .csv file with dry weight data, leaf CN data, and
# fluorescence data
#####################################################################
dry.wgt <- read.csv("../data/2021NxI_biomass_TLA.csv")
dry.wgt$id <- toupper(dry.wgt$id)

leaf.cn <- read.csv("../data/2021NxI_leafcn.csv")
leaf.cn$id <- toupper(leaf.cn$id)

fluorescence <- read.csv("../data/2021NxI_fluorescence.csv")
fluorescence$id <- toupper(fluorescence$id)

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
  full_join(dry.wgt) %>%
  group_by(id) %>%
  mutate(total.biomass = ifelse(is.na(root.biomass), 
                                NA, sum(focal.biomass, leaf.biomass, 
                                        stem.biomass, root.biomass, 
                                        nodule.biomass, na.rm = TRUE))) %>%
  full_join(leaf.cn) %>%
  full_join(fluorescence) %>%
  group_by(id) %>%
  summarize_all(mean) %>%
  separate(col = "id", 
           sep = "(_*)[_]_*",
           into = c("rep", "n.trt", "inoc"),
           remove = FALSE) %>%
  mutate(rep = gsub("R", "", rep),
                rep = str_pad(rep, width = 2, side = "left", pad = "0"),
         id = tolower(id),
         total.leaf.area = focal.area + total.leaf.area) %>%
  arrange(rep) %>%
  mutate(block = rep(1:4, each = 16),
         sla = (focal.area / focal.biomass), # SLA is in cm^2 g^-1
         narea = (leaf.n/100) / sla * 10000,
         leaf.cn = leaf.c / leaf.n)

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
  mutate(rep = gsub("r", "", rep),
         rep = str_pad(rep, width = 2, side = "left", pad = "0"),
         n.trt = toupper(n.trt),
         inoc = toupper(inoc),
         keep.row = ifelse(lag(CO2_r > 1501, n = 1L),"no","yes"),
         keep.row = tidyr::replace_na(keep.row, "yes")) %>%
  data.frame()
aci.merged

## Remove rows based on A/Ci fits
aci.merged$keep.row[c(31, 38, 90, 97, 148, 155, 163, 170, 214,
                      222, 229, 267, 274, 276, 282, 408, 401,
                      431, 444, 453, 446, 468, 506, 513, 519,
                      526, 698, 705, 749, 755, 802, 809, 815,
                      871, 872, 922, 929, 937, 938)] <- "no"

aci.temp <- aci.merged %>%
  filter(keep.row == "yes") %>%
  group_by(id) %>%
  summarize(leaf.temp = mean(TleafEB, na.rm = TRUE))

#####################################################################
# Extract A400, Ci:Ca, gsw values from each ID
#####################################################################
a.gs <- aci.merged %>%
  filter(keep.row == "yes") %>%
  group_by(id) %>%
  filter(CO2_r > 350 & CO2_r < 425) %>%
  filter(row_number() == 1) %>%
  dplyr::select(id, rep, n.trt, inoc, machine, A, Ci, Ca, gsw) %>%
  mutate(iwue = A / gsw,
         ci.ca = Ci / Ca) %>%
  data.frame()
a.gs

#####################################################################
# Run A/Ci curves with TPU limitation
#####################################################################
aci.tpu <- aci.merged %>%
  filter(keep.row == "yes") %>%
  fitacis(group = "id",
          varnames = list(ALEAF = "A",
                          Tleaf = "TleafEB",
                          Ci = "Ci",
                          PPFD = "Qin",
                          Rd = "rd"),
          fitTPU = TRUE,
          useRd = TRUE,
          Tcorrect = FALSE)

## Check plots for fit and check if any points need to be removed
# plot(aci.tpu[[4]])
# aci.tpu[[4]]

## Remove r4_hn_ni from list (no fit because no Rd value)
aci.tpu$r4_hn_ni <- NULL

## Do fitaci fxn for r4_hn_ni with useRd = FALSE
r4_hn_ni <- fitaci(subset(aci.merged, id == "r4_hn_ni"),
                   varnames = list(ALEAF = "A",
                                   Tleaf = "Tair",
                                   Ci = "Ci",
                                   PPFD = "Qin"),
                   fitTPU = TRUE,
                   useRd = FALSE,
                   Tcorrect = FALSE)

## Check model fit for r4_hn_ni
plot(r4_hn_ni)

## Create data frame by transposing r4_hn_ni coefficients, to be
## merged back into larger fitacis list
aci.r4_hn_ni <- data.frame(id = "r4_hn_ni",
                           t(coef(r4_hn_ni)))

## Check that transposing looks right. Should have one column with id
## on left column with column for each of Vcmax, Jmax, Rd, and TPU
head(aci.r4_hn_ni)

## Extract coefficients and separate id into rep, n.trt, and inoc.
## Also, merge r4_hn_ni coefficients. Add leaf temp for standardizing
## Vcmax and Jmax to 25 deg C
aci.coef <- aci.tpu %>%
  coef() %>%
  full_join(aci.r4_hn_ni) %>%
  dplyr::select(id,
                vcmax = Vcmax,
                jmax = Jmax,
                TPU) %>%
  separate(col = "id",
         sep = "(_*)[_]_*",
         into = c("rep", "n.trt", "inoc"),
         remove = FALSE) %>%
  mutate(rep = gsub("r", "", rep),
         rep = str_pad(rep, width = 2, side = "left", pad = "0"),
         n.trt = toupper(n.trt),
         inoc = toupper(inoc)) %>%
  arrange(rep) %>%
  mutate(block = rep(1:4, each = 16)) %>%
  full_join(aci.temp) %>%
  full_join(resp.merged) %>%
  full_join(a.gs) %>%
  full_join(leaf.traits) %>%
  data.frame()

## Create data frame with id and machine, join with aci.coef file
aci.coef <- aci.merged %>%
  group_by(id) %>%
  summarize(machine = unique(machine)) %>%
  left_join(aci.coef) %>%
  data.frame()

## Check aci.coef data frame. Should have id, n.trt, inoc, machine, 
## model coefficients, and block as designated columns
head(aci.coef)

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
# Standardize Vcmax and Jmax to 25 deg C
#####################################################################
## NOTE: tGrow is set to 30 arbitrarily, will be replaced by mean temp
## with HOBO data once experiment is taken down
aci.coef <- aci.coef %>%
  group_by(id) %>%
  mutate(vcmax25 = temp_standardize(estimate = vcmax,
                                    estimate.type = "Vcmax",
                                    standard.to = 25,
                                    tLeaf = leaf.temp,
                                    tGrow = tGrow),
         jmax25 = temp_standardize(estimate = jmax,
                                   estimate.type = "Jmax",
                                   standard.to = 25,
                                   tLeaf = leaf.temp,
                                   tGrow = tGrow),
         rd25.vcmax25 = rd25 / vcmax25, # Rd is temp standardized, so using Vcmax25
         jmax25.vcmax25 = jmax25 / vcmax25,
         vcmax.gs = vcmax / gsw,
         narea.gs = narea / gsw,
         pnue = A / narea) %>% # gs is not temp standardized, so using Vcmax
  dplyr::select(id, rep, n.trt, inoc, block, machine, anet = A, vcmax, vcmax25, 
                jmax, jmax25, rd, rd25, TPU, rd25.vcmax25, jmax25.vcmax25, gsw, 
                ci.ca, pnue, iwue, vcmax.gs, narea.gs, sla, focal.area, 
                focal.biomass, leaf.n, leaf.cn, narea, everything(), 
                -leaf.temp, -tleaf, -(ambient.humidity:ambient.temp)) %>%
  dplyr::rename_all(tolower) %>%
  data.frame()

head(aci.coef)

# #####################################################################
# # Run A/Ci curves without TPU limitation
# #####################################################################
# aci.notpu <- aci.merged %>%
#   filter(keep.row == "yes") %>%
#   fitacis(group = "id",
#           varnames = list(ALEAF = "A",
#                           Tleaf = "TleafEB",
#                           Ci = "Ci",
#                           PPFD = "Qin",
#                           Rd = "resp"),
#           fitTPU = FALSE,
#           useRd = TRUE,
#           Tcorrect = FALSE)
# 
# ## Remove r4_hn_ni from list (no fit because no Rd value)
# aci.notpu$r4_hn_ni <- NULL
# 
# ## Do fitaci fxn for r4_hn_ni with useRd = FALSE
# r4_hn_ni.notpu <- fitaci(subset(aci.merged, id == "r4_hn_ni"),
#                    varnames = list(ALEAF = "A",
#                                    Tleaf = "TleafEB",
#                                    Ci = "Ci",
#                                    PPFD = "Qin"),
#                    fitTPU = FALSE,
#                    useRd = FALSE,
#                    Tcorrect = FALSE)
# 
# ## Check model fit for r4_hn_ni
# plot(r4_hn_ni.notpu)
# 
# ## Create data frame by transposing r4_hn_ni coefficients, to be
# ## merged back into larger fitacis list
# aci.r4_hn_ni.notpu <- data.frame(id = "r4_hn_ni",
#                            t(coef(r4_hn_ni.notpu)))
# 
# ## Check that transposing looks right. Should have one column with id
# ## on left column with column for each of Vcmax, Jmax, Rd, and TPU
# head(aci.r4_hn_ni.notpu)
# 
# ## Extract coefficients and separate id into rep, n.trt, and inoc.
# ## Also, merge r4_hn_ni coefficients
# aci.coef <- aci.notpu %>%
#   coef() %>%
#   full_join(aci.r4_hn_ni.notpu) %>%
#   full_join(aci.coef) %>%
#   dplyr::select(id, rep, n.trt,
#                 inoc, block, machine,
#                 Vcmax.noTPU = Vcmax, 
#                 Jmax.noTPU = Jmax, 
#                 Rd,
#                 Vcmax.TPU,
#                 Jmax.TPU,
#                 TPU.TPU) %>%
#   mutate(Rd.Vcmax.noTPU = Rd / Vcmax.noTPU,
#          Rd.Vcmax.TPU = Rd / Vcmax.TPU) %>%
#   full_join(leaf.area) 

## Check aci.coef data frame. Should have coefficients with and without
## TPU fit. Rd is the same across all fits, so only including single Rd
## value. Includes Rd:Vcmax with and without TPU and focal leaf area. 
## Should also include concatenated id, and columns for rep, n.trt, inoc,
## block number, and machine
#head(aci.coef)

## Write .csv file for leaf trait data
write.csv(aci.coef, "../data/2021NxI_trait_data.csv", row.names = FALSE)

