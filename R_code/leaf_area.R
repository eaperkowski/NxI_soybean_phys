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

## Separate id into rep, n.trt, and inoc and rename leaf area column
## Then, remove R from rep name a add a leading zero for all single 
## digit reps
leaf.area <- leaf.area %>%
  separate(col = "sample", 
           sep = "(_*)[_]_*",
           into = c("rep", "n.trt", "inoc"),
           remove = FALSE) %>%
  mutate(rep = gsub("R", "", rep),
                rep = str_pad(rep, width=2, side="left", pad="0"),
         sample = tolower(sample)) %>%
  dplyr::rename(focal.area = total.leaf.area,
                id = sample) %>%
  arrange(rep) %>%
  mutate(block = rep(1:4, each = 16))

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
  summarize(resp = mean(A, na.rm = TRUE)) %>%
  mutate(resp = abs(resp))
head(resp.merged)

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
## function. Add respiration mean values to each ID
aci.merged <- df.aci %>%
  merge_all() %>%
  group_by(id) %>%
  select(id, A, Ci, gsw, 
         CO2_s,	CO2_r,	H2O_s,	H2O_r,
         Qin, VPDleaf, Flow,	Tair,	Tleaf) %>%
  arrange(id) %>%
  left_join(resp.merged, by = "id") %>%
  separate(col = "id",
           sep = "(_*)[_]_*",
           into = c("rep", "n.trt", "inoc"),
           remove = FALSE) %>%
  mutate(rep = gsub("r", "", rep),
         rep = str_pad(rep, width = 2, side = "left", pad = "0"),
         n.trt = toupper(n.trt),
         inoc = toupper(inoc)) %>%
  data.frame()

head(aci.merged)

#####################################################################
# Run A/Ci curves with TPU limitation
#####################################################################
aci.tpu <- fitacis(aci.merged, 
                   group = "id",
                   varnames = list(ALEAF = "A",
                                   Tleaf = "Tleaf",
                                   Ci = "Ci",
                                   PPFD = "Qin",
                                   Rd = "resp"),
                   fitTPU = TRUE,
                   useRd = TRUE,
                   Tcorrect = FALSE)

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
plot(r4_hn_ni)
aci.r4_hn_ni <- data.frame(id = "r4_hn_ni",
                           t(coef(r4_hn_ni)))
head(aci.r4_hn_ni)

## Extract coefficients and separate id into rep, n.trt, and inoc.
## Also, merge separate r4_hn_ni coefficients
aci.coef <- aci.tpu %>%
  coef() %>%
  select(id, Vcmax, Jmax, Rd, TPU) %>%
  full_join(aci.r4_hn_ni) %>%
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
  data.frame() %>%
  left_join(leaf.area)

aci.coef

write.csv(aci.coef, "../data/trait_data.csv")


