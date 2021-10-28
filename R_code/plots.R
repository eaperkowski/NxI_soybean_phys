## Libraries
library(tidyverse)
library(dplyr)
library(ggpubr)

## Central figure theme
pubtheme <- theme_bw() +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, face = "bold"),
        panel.border = element_rect(size = 3, fill = NA),
        axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 18, face = "bold"),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 14, face = "bold"),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor = element_blank())

## Load data and add treatment column
data <- read.csv("../data/trait_data.csv",
                 na.strings = "NA")
data <- data %>%
  unite("treatment", n.trt:inoc, remove = "FALSE") %>%
  mutate(treatment = factor(treatment, levels = c("LN_NI", "HN_NI",
                                                  "LN_YI", "HN_YI")))

## Load compact letters data and add treatment column
comp.letters.full <- read.csv("../data/comp.letters.full.csv")
comp.letters.full <- comp.letters.full %>%
  unite("treatment", n.trt:inoc, remove = "FALSE") %>%
  mutate(treatment = factor(treatment, levels = c("LN_NI", "HN_NI",
                                                  "LN_YI", "HN_YI")),
         compact = tolower(compact))

comp.letters.soiln <- read.csv("../data/comp.letters.soil.n.csv")
comp.letters.soiln <- comp.letters.soiln %>%
  unite("treatment", n.trt, remove = "FALSE") %>%
  mutate(treatment = factor(treatment, levels = c("LN", "HN")),
         compact = tolower(compact))

## Add colorblind friendly palette
cbbPalette <- c("#DDAA33", "#BB5566", "#004488", "#BBBBBB")

# Remove outliers from Bonferroni tests
data$jmax25.vcmax25[c(46, 49)] <- NA
data$rd25[c(35, 63)] <- NA
data$rd25.vcmax25[35] <- NA
data$ci.ca[c(23)] <- NA

##########################################################################
## A400
##########################################################################
a400 <- ggplot(data = data, 
                  aes(x = treatment,
                      y = a, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(aes(shape = inoc),
              width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.full, variable == "a400"),
            aes(y = 19, 
                label = compact), 
            fontface = "bold", size = 5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 20),
                     breaks = seq(0, 20, 5)) +
  labs(x = NULL,
       y = expression(bold("A"["400"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
a400

a400.soil <- ggplot(data = data, 
                    aes(x = factor(n.trt,
                                   levels = c("LN", "HN")),
                        y = a, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(aes(shape = inoc), width = 0.05, size = 3, 
              alpha = 0.25, fill = "black") +
  geom_text(data = subset(comp.letters.soiln, 
                          variable == "a400"),
            aes(x = n.trt, y = 19, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 20),
                     breaks = seq(0, 20, 5)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme
a400.soil

##########################################################################
## Vcmax25
##########################################################################
vcmax <- ggplot(data = data, 
                aes(x = treatment,
                    y = vcmax, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.full, variable == "vcmax25"),
            aes(y = 125, 
                label = compact), 
            fontface = "bold", size = 5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(30, 130),
                     breaks = seq(30, 130, 25)) +
  labs(x = NULL,
       y = expression(bold("V"["cmax25"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
vcmax

vcmax.soil <- ggplot(data = data, 
                     aes(x = factor(n.trt,
                                    levels = c("LN", "HN")),
                         y = vcmax, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.soiln, variable == "vcmax25"),
            aes(x = n.trt, y = 125, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(30, 130),
                     breaks = seq(30, 130, 25)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme
vcmax.soil
##########################################################################
## Jmax25
##########################################################################
jmax <- ggplot(data = data, aes(x = treatment,
                                y = jmax25, 
                                fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters.full, 
                          variable == "jmax25"),
            aes(y = 125, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(30, 130),
                     breaks = seq(30, 130, 25)) +
  labs(x = NULL,
       y = expression(bold("J"["max25"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
jmax

jmax.soil <- ggplot(data = data, aes(x = factor(n.trt,
                                              levels = c("LN", "HN")),
                                   y = jmax25, 
                                   fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.soiln, variable == "jmax25"),
            aes(y = 125, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(30, 130),
                     breaks = seq(30, 130, 25)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme
jmax.soil

##########################################################################
## Jmax:Vcmax25
##########################################################################
jmax.vcmax <- ggplot(data = data, aes(x = treatment,
                                      y = jmax25.vcmax25,
                                      fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters.full, variable == "jmax25.vcmax25"),
            aes(y = 1.55, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0.8, 1.6),
                     breaks = seq(0.8, 1.6, 0.2)) +
  labs(x = NULL,
       y = expression(bold("J"["max25"]~": V"["cmax25"]))) +
  guides(fill = "none") +
  pubtheme
jmax.vcmax

jmax.vcmax.soil <- ggplot(data = data,
                          aes(x = factor(n.trt,
                                         levels = c("LN", "HN")),
                              y = jmax25.vcmax25,
                              fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.soiln, variable == "jmax25.vcmax25"),
            aes(y = 1.55, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0.8, 1.6),
                     breaks = seq(0.8, 1.6, 0.2)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme
jmax.vcmax.soil

##########################################################################
## Rd 
##########################################################################
rd <- ggplot(data = data, aes(x = treatment,
                              y = rd25, 
                              fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters.full, variable == "rd25"),
            aes(y = 1.5, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 1.6),
                     breaks = seq(0, 1.6, 0.4)) +
  labs(x = "Soil N fertilization",
       y = expression(bold("R"["d25"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
rd

rd.soil <- ggplot(data = data, aes(x = factor(n.trt,
                                              levels = c("LN", "HN")),
                                   y = rd25, 
                                   fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.soiln, variable == "rd25"),
            aes(y = 1.5, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 1.6),
                     breaks = seq(0, 1.6, 0.4)) +
  labs(x = "Soil N fertilization",
       y = NULL) +
  guides(fill = "none") +
  pubtheme
rd.soil

##########################################################################
## Rd25:Vcmax25
##########################################################################
rd.vcmax <- ggplot(data = data, aes(x = treatment,
                                    y = rd25.vcmax25,
                                    fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters.full, variable == "rd.vcmax"),
            aes(y = 0.028, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 0.03),
                     breaks = seq(0, 0.03, 0.01)) +
  labs(x = NULL,
       y = expression(bold("R"["d25"]~": V"["cmax25"]))) +
  guides(fill = "none") +
  pubtheme
rd.vcmax

rd.vcmax.soil <- ggplot(data = data, aes(x = factor(n.trt, 
                                                    levels = c("LN", "HN")),
                                         y = rd25.vcmax25,
                                         fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters.soiln, variable == "rd25.vcmax25"),
            aes(y = 0.028, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 0.03),
                     breaks = seq(0, 0.03, 0.01)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme
rd.vcmax.soil

##########################################################################
## Gs
##########################################################################
gs <- ggplot(data = data, aes(x = treatment,
                                 y = gsw,
                                 fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters.full, variable == "gs"),
            aes(y = 0.425, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 0.45),
                     breaks = seq(0, 0.45, 0.15)) +
  labs(x = NULL,
       y = expression(bold("g"["s400"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
gs

gs.soil <- ggplot(data = data, aes(x = factor(n.trt, 
                                              levels = c("LN", "HN")),
                                   y = gsw,
                                   fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.soiln, variable == "gs"),
            aes(y = 0.425, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 0.45),
                     breaks = seq(0, 0.45, 0.15)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme
gs.soil 

##########################################################################
## Ci:Ca
##########################################################################
cica <- ggplot(data = data,  aes(x = treatment,
                                 y = ci.ca,
                                 fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters.full, variable == "ci.ca"),
            aes(y = 0.95, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0.4, 1),
                     breaks = seq(0.4, 1, 0.2)) +
  labs(x = NULL,
       y = expression(bold("C"["i"]~": C"["a"]))) +
  guides(fill = "none") +
  pubtheme
cica

cica.soil <- ggplot(data = data, aes(x = factor(n.trt, 
                                                levels = c("LN", "HN")),
                                     y = ci.ca,
                                     fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.soiln, variable == "ci.ca"),
            aes(y = 0.95, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0.4, 1),
                     breaks = seq(0.4, 1, 0.2)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme
cica.soil

##########################################################################
## iWUE
##########################################################################
iwue <- ggplot(data = data, aes(x = treatment,
                                y = iwue,
                                fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters.full, variable == "iwue"),
            aes(y = 150, 
                label = compact), 
            fontface = "bold", size = 5, hjust = 0.575) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 160),
                     breaks = seq(0, 160, 40)) +
  labs(x = NULL,
       y = expression(bold("A"["400"]~": g"["s400"]~"(μmol mol"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
iwue

iwue.soil <- ggplot(data = data, aes(x = factor(n.trt, 
                                                levels = c("LN", "HN")),
                                     y = iwue,
                                     fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.soiln, variable == "iwue"),
            aes(y = 150, 
                label = compact), 
            fontface = "bold", size = 5, hjust = 0.575) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm N", "630 ppm N")) +
  scale_y_continuous(limits = c(0, 160),
                     breaks = seq(0, 160, 40)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme
iwue.soil

##########################################################################
## Vcmax:gs
##########################################################################
vcmax.gs <- ggplot(data = data, aes(x = treatment,
                                    y = vcmax.gs,
                                    fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters.full, variable == "vcmax.gs"),
            aes(y = 1500, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 1600),
                     breaks = seq(0, 1600, 400)) +
  labs(x = NULL,
       y = expression(bold("V"["cmax"]~":g"["s400"]~"(mol mol"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
vcmax.gs

vcmax.gs.soil <- ggplot(data = data, aes(x = factor(n.trt, 
                                                    levels = c("LN", "HN")),
                                         y = vcmax.gs,
                                         fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.soiln, variable == "vcmax.gs"),
            aes(y = 1500, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm N", "630 ppm N")) +
  scale_y_continuous(limits = c(0, 1600),
                     breaks = seq(0, 1600, 400)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme
vcmax.gs.soil

##########################################################################
## SLA
##########################################################################
sla <- ggplot(data = data, aes(x = treatment,
                               y = sla,
                               fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters.full, variable == "sla"),
            aes(y = 7.75, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(4, 8),
                     breaks = seq(4, 8, 2)) +
  labs(x = NULL,
       y = expression(bold("SLA (cm"^"-2"~"g"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
sla

sla.soil <- ggplot(data = data, aes(x = factor(n.trt, 
                                               levels = c("LN", "HN")),
                                    y = sla, 
                                    fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.soiln, variable == "sla"),
            aes(y = 7.75, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm N", "630 ppm N")) +
  scale_y_continuous(limits = c(4, 8),
                     breaks = seq(4, 8, 2)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme
sla.soil

##########################################################################
## Focal area
##########################################################################
fa <- ggplot(data = data, aes(x = treatment,
                              y = focal.area,
                              fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters.full, variable == "focal.area"),
            aes(y = 85, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 90),
                     breaks = seq(0, 90, 30)) +
  labs(x = NULL,
       y = expression(bold("Focal area (cm"^"2"~")"))) +
  guides(fill = "none") +
  pubtheme
fa

fa.soil <- ggplot(data = data, aes(x = factor(n.trt, 
                                              levels = c("LN", "HN")),
                                   y = focal.area, 
                                   fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.soiln, variable == "focal.area"),
            aes(y = 85, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 90),
                     breaks = seq(0, 90, 30)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme +
  theme(strip.text = element_blank())
fa.soil

##########################################################################
## Focal biomass
##########################################################################
fbio <- ggplot(data = data, aes(x = treatment, 
                                y = dry.biomass,
                                fill = n.trt)) +
  geom_boxplot() +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.full, variable == "focal.biomass"),
            aes(y = 0.145, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0.05, 0.15),
                     breaks = seq(0.05, 0.15, 0.05)) +
  labs(x = "Soil N fertilization",
       y = "Dry leaf biomass (g)") +
  guides(fill = "none") +
  pubtheme
fbio

fbio.soil <- ggplot(data = data, aes(x = factor(n.trt, 
                                                levels = c("LN", "HN")),
                                     y = dry.biomass, 
                                     fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters.soiln, variable == "focal.biomass"),
            aes(y = 0.145, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0.05, 0.15),
                     breaks = seq(0.05, 0.15, 0.05)) +
  labs(x = "Soil N fertilization",
       y = NULL) +
  guides(fill = "none") +
  pubtheme
fbio.soil

##########################################################################
## Figure S1 (?): Leaf morphology
##########################################################################
a <- ggarrange(sla, sla.soil, fa, fa.soil, fbio, fbio.soil,
          nrow = 3, ncol = 2, common.legend = TRUE, widths = c(1.5,1),
          align = "hv", legend = "right", 
          labels = "AUTO", 
          font.label = list(size = 18, face = "bold"))
a
ggsave(filename = "../docs/figs/2021NxI_soy_fig1_leafMorph.png",
       a,
       width = 9,
       height = 10,
       units = "in",
       dpi = "retina")

##########################################################################
## Figure 1: Leaf photosynthesis and respiration
##########################################################################
b <- ggarrange(a400.soil, vcmax.soil, jmax.soil, rd.soil,
               nrow = 2, ncol = 2, common.legend = TRUE, 
               align = "hv", labels = "AUTO", 
               font.label = list(size = 18, face = "bold"))
b
ggsave(filename = "../docs/figs/2021NxI_soy_fig1_leafPhotoResp.png",
       b,
       width = 9,
       height = 7,
       units = "in",
       dpi = "retina")

##########################################################################
## Figure 2: Leaf biochemistry
##########################################################################
c <- ggarrange(vcmax.soil, jmax.soil,
              jmax.vcmax.soil, 
              nrow = 3, ncol = 2, common.legend = TRUE, 
              align = "hv", legend = "right", 
              labels = "AUTO", widths = c(1.5, 1),
              font.label = list(size = 18, face = "bold"))
c

ggsave(filename = "../docs/figs/2021NxI_soy_fig2_leafBiochem.png",
       c,
       width = 9,
       height = 12,
       units = "in",
       dpi = "retina")

##########################################################################
## Figure 3: Water usage
##########################################################################
d <- 



##########################################################################
## Figure 4: PNUE/WUE
##########################################################################



##########################################################################
## Figure 5: Whole plant measures
##########################################################################