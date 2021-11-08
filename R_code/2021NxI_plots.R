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
        axis.title = element_text(size = 16, face = "bold"),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 14, face = "bold"),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor = element_blank())

## Load data and add treatment column
data <- read.csv("../data/2021NxI_trait_data.csv",
                 na.strings = "NA")
data <- data %>%
  unite("treatment", n.trt:inoc, remove = "FALSE") %>%
  mutate(treatment = factor(treatment, levels = c("LN_NI", "HN_NI",
                                                  "LN_YI", "HN_YI")))

## Load compact letters data and add treatment column
comp.letters <- read.csv("../data/2021NxI_compact_letters.csv",
                         na.strings = "<NA>")

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
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "a400" & comparison == "full"),
            aes(y = 19, label = compact), 
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
  geom_jitter(width = 0.05, size = 2, alpha = 0.5, fill = "black") +
  geom_text(data = subset(comp.letters, 
                          variable == "a400" & comparison == "n.trt"),
            aes(y = 19, label = compact), 
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

a400.inoc <- ggplot(data = data, 
                    aes(x = factor(inoc,
                                   levels = c("NI", "YI")),
                        y = a, fill = inoc)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5, fill = "black") +
  geom_text(data = subset(comp.letters, 
                          variable == "a400" & comparison == "inoc"),
            aes(y = 19, label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("Not inoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 20),
                     breaks = seq(0, 20, 5)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme
a400.inoc

##########################################################################
## Vcmax25
##########################################################################
vcmax <- ggplot(data = data, 
                aes(x = treatment,
                    y = vcmax, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "vcmax25" & comparison == "full"),
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
  pubtheme +
  theme(axis.title.y = element_text(size = 16))
vcmax

vcmax.soil <- ggplot(data = data, 
                     aes(x = factor(n.trt,
                                    levels = c("LN", "HN")),
                         y = vcmax, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "vcmax25" & comparison == "n.trt"),
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
  pubtheme +
  theme(axis.title.y = element_text(16))
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
  geom_text(data = subset(comp.letters, 
                          variable == "jmax25" & comparison == "full"),
            aes(y = 125, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(30, 130),
                     breaks = seq(30, 130, 25)) +
  labs(x = "Soil nitrogen fertilization",
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
  geom_text(data = subset(comp.letters, 
                          variable == "jmax25" & comparison == "n.trt"),
            aes(y = 125, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(30, 130),
                     breaks = seq(30, 130, 25)) +
  labs(x = "Soil nitrogen fertilization",
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
  geom_text(data = subset(comp.letters, 
                          variable == "jmax25.vcmax25" & comparison == "full"),
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
  geom_text(data = subset(comp.letters, variable == "rd25" & comparison == "full"),
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
  geom_text(data = subset(comp.letters, 
                          variable == "rd25" & comparison == "n.trt"),
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
  geom_text(data = subset(comp.letters, 
                          variable == "rd.vcmax" & comparison == "full"),
            aes(y = 0.028, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 0.03),
                     breaks = seq(0, 0.03, 0.01)) +
  labs(x = "Soil nitrogen fertilization",
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
  geom_text(data = subset(comp.letters, 
                          variable == "rd25.vcmax25" & comparison == "n.trt"),
            aes(y = 0.028, 
                label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 0.03),
                     breaks = seq(0, 0.03, 0.01)) +
  labs(x = "Soil nitrogen fertilization",
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
  geom_text(data = subset(comp.letters, 
                          variable == "gs" & comparison == "full"),
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
  geom_text(data = subset(comp.letters,
                          variable == "gs" & comparison == "n.trt"),
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
  geom_text(data = subset(comp.letters, 
                          variable == "ci.ca" & comparison == "full"),
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
  geom_text(data = subset(comp.letters, 
                          variable == "ci.ca" & comparison == "n.trt"),
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
  geom_text(data = subset(comp.letters, 
                          variable == "iwue" & comparison == "full"),
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
  geom_text(data = subset(comp.letters, 
                          variable == "iwue" & comparison == "n.trt"),
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
  geom_text(data = subset(comp.letters,
                          variable == "vcmax.gs" & comparison == "full"),
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
  geom_text(data = subset(comp.letters,
                          variable == "vcmax.gs" & comparison == "n.trt"),
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
  geom_text(data = subset(comp.letters, variable == "sla" & comparison == "full"),
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
  geom_text(data = subset(comp.letters, variable == "sla" & comparison == "n.trt"),
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
  geom_text(data = subset(comp.letters, 
                          variable == "focal.area" & comparison == "full"),
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
  geom_text(data = subset(comp.letters, variable == "focal.area" & comparison == "full"),
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
  geom_text(data = subset(comp.letters,
                          variable == "focal.biomass" & comparison == "full"),
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
  geom_text(data = subset(comp.letters,
                          variable == "focal.biomass" & comparison == "n.trt"),
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
## Total leaf area
##########################################################################
tla <- ggplot(data = data, 
               aes(x = treatment,
                   y = total.leaf.area, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "total.leaf.area" & comparison == "full"),
            aes(y = 1600, label = compact), 
            fontface = "bold", size = 5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 1700),
                     breaks = seq(0, 1600, 400)) +
  labs(x = "Soil nitrogen fertilization",
       y = expression(bold("Total leaf area (cm"^"2"~")"))) +
  guides(fill = "none") +
  pubtheme
tla

tla.soil <- ggplot(data = data, 
                    aes(x = factor(n.trt,
                                   levels = c("LN", "HN")),
                        y = total.leaf.area, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5, fill = "black") +
  geom_text(data = subset(comp.letters, 
                          variable == "total.leaf.area" & comparison == "n.trt"),
            aes(y = 1600, label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 1700),
                     breaks = seq(0, 1600, 400)) +
  labs(x = "Soil nitrogen fertilization",
       y = NULL) +
  guides(fill = "none") +
  pubtheme
tla.soil

a400.inoc <- ggplot(data = data, 
                    aes(x = factor(inoc,
                                   levels = c("NI", "YI")),
                        y = total.leaf.area, fill = inoc)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5, fill = "black") +
  geom_text(data = subset(comp.letters, 
                          variable == "total.leaf.area" & comparison == "inoc"),
            aes(y = 1600, label = compact), 
            fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("Not inoculated", "Inoculated")) +
  scale_y_continuous(limits = c(0, 1700),
                     breaks = seq(0, 1600, 400)) +
  labs(x = NULL,
       y = NULL) +
  guides(fill = "none") +
  pubtheme
a400.inoc

##########################################################################
## Figure S1 (?): Leaf morphology
##########################################################################
a <- ggarrange(sla, sla.soil, fa, fa.soil, fbio, fbio.soil,
          nrow = 3, ncol = 2, common.legend = TRUE, widths = c(1.5,1),
          align = "hv", legend = "right", 
          labels = "AUTO", 
          font.label = list(size = 18, face = "bold"))
a
ggsave(filename = "../docs/figs/2021NxI_soy_figS1_leafMorph.png",
       a,
       width = 9,
       height = 10,
       units = "in",
       dpi = "retina")

##########################################################################
## Figure 1: Leaf photosynthesis and respiration
##########################################################################
b <- ggarrange(a400, a400.soil, rd, rd.soil, widths = c(1.5,1),
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
c <- ggarrange(vcmax, vcmax.soil, jmax, jmax.soil,
              nrow = 2, ncol = 2, common.legend = TRUE, 
              align = "hv", legend = "right", 
              labels = "AUTO", widths = c(1.5,1),
              font.label = list(size = 18, face = "bold"))
c

ggsave(filename = "../docs/figs/2021NxI_soy_fig2_leafBiochem.png",
       c,
       width = 9,
       height = 7,
       units = "in",
       dpi = "retina")

##########################################################################
## Figure 3: Rd:Vcmax
##########################################################################
d <- ggarrange(rd.vcmax, rd.vcmax.soil,
               nrow = 1, ncol = 2, common.legend = TRUE, 
               align = "hv", legend = "right", 
               labels = "AUTO", widths = c(1.5,1),
               font.label = list(size = 18, face = "bold"))
d

ggsave(filename = "../docs/figs/2021NxI_soy_fig3_rdvcmax.png",
       d,
       width = 9,
       height = 3.5,
       units = "in",
       dpi = "retina")
##########################################################################
## Figure 3: Water usage
##########################################################################
e <- ggarrange(gs, gs.soil, iwue, iwue.soil, cica, cica.soil,
               nrow = 3, ncol = 2, common.legend = TRUE, widths = c(1.5,1),
               align = "hv", legend = "right", 
               labels = "AUTO", 
               font.label = list(size = 18, face = "bold"))


e
ggsave(filename = "../docs/figs/2021NxI_soy_fig3_waterUsage.png",
       e,
       width = 9,
       height = 11,
       units = "in",
       dpi = "retina")

##########################################################################
## Figure 4: PNUE/WUE
##########################################################################


##########################################################################
## Figure 5: Whole plant measures
##########################################################################

f <- ggarrange(tla, tla.soil,
               nrow = 1, ncol = 2, common.legend = TRUE, 
               align = "hv", legend = "right", 
               labels = "AUTO", widths = c(1.5,1),
               font.label = list(size = 18, face = "bold"))
f

ggsave(filename = "../docs/figs/2021NxI_soy_fig5_totalLeafArea.png",
       f,
       width = 9,
       height = 3.5,
       units = "in",
       dpi = "retina")
