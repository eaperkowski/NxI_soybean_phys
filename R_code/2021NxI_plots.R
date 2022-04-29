## Libraries
library(tidyverse)
library(dplyr)
library(ggpubr)
library(patchwork)
library(gghalves)

## Central figure theme
pubtheme <- theme_bw() +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 18, face = "bold"),
        panel.border = element_rect(size = 3, fill = NA),
        axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 20, face = "bold"),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, face = "bold"),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor = element_blank())

## Load data and add treatment column
data <- read.csv("../data/2021NxI_trait_data.csv", na.strings = "NA")
data <- data %>%
  unite("treatment", n.trt:inoc, remove = "FALSE") %>%
  mutate(treatment = factor(treatment, levels = c("ln_ni", "ln_yi",
                                                  "hn_ni", "hn_yi")),
         nod.root.biomass = nodule.biomass / root.biomass)

## Load compact letters data and add treatment column
comp.letters <- read.csv("../data/2021NxI_compact_letters.csv",
                         na.strings = "<NA>")


## Add colorblind friendly palette
cbbPalette <- c("#DDAA33", "#BB5566", "#004488", "#BBBBBB")

# Remove outliers from Bonferroni tests
data$narea[39] <- NA
data$jmax25.vcmax25[c(46, 49)] <- NA
data$rd25[c(35, 63)] <- NA
data$rd25.vcmax25[c(35)] <- NA
data$ci.ca[c(23, 59)] <- NA
data$pnue[39] <- NA
data$narea.gs[39] <- NA

##########################################################################
## Nmass
##########################################################################
nmass <- ggplot(data = data, aes(x = treatment, y = leaf.n, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "nmass" & comparison == "full"),
            aes(y = 8, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(2, 8.4), breaks = seq(2, 8, 2)) +
  labs(x = NULL,
       y = expression(bold("N"["mass"]~"(g g"^"-1"~")")),
       fill = "Inoculation status") +
  pubtheme +
  theme(axis.text.x = element_blank())
nmass

##########################################################################
## SLA
##########################################################################
sla <- ggplot(data = data, aes(x = treatment, y = sla, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, variable == "sla" & comparison == "full"),
            aes(y = 750, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(450, 770), breaks = seq(450, 750, 100)) +
  labs(x = NULL,
       y = expression(bold("SLA (cm"^"-2"~"g"^"-1"~")")),
       fill = "Inoculation status") +
  pubtheme
sla

##########################################################################
## Narea
##########################################################################
narea <- ggplot(data = data, aes(x = treatment, y = narea, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "narea" & comparison == "full"),
            aes(y = 1.8, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0.6, 1.8), breaks = seq(0.6, 1.8, 0.4)) +
  labs(x = NULL,
       y = expression(bold("N"["area"]~"(g m"^"-2"~")")),
       fill = "Inoculation status") +
  pubtheme
narea

##########################################################################
## A400
##########################################################################
a400 <- ggplot(data = data, aes(x = treatment, y = anet, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "a400" & comparison == "full"),
            aes(y = 20, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = NULL,
       y = expression(bold("A"["net"]~"(μmol m"^"-2"~"s"^"-1"~")")),
       fill = "Inoculation status") +
  pubtheme +
  theme(axis.text.x = element_blank())
a400

##########################################################################
## Vcmax25
##########################################################################
vcmax <- ggplot(data = data, aes(x = treatment, y = vcmax25, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "vcmax25" & comparison == "full"),
            aes(y = 130, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(30, 130), breaks = seq(30, 130, 25)) +
  labs(x = NULL,
       y = expression(bold("V"["cmax25"]~"(μmol m"^"-2"~"s"^"-1"~")")),
       fill = "Inoculation status") +
  pubtheme
vcmax

##########################################################################
## Jmax25
##########################################################################
jmax <- ggplot(data = data, aes(x = treatment, y = jmax25, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "jmax25" & comparison == "full"),
            aes(y = 130, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(30, 130), breaks = seq(30, 130, 25)) +
  labs(x = NULL,
       y = expression(bold("J"["max25"]~"(μmol m"^"-2"~"s"^"-1"~")")),
       fill = "Inoculation status") +
  pubtheme
jmax

##########################################################################
## Jmax:Vcmax25
##########################################################################
jmax.vcmax <- ggplot(data = data, aes(x = treatment, 
                                      y = jmax25.vcmax25, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "jmax25.vcmax25" & comparison == "full"),
            aes(y = 1.6, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0.8, 1.6), breaks = seq(0.8, 1.6, 0.2)) +
  labs(x = NULL,
       y = expression(bold("J"["max25"]~": V"["cmax25"])),
       fill = "Inoculation status") +
  pubtheme
jmax.vcmax

##########################################################################
## Rd 
##########################################################################
rd <- ggplot(data = data, aes(x = treatment, y = rd25, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, variable == "rd25" & comparison == "full"),
            aes(y = 2, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  labs(x = NULL,
       y = expression(bold("R"["d25"]~"(μmol m"^"-2"~"s"^"-1"~")")),
       fill = "Inoculation status") +
  pubtheme +
  theme(axis.text.x = element_blank())
rd

##########################################################################
## Rd25:Vcmax25
##########################################################################
rd.vcmax <- ggplot(data = data, aes(x = treatment, y = rd25.vcmax25, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "rd.vcmax" & comparison == "full"),
            aes(y = 0.03, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 0.03), breaks = seq(0, 0.03, 0.01)) +
  labs(x = "Soil nitrogen fertilization",
       y = expression(bold("R"["d25"]~": V"["cmax25"])),
       fill = "Inoculation status") +
  pubtheme
rd.vcmax

##########################################################################
## Stomatal conductance
##########################################################################
gs <- ggplot(data = data, aes(x = treatment, y = gsw, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "gs" & comparison == "full"),
            aes(y = 0.45, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 0.45), breaks = seq(0, 0.45, 0.15)) +
  labs(x = NULL,
       y = expression(bold("g"["s400"]~"(μmol m"^"-2"~"s"^"-1"~")")),
       fill = "Inoculation status") +
  pubtheme
gs

##########################################################################
## Ci:Ca
##########################################################################
cica <- ggplot(data = data, aes(x = treatment, y = ci.ca, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "ci.ca" & comparison == "full"),
            aes(y = 1, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, 0.2)) +
  labs(x = NULL,
       y = expression(bold("C"["i"]~": C"["a"])),
       fill = "Inoculation status") +
  pubtheme

##########################################################################
## PNUE
##########################################################################
pnue <- ggplot(data = data, aes(x = treatment, y = pnue, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "pnue" & comparison == "full"),
            aes(y = 20, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = NULL,
       y = expression(bold("PNUE (μmol CO"["2"]~" gN"^"-1"~"s"^"-1"~")")),
       fill = "Inoculation status") +
  pubtheme +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 18))
pnue

##########################################################################
## iWUE
##########################################################################
iwue <- ggplot(data = data, aes(x = treatment, y = iwue, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "iwue" & comparison == "full"),
            aes(y = 160, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
  labs(x = NULL,
       y = expression(bold("iWUE (μmol CO"["2"]~"mol"^"-1"~"H"["2"]~"O)")),
       fill = "Inoculation status") +
  pubtheme +
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(size = 18))
iwue

##########################################################################
## Narea:gs
##########################################################################
narea.gs <- ggplot(data = data, aes(x = treatment, y = narea.gs, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters,
                          variable == "narea.gs" & comparison == "full"),
            aes(y = 20, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = NULL,
       y = expression(bold("N"["area"]~":g"["s"]~"(gN s mol"^"-1"~"H"["2"]~"O)")),
       fill = "Inoculation status") +
  pubtheme +
  theme(axis.title.y = element_text(size = 18))
narea.gs

##########################################################################
## Vcmax:gs
##########################################################################
vcmax.gs <- ggplot(data = data, aes(x = treatment, y = vcmax.gs, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters,
                          variable == "vcmax.gs" & comparison == "full"),
            aes(y = 1600, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 1600), breaks = seq(0, 1600, 400)) +
  labs(x = NULL,
       y = expression(bold("V"["cmax"]~":g"["s"]~"(μmol CO"["2"]~"mol"^"-1"~"H"["2"]~"O)")),
       fill = "Inoculation status") +
  pubtheme +
  theme(axis.title.y = element_text(size = 18))
vcmax.gs

##########################################################################
## Total leaf area
##########################################################################
tla <- ggplot(data = data, aes(x = treatment, y = total.leaf.area, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "total.leaf.area" & comparison == "full"),
            aes(y = 1800, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 1800), breaks = seq(0, 1800, 600)) +
  labs(x = NULL,
       y = expression(bold("Total leaf area (cm"^"2"~")")),
       fill = "Inoculation status") +
  pubtheme
tla

##########################################################################
## Total biomass
##########################################################################
tbio <- ggplot(data = data, aes(x = treatment, y = total.biomass, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, 
                          variable == "total.biomass" & comparison == "full"),
            aes(y = 12, label = compact), fontface = "bold", size = 7) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated", 
                                                    "Inoculated")) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 3)) +
  labs(x = NULL,
       y = "Whole plant biomass (g)",
       fill = "Inoculation status") +
  pubtheme
tbio

##########################################################################
## Structural carbon costs to acquire nitrogen
##########################################################################
ncost <- ggplot(data = data, aes(x = treatment, y = n.cost, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, variable == "ncost" & comparison == "full"),
            aes(y = 20, label = .group), fontface = "bold", size = 7) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm","630 ppm")) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated",
                                                    "Inoculated")) +
  labs(x = NULL,
       y = expression(bold("N"["cost"]~"(g C g"^"-1"~"N)")),
       fill = "Inoculation status") +
  pubtheme
ncost

##########################################################################
## Belowground carbon figure
##########################################################################
bgc.plot <- ggplot(data = data, aes(x = treatment, y = bg.total.c, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, variable == "bgc" & comparison == "full"),
            aes(y = 1.8, label = .group), fontface = "bold", size = 7) +
  scale_y_continuous(limits = c(0, 1.85), breaks = seq(0, 1.8, 0.6)) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm","630 ppm")) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated",
                                                    "Inoculated")) +
  labs(x = NULL,
       y = expression(bold("C"["bg"]~"(g C)")),
       fill = "Inoculation status") +
  pubtheme +
  theme(axis.text.x = element_blank())
bgc.plot

##########################################################################
## Whole plant nitrogen
##########################################################################
wpn.plot <- ggplot(data = data, aes(x = treatment, y = wp.total.n, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, variable == "wpn" & comparison == "full"),
            aes(y = 0.3, label = .group), fontface = "bold", size = 7) +
  scale_y_continuous(limits = c(0, 0.32), breaks = seq(0, 0.3, 0.1)) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated",
                                                    "Inoculated")) +
  labs(x = NULL,
       y = expression(bold("N"["wp"]~"(g N)")),
       fill = "Inoculation status") +
  pubtheme
wpn.plot

##########################################################################
## Root nodule biomass:root biomass figure
##########################################################################
nodroot.plot <- ggplot(data = data, aes(x = treatment, 
                                        y = nod.root.biomass, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, variable == "nodroot" & comparison == "full"),
            aes(y = 0.1, label = .group), fontface = "bold", size = 7) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.025)) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated",
                                                    "Inoculated")) +
  labs(x = NULL,
       y = "Nodule biomass: root biomass",
       fill = "Inoculation status") +
  pubtheme
nodroot.plot

##########################################################################
## Root nodule biomass figure
##########################################################################
nod.plot <- ggplot(data = data, aes(x = treatment, 
                                    y = nodule.biomass,
                                    fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, variable == "nod" & comparison == "full"),
            aes(y = 0.1, label = .group), fontface = "bold", size = 7) +
  scale_y_continuous(limits = c(0, 0.105), breaks = seq(0, 0.1, 0.025)) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated",
                                                    "Inoculated")) +
  labs(x = NULL,
       y = "Nodule biomass (g)",
       fill = "Inoculation status") +
  pubtheme +
  theme(axis.text.x = element_blank())
nod.plot

##########################################################################
## Root biomass figure
##########################################################################
root.plot <- ggplot(data = data, aes(x = treatment, 
                                     y = root.biomass,
                                     fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, variable == "root" & comparison == "full"),
            aes(y = 4, label = .group), fontface = "bold", size = 7) +
  scale_y_continuous(limits = c(0, 4.2), breaks = seq(0, 4.2, 1)) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630 ppm", "630 ppm")) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated",
                                                    "Inoculated")) +
  labs(x = NULL,
       y = "Root biomass (g)",
       fill = "Inoculation status") +
  pubtheme
root.plot

##########################################################################
## BVR figure
##########################################################################
bvr.plot <- ggplot(data = data, aes(x = treatment,
                                    y = bvr, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75) +
  geom_jitter(width = 0.05, size = 3, alpha = 0.5, show.legend = FALSE) +
  geom_text(data = subset(comp.letters, variable == "bvr" & comparison == "full"),
            aes(y = 2, label = .group), fontface = "bold", size = 7) +
  geom_hline(yintercept = 1, size = 1.5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  scale_x_discrete(labels = c("70 ppm", "70 ppm", "630ppm", "630 ppm")) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated",
                                                    "Inoculated")) +
  labs(x = "Soil nitrogen fertilization",
       y = expression(bold("Biomass: pot volume (g L"^"-1"~")")),
       fill = "Inoculation status") +
  pubtheme
bvr.plot

##########################################################################
## Figure 1: Leaf nitrogen allocation
##########################################################################
fig1 <- narea + (nmass / sla) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 20, face = "bold"))


ggarrange(fig1) %>%
  annotate_figure(bottom = text_grob(label = "Soil nitrogen fertilization",
                                     face = "bold", size = 20, hjust = 0.75)) %>%
  ggexport(filename = "../docs/figs/2021NxI_soy_fig1_leafN.png", 
           width = 9000, height = 4000, res = 600)

##########################################################################
## Figure 2: Leaf photosynthesis, respiration, Vcmax, Jmax
##########################################################################
fig2 <- a400 + rd + vcmax + jmax +
  plot_layout(ncol = 2, nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 20, face = "bold"))


ggarrange(fig2) %>%
  annotate_figure(bottom = text_grob(expression(
    bold("Soil nitrogen fertilization")), size = 20)) %>%
  ggexport(filename = "../docs/figs/2021NxI_soy_fig2_leafPhoto.png", 
           width = 9000, height = 5000, res = 600)

##########################################################################
## Figure 3: PNUE, iWUE (to be replaced with chi), Narea:gs, Vcmax:gs
##########################################################################
fig3 <- pnue + iwue + narea.gs + vcmax.gs + 
  plot_layout(ncol = 2, nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 20, face = "bold"))

ggarrange(fig3) %>%
  annotate_figure(bottom = text_grob(expression(
                    bold("Soil nitrogen fertilization")), size = 20)) %>%
  ggexport(filename = "../docs/figs/2021NxI_soy_fig3_pnueiwue.png", 
           width = 9000, height = 5200, res = 600)

##########################################################################
## Figure 4: Structural carbon costs to acquire nitrogen
##########################################################################
fig4 <- ncost + (bgc.plot / wpn.plot) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 20, face = "bold"))


ggarrange(fig4) %>%
  annotate_figure(bottom = text_grob(label = "Soil nitrogen fertilization",
                                     face = "bold", size = 20, hjust = 0.75)) %>%
  ggexport(filename = "../docs/figs/2021NxI_soy_fig4_ncost.png", 
           width = 9000, height = 4000, res = 600)

##########################################################################
## Figure 5: Whole plant measures
##########################################################################
fig5 <- tla + tbio +
  plot_layout(ncol = 2, nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 20, face = "bold"))


ggarrange(fig5) %>%
  annotate_figure(bottom = text_grob(expression(
    bold("Soil nitrogen fertilization")), size = 20)) %>%
  ggexport(filename = "../docs/figs/2021NxI_soy_fig4_tla.png", 
           width = 9000, height = 3000, res = 600)


##########################################################################
## Figure 6: Root nodulation
##########################################################################
fig6 <- nodroot.plot + (nod.plot / root.plot) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 20, face = "bold"))



ggarrange(fig6) %>%
  annotate_figure(bottom = text_grob(label = "Soil nitrogen fertilization",
                                     face = "bold", size = 20, hjust = 0.75)) %>%
  ggexport(filename = "../docs/figs/2021NxI_soy_fig6_nodulation.png", 
           width = 9000, height = 4000, res = 600)

##########################################################################
## Figure S1: BVR
##########################################################################
ggarrange(bvr.plot, nrow = 1, ncol = 1, common.legend = TRUE, align = "hv", 
          legend = "right", font.label = list(size = 18, face = "bold")) %>%
  ggexport(filename = "../docs/figs/2021NxI_soy_figS1_bvr.png", 
           width = 4000, height = 3000, res = 600)

