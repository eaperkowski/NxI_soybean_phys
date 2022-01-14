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
        axis.title = element_text(size = 15, face = "bold"),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 14, face = "bold"),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor = element_blank())

## Load data and add treatment column
data <- read.csv("../data/2021NxI_trait_data.csv", na.strings = "NA")
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
## Nmass
##########################################################################
nmass <- ggplot(data = data, aes(x = treatment, y = leaf.n, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "nmass" & comparison == "full"),
            aes(y = 9, label = compact), fontface = "bold", size = 5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
  labs(x = NULL,
       y = expression(bold("N"["mass"]~"(g g"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
nmass

##########################################################################
## SLA
##########################################################################
sla <- ggplot(data = data, aes(x = treatment, y = sla, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters, variable == "sla" & comparison == "full"),
            aes(y = 7.75, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(4, 8), breaks = seq(4, 8, 1)) +
  labs(x = NULL,
       y = expression(bold("SLA (cm"^"-2"~"g"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
sla

##########################################################################
## Narea
##########################################################################
narea <- ggplot(data = data, aes(x = treatment, y = narea, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "narea" & comparison == "full"),
            aes(y = 1.85, label = compact), fontface = "bold", size = 5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  labs(x = "Soil nitrogen fertilization",
       y = expression(bold("N"["area"]~"(g m"^"-2"~")"))) +
  guides(fill = "none") +
  pubtheme
narea

##########################################################################
## A400
##########################################################################
a400 <- ggplot(data = data, aes(x = treatment, y = a, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "a400" & comparison == "full"),
            aes(y = 19, label = compact), fontface = "bold", size = 5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = NULL,
       y = expression(bold("A"["400"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
a400

##########################################################################
## Vcmax25
##########################################################################
vcmax <- ggplot(data = data, aes(x = treatment, y = vcmax, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "vcmax25" & comparison == "full"),
            aes(y = 125, label = compact), fontface = "bold", size = 5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(30, 130), breaks = seq(30, 130, 25)) +
  labs(x = NULL,
       y = expression(bold("V"["cmax25"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
vcmax

##########################################################################
## Jmax25
##########################################################################
jmax <- ggplot(data = data, aes(x = treatment, y = jmax25, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "jmax25" & comparison == "full"),
            aes(y = 125, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(30, 130), breaks = seq(30, 130, 25)) +
  labs(x = NULL,
       y = expression(bold("J"["max25"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
jmax

##########################################################################
## Jmax:Vcmax25
##########################################################################
jmax.vcmax <- ggplot(data = data, aes(x = treatment, y = jmax25.vcmax25, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "jmax25.vcmax25" & comparison == "full"),
            aes(y = 1.55, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0.8, 1.6), breaks = seq(0.8, 1.6, 0.2)) +
  labs(x = NULL,
       y = expression(bold("J"["max25"]~": V"["cmax25"]))) +
  guides(fill = "none") +
  pubtheme
jmax.vcmax

##########################################################################
## Rd 
##########################################################################
rd <- ggplot(data = data, aes(x = treatment, y = rd25, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters, variable == "rd25" & comparison == "full"),
            aes(y = 1.5, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 1.6), breaks = seq(0, 1.6, 0.4)) +
  labs(x = NULL,
       y = expression(bold("R"["d25"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
rd

##########################################################################
## Rd25:Vcmax25
##########################################################################
rd.vcmax <- ggplot(data = data, aes(x = treatment, y = rd25.vcmax25, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "rd.vcmax" & comparison == "full"),
            aes(y = 0.028, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 0.03), breaks = seq(0, 0.03, 0.01)) +
  labs(x = "Soil nitrogen fertilization",
       y = expression(bold("R"["d25"]~": V"["cmax25"]))) +
  guides(fill = "none") +
  pubtheme
rd.vcmax

##########################################################################
## Stomatal conductance
##########################################################################
gs <- ggplot(data = data, aes(x = treatment, y = gsw, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "gs" & comparison == "full"),
            aes(y = 0.425, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 0.45), breaks = seq(0, 0.45, 0.15)) +
  labs(x = NULL,
       y = expression(bold("g"["s400"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
gs

##########################################################################
## Ci:Ca
##########################################################################
cica <- ggplot(data = data, aes(x = treatment,y = ci.ca, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "ci.ca" & comparison == "full"),
            aes(y = 0.95, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, 0.2)) +
  labs(x = NULL,
       y = expression(bold("C"["i"]~": C"["a"]))) +
  guides(fill = "none") +
  pubtheme
cica

##########################################################################
## PNUE
##########################################################################
pnue <- ggplot(data = data, aes(x = treatment, y = pnue, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "pnue" & comparison == "full"),
            aes(y = 19, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = NULL,
       y = expression(bold("PNUE (μmol CO"["2"]~" gN"^"-1"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
pnue

##########################################################################
## iWUE
##########################################################################
iwue <- ggplot(data = data, aes(x = treatment, y = iwue, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "iwue" & comparison == "full"),
            aes(y = 150, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40)) +
  labs(x = NULL,
       y = expression(bold("A"["400"]~": g"["s400"]~"(μmol mol"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme
iwue

##########################################################################
## Narea:gs
##########################################################################
narea.gs <- ggplot(data = data, aes(x = treatment, y = narea.gs, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters,
                          variable == "narea.gs" & comparison == "full"),
            aes(y = 19, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(x = NULL,
       y = expression(bold("N"["area"]~":g"["s"]~"(gN s mol"^"-1"~"H"["2"]~"O)"))) +
  guides(fill = "none") +
  pubtheme
narea.gs

##########################################################################
## Vcmax:gs
##########################################################################
vcmax.gs <- ggplot(data = data, aes(x = treatment, y = vcmax.gs, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters,
                          variable == "vcmax.gs" & comparison == "full"),
            aes(y = 1550, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 1600), breaks = seq(0, 1600, 400)) +
  labs(x = NULL,
       y = expression(bold("V"["cmax"]~":g"["s"]~"(μmol CO"["2"]~"mol"^"-1"~"H"["2"]~"O)"))) +
  guides(fill = "none") +
  pubtheme
vcmax.gs

##########################################################################
## Focal area
##########################################################################
fa <- ggplot(data = data, aes(x = treatment, y = focal.area, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "focal.area" & comparison == "full"),
            aes(y = 85, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, 30)) +
  labs(x = NULL,
       y = expression(bold("Focal area (cm"^"2"~")"))) +
  guides(fill = "none") +
  pubtheme
fa

##########################################################################
## Focal biomass
##########################################################################
fbio <- ggplot(data = data, aes(x = treatment, y = dry.biomass, fill = n.trt)) +
  geom_boxplot() +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters,
                          variable == "focal.biomass" & comparison == "full"),
            aes(y = 0.145, label = compact), fontface = "bold", size = 5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0.05, 0.15), breaks = seq(0.05, 0.15, 0.05)) +
  labs(x = "Soil N fertilization",
       y = "Dry leaf biomass (g)") +
  guides(fill = "none") +
  pubtheme
fbio

##########################################################################
## Total leaf area
##########################################################################
tla <- ggplot(data = data, aes(x = treatment, y = total.leaf.area, fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, 
                          variable == "total.leaf.area" & comparison == "full"),
            aes(y = 1600, label = compact), fontface = "bold", size = 5) +
  geom_vline(xintercept = 2.5, linetype = "dashed", size = 1.5) +
  scale_fill_manual(values = cbbPalette) +
  scale_x_discrete(labels = c("70 ppm", "630 ppm", "70 ppm", "630 ppm")) +
  scale_y_continuous(limits = c(0, 1700), breaks = seq(0, 1600, 400)) +
  labs(x = "Soil nitrogen fertilization",
       y = expression(bold("Total leaf area (cm"^"2"~")"))) +
  guides(fill = "none") +
  pubtheme
tla

##########################################################################
## Figure 1: Leaf nitrogen allocation
##########################################################################

ggarrange(nmass, sla, narea, ncol = 1, nrow = 3, common.legend = TRUE,
          align = "hv", legend = "right", labels = "AUTO",
          font.label = list(size = 18, face = "bold")) %>%
  ggexport(filename = "../docs/figs/2021NxI_soy_fig1_leafN.png", 
           width = 3000, height = 5000, res = 600)

##########################################################################
## Figure 2: Leaf photosynthesis, respiration, Vcmax, Jmax
##########################################################################
ggarrange(a400, rd, vcmax, jmax, nrow = 2, ncol = 2, common.legend = TRUE, 
          align = "hv", labels = "AUTO", font.label = list(size = 18, 
                                                           face = "bold")) %>%
  annotate_figure(bottom = text_grob(expression(
    bold("Soil nitrogen fertilization")), size = 15)) %>%
  ggexport(filename = "../docs/figs/2021NxI_soy_fig2_leafPhoto.png", 
           width = 5500, height = 5000, res = 600)

##########################################################################
## Figure 3: PNUE, iWUE (to be replaced with chi), Narea:gs, Vcmax:gs
##########################################################################
ggarrange(pnue, iwue, narea.gs, vcmax.gs, nrow = 2, ncol = 2, 
          common.legend = TRUE, align = "hv", legend = "right", 
          labels = "AUTO", font.label = list(size = 18, face = "bold")) %>%
  annotate_figure(bottom = text_grob(expression(
                    bold("Soil nitrogen fertilization")), size = 15)) %>%
  ggexport(filename = "../docs/figs/2021NxI_soy_fig3_pnueiwue.png", 
           width = 5500, height = 5000, res = 600)

##########################################################################
## Figure 4: Whole plant measures
##########################################################################
ggarrange(tla, nrow = 1, ncol = 1, common.legend = TRUE, align = "hv", 
          legend = "right", labels = "AUTO", 
          font.label = list(size = 18, face = "bold")) %>%
  ggexport(filename = "../docs/figs/2021NxI_soy_fig4_tla.png", 
           width = 3500, height = 3000, res = 600)
