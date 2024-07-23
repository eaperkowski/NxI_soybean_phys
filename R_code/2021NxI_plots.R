## Libraries
library(tidyverse)
library(dplyr)
library(ggpubr)
library(patchwork)

## Central figure theme
pubtheme <- theme_bw(base_size = 18) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(size = 1.5, fill = NA),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.title = element_text(face = "bold"),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor.y = element_blank(),
        legend.text.align = 0)

## Load data and add treatment column
data <- read.csv("../data/2021NxI_trait_data.csv", na.strings = "NA") %>%
  unite("treatment", n.trt:inoc, remove = FALSE) %>%
  filter(inoc == "yi" | (inoc == "ni" & nodule.biomass == 0)) %>%
  mutate(treatment = factor(treatment, levels = c("70_ni", "70_yi",
                                                  "630_ni", "630_yi")),
         n.trt = factor(n.trt, levels = c(70, 630)))

## Load compact letters data and add treatment column
comp.letters <- read.csv("../data/2021NxI_compact_letters.csv",
                         na.strings = "<NA>") %>%
  mutate(treatment = factor(treatment, levels = c("70_ni", "70_yi",
                                                  "630_ni", "630_yi")),
         n.trt = factor(n.trt, levels = c(70, 630)))

## Add colorblind friendly palette
cbbPalette <- c("#DDAA33", "#BB5566", "#004488", "#BBBBBB")

##########################################################################
## Carbon costs to acquire nitrogen
##########################################################################
ncost.plot <- ggplot(data = data, 
                     aes(x = n.trt, y = n.cost, fill = inoc)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(size = 0.75, outlier.shape = NA) +
  geom_point(size = 3, alpha = 0.5, show.legend = FALSE, shape = 21,
              position = position_jitterdodge(jitter.width = 0.05, 
                                              dodge.width = 0.75)) +
  geom_text(data = subset(comp.letters, variable == "ncost"),
            aes(y = 12, label = .group), 
            fontface = "bold", size = 7,
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(0, 12), breaks = seq(0, 12, 4)) +
  scale_x_discrete(labels = c("low N", "high N")) +
  scale_fill_manual(values = cbbPalette, labels = c("Not inoculated",
                                                    "Inoculated")) +
  labs(x = "Nitrogen fertilization treatment",
       y = expression(bold("Carbon cost to acquire nitrogen (gC gN"^"-1"~")")),
       fill = "Inoculation status") +
  pubtheme + theme(axis.title.y = element_text(size = 14))
ncost.plot

##########################################################################
## Belowground carbon figure
##########################################################################
bgc.plot <- ggplot(data = data, aes(x = n.trt, y = bg.total.c, fill = inoc)) +
  stat_boxplot(linewidth = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(size = 0.75, outlier.shape = NA) +
  geom_point(size = 3, alpha = 0.5, show.legend = FALSE, shape = 21,
             position = position_jitterdodge(jitter.width = 0.05, 
                                             dodge.width = 0.75)) +
  geom_text(data = subset(comp.letters, variable == "bgc"),
            aes(y = 1.6, label = .group), 
            fontface = "bold", size = 7,
            position = position_dodge(width = 0.75)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = seq(0, 1.6, 0.4)) +
  scale_x_discrete(labels = c("low N", "high N")) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("Not inoculated", "Inoculated")) +
  labs(x = "Nitrogen fertilization treatment",
       y = expression(bold("Belowground biomass carbon (gC)")),
       fill = "Inoculation status") +
  pubtheme + theme(axis.title.y = element_text(size = 14))
bgc.plot

##########################################################################
## Whole plant nitrogen
##########################################################################
wpn.plot <- ggplot(data = data, aes(x = n.trt, y = wp.total.n, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(size = 0.75, outlier.shape = NA) +
  geom_point(size = 3, alpha = 0.5, show.legend = FALSE, shape = 21,
              position = position_jitterdodge(jitter.width = 0.05, 
                                              dodge.width = 0.75)) +
  geom_text(data = subset(comp.letters, variable == "wpn"),
            aes(y = 0.28, label = .group), 
            fontface = "bold", size = 7,
            position = position_dodge(0.75)) +
  scale_y_continuous(limits = c(0, 0.28), breaks = seq(0, 0.28, 0.07)) +
  scale_x_discrete(labels = c("low N", "high N")) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("Not inoculated", "Inoculated")) +
  labs(x = "Nitrogen fertilization treatment",
       y = expression(bold("Total nitrogen biomass (gN)")),
       fill = "Inoculation status") +
  pubtheme + theme(axis.title.y = element_text(size = 14))
wpn.plot

##########################################################################
## Total leaf area
##########################################################################
tla.plot <- ggplot(data = data, aes(x = n.trt, y = total.leaf.area, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(size = 0.75, outlier.shape = NA) +
  geom_point(size = 3, alpha = 0.5, show.legend = FALSE, shape = 21,
             position = position_jitterdodge(jitter.width = 0.05, 
                                             dodge.width = 0.75)) +
  geom_text(data = subset(comp.letters, 
                          variable == "total.leaf.area"),
            aes(y = 1500, label = .group), fontface = "bold", size = 7,
            position = position_dodge(0.75)) +
  scale_y_continuous(limits = c(300, 1500), breaks = seq(300, 1500, 300)) +
  scale_x_discrete(labels = c("low N", "high N")) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("Not inoculated", "Inoculated")) +
  labs(x = "Nitrogen fertilization treatment",
       y = expression(bold("Total leaf area (cm"^"2"~")")),
       fill = "Inoculation status") +
  pubtheme
tla.plot

##########################################################################
## Total biomass
##########################################################################
tbio.plot <- ggplot(data = data, aes(x = n.trt, y = total.biomass, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(size = 0.75, outlier.shape = NA) +
  geom_point(size = 3, alpha = 0.5, show.legend = FALSE, shape = 21,
             position = position_jitterdodge(jitter.width = 0.05, 
                                             dodge.width = 0.75)) +
  geom_text(data = subset(comp.letters, 
                          variable == "total.biomass"),
            aes(y = 8, label = .group), fontface = "bold", size = 7,
            position = position_dodge(0.75)) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2)) +
  scale_x_discrete(labels = c("low N", "high N")) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("Not inoculated", "Inoculated")) +
  labs(x = "Nitrogen fertilization treatment",
       y = "Total biomass (g)",
       fill = "Inoculation status") +
  pubtheme
tbio.plot

##########################################################################
## Root nodule biomass:root biomass figure
##########################################################################
nodroot.plot <- ggplot(data = subset(data, inoc == "yi"), 
                       aes(x = n.trt, y = nod.root.biomass)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75, outlier.shape = NA, fill = cbbPalette[2], width = 0.375) +
  geom_point(size = 3, alpha = 0.5, show.legend = FALSE, shape = 21,
             position = position_jitter(width = 0.05),
             fill = cbbPalette[2]) +
  geom_text(data = subset(comp.letters, variable == "nodroot"),
            aes(y = 0.1, label = .group), fontface = "bold", size = 7) +
  scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, 0.025)) +
  scale_x_discrete(labels = c("low N", "high N")) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("Not inoculated", "Inoculated")) +
  labs(x = "Nitrogen fertilization treatment",
       y = "Nodule biomass: root biomass",
       fill = "Inoculation status") +
  pubtheme
nodroot.plot

##########################################################################
## Root nodule biomass figure
##########################################################################
nod.plot <- ggplot(data = subset(data, inoc == "yi"), 
                   aes(x = n.trt, y = nodule.biomass)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25) +
  geom_boxplot(size = 0.75, outlier.shape = NA, fill = cbbPalette[2], width = 0.375) +
  geom_point(size = 3, alpha = 0.5, show.legend = FALSE, shape = 21,
             position = position_jitter(width = 0.05),
             fill = cbbPalette[2]) +
  geom_text(data = subset(comp.letters, variable == "nod"),
            aes(y = 0.08, label = .group), fontface = "bold", size = 7) +
  scale_y_continuous(limits = c(0, 0.08), breaks = seq(0, 0.08, 0.02)) +
  scale_x_discrete(labels = c("low N", "high N")) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("Not inoculated", "Inoculated")) +
  labs(x = "Nitrogen fertilization treatment",
       y = "Nodule biomass (g)",
       fill = "Inoculation status") +
  pubtheme
nod.plot

##########################################################################
## Root biomass figure
##########################################################################
root.plot <- ggplot(data = data, aes(x = n.trt, y = root.biomass, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(size = 0.75, outlier.shape = NA) +
  geom_point(size = 3, alpha = 0.5, show.legend = FALSE, shape = 21,
             position = position_jitterdodge(jitter.width = 0.05, 
                                             dodge.width = 0.75)) +
  geom_text(data = subset(comp.letters, variable == "root"),
            aes(y = 3, label = .group), fontface = "bold", size = 7,
            position = position_dodge(0.75)) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_x_discrete(labels = c("low N", "high N")) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("Not inoculated", "Inoculated")) +
  labs(x = "Nitrogen fertilization treatment",
       y = "Root biomass (g)",
       fill = "Inoculation status") +
  pubtheme
root.plot

##########################################################################
## BVR figure
##########################################################################
bvr.plot <- ggplot(data = data, aes(x = n.trt, y = bvr, fill = inoc)) +
  stat_boxplot(size = 0.75, geom = "errorbar", width = 0.25, 
               position = position_dodge(width = 0.75)) +
  geom_boxplot(size = 0.75, outlier.shape = NA) +
  geom_point(size = 3, alpha = 0.5, show.legend = FALSE, shape = 21,
             position = position_jitterdodge(jitter.width = 0.05, 
                                             dodge.width = 0.75)) +
  geom_text(data = subset(comp.letters, variable == "bvr"),
            aes(y = 2, label = .group), fontface = "bold", size = 7) +
  geom_hline(yintercept = 1, size = 1.5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 2), breaks = seq(0, 2, 0.5)) +
  scale_x_discrete(labels = c("low N", "high N")) +
  scale_fill_manual(values = cbbPalette, 
                    labels = c("Not inoculated", "Inoculated")) +
  labs(x = "Nitrogen fertilization treatment",
       y = expression(bold("Biomass: pot volume (g L"^"-1"~")")),
       fill = "Inoculation status") +
  pubtheme
bvr.plot

##########################################################################
## Figure 1: Structural carbon costs to acquire nitrogen
##########################################################################
png("../docs/figs/2021NxI_fig1_ncost.png",
    width = 13, height = 10, units = 'in', res = 600)
ggarrange(ncost.plot, bgc.plot, wpn.plot, ncol = 2, nrow = 2,
          align = "hv", common.legend = T, legend = "right",
          labels = "AUTO", vjust = 1,
          font.label = list(size = 18, face = "bold"))
dev.off()

##########################################################################
## Figure 2: Whole-plant growth
##########################################################################
png("../docs/figs/2021NxI_fig2_tla.png",
    width = 12, height = 4.5, units = 'in', res = 600)
ggarrange(tla.plot, tbio.plot, ncol = 2, 
          common.legend = TRUE, align = "hv", legend = "right", 
          labels = "AUTO",
          font.label = list(size = 18, face = "bold"), widths = 1)
dev.off()

##########################################################################
## Figure 3: N fixation
##########################################################################
png("../docs/figs/2021NxI_fig3_nfix.png",
    width = 13, height = 10, units = 'in', res = 600)
ggarrange(nodroot.plot, nod.plot, root.plot, ncol = 2, nrow = 2,
          align = "hv", common.legend = T, legend = "right",
          labels = "AUTO", vjust = 1,
          font.label = list(size = 18, face = "bold"))
dev.off()

##########################################################################
## Figure S1: BVR
##########################################################################
png("../docs/figs/2021NxI_figS1_bvr.png",
    width = 9, height = 6, units = 'in', res = 600)
bvr.plot
dev.off()


