## Libraries
library(tidyverse)
library(dplyr)
library(ggpubr)

## Central figure theme
pubtheme <- theme_bw() +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        panel.border = element_rect(size = 3, fill = NA),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18, face = "bold"),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 14, face = "bold"),
        axis.ticks.length = unit(0.25, "cm"))

## Load data
data <- read.csv("../data/trait_data.csv",
                 na.strings = "NA")

## Add colorblind friendly palette
cbbPalette <- c("#DDAA33", "#BB5566", "#004488", "#BBBBBB")

# Remove outliers from Bonferroni tests
data$Jmax25.Vcmax25[c(17, 62)] <- NA
data$Rd.TPU[c(2, 7, 31)] <- NA
data$Rd.Vcmax[c(7, 31)] <- NA

data <- data %>%
  unite("treatment", n.trt:inoc, remove = FALSE)

facet.labels <- c("Not inoculated", "Inoculated")
names(facet.labels) <- c("NI", "YI")

##########################################################################
## A400
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = A,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(5, 20),
                     breaks = seq(5, 20, 5)) +
  labs(x = "Treatment",
       y = expression(bold("A"["400"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Vcmax25 with TPU
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = Vcmax25.TPU,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(30, 110),
                     breaks = seq(30, 110, 20)) +
  labs(x = "Treatment",
       y = expression(bold("V"["cmax25"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Jmax25 with TPU
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = Jmax25.TPU,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(30, 110),
                     breaks = seq(30, 110, 20)) +
  labs(x = "Treatment",
       y = expression(bold("J"["max25"]~"μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Jmax:Vcmax25 with TPU
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = Jmax25.Vcmax25,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0.8, 1.6),
                     breaks = seq(0.8, 1.6, 0.2)) +
  labs(x = "Treatment",
       y = expression(bold("J"["max25"]~": V"["cmax25"]))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Rd 
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = Rd.TPU,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0, 1.5),
                     breaks = seq(0, 1.5, 0.5)) +
  labs(x = "Treatment",
       y = expression(bold("R"["d"]~": V"["cmax25"]))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Rd:Vcmax (Vcmax is not standardized since Rd is not temp standardized)
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = Rd.Vcmax,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0, 0.025),
                     breaks = seq(0, 0.025, 0.005)) +
  labs(x = "Treatment",
       y = expression(bold("R"["d"]~": V"["cmax25"]))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Gs
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = gsw,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0, 0.4),
                     breaks = seq(0, 0.4, 0.1)) +
  labs(x = "Treatment",
       y = expression(bold("R"["d"]~": V"["cmax25"]))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Ci:Ca
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = ci.ca,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0.4, 0.9),
                     breaks = seq(0.4, 0.9, 0.1)) +
  labs(x = "Treatment",
       y = expression(bold("R"["d"]~": V"["cmax25"]))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## iWUE
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = iwue,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0, 150),
                     breaks = seq(0, 150, 50)) +
  labs(x = "Treatment",
       y = expression(bold("R"["d"]~": V"["cmax25"]))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Vcmax:gs
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = Vcmax.gs,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0, 1500),
                     breaks = seq(0, 1500, 500)) +
  labs(x = "Treatment",
       y = expression(bold("R"["d"]~": V"["cmax25"]))) +
  guides(fill = "none") +
  pubtheme




