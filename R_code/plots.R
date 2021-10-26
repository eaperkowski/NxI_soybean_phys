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
comp.letters <- read.csv("../data/comp.letters.csv", strip.white = TRUE)

## Add colorblind friendly palette
cbbPalette <- c("#DDAA33", "#BB5566", "#004488", "#BBBBBB")

# Remove outliers from Bonferroni tests
data$jmax25.vcmax25[c(46, 49)] <- NA
data$rd25[c(35, 63)] <- NA
data$rd25.vcmax25[35] <- NA
data$ci.ca[c(23)] <- NA

# Add facet labels
facet.labels <- c("Not inoculated", "Inoculated")
names(facet.labels) <- c("NI", "YI")

##########################################################################
## A400
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = a,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(5, 20),
                     breaks = seq(5, 20, 5)) +
  geom_text(data = subset(comp.letters, variable == "a400"),
            aes(y = 19, 
                label = compact), 
            fontface = "bold", size = 5) +
  labs(x = "Treatment",
       y = expression(bold("A"["400"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Vcmax25 with TPU
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = vcmax25,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2) +
  geom_text(data = subset(comp.letters, variable == "vcmax25"),
            aes(y = 105, 
                label = compact), 
            fontface = "bold", size = 5) +
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
                        y = jmax25,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2) +
  geom_text(data = subset(comp.letters, variable == "jmax25"),
            aes(y = 110, 
                label = compact), 
            fontface = "bold", size = 5) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(30, 110),
                     breaks = seq(30, 110, 20)) +
  labs(x = "Treatment",
       y = expression(bold("J"["max25"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Jmax:Vcmax25 with TPU
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, 
                                   levels = c("LN", "HN")),
                        y = jmax25.vcmax25,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, variable == "a400"),
            aes(y = 1.55, 
                label = compact), 
            fontface = "bold", size = 5) +
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
                        y = rd25,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, variable == "rd"),
            aes(y = 1.5, 
                label = compact), 
            fontface = "bold", size = 5) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0, 1.6),
                     breaks = seq(0, 1.6, 0.4)) +
  labs(x = "Treatment",
       y = expression(bold("R"["d"]~"(μmol m"^"-2"~"s"^"-1"~")"))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Rd:Vcmax (Vcmax is not standardized since Rd is not temp standardized)
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = rd25.vcmax25,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, variable == "rd.vcmax"),
            aes(y = 0.028, 
                label = compact), 
            fontface = "bold", size = 5) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0, 0.03),
                     breaks = seq(0, 0.03, 0.01)) +
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
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, variable == "gsw"),
            aes(y = 0.425, 
                label = compact), 
            fontface = "bold", size = 5) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0, 0.45),
                     breaks = seq(0, 0.45, 0.15)) +
  labs(x = "Treatment",
       y = expression(bold("g"["s"]~"(mol mol"["-1"]~")"))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Ci:Ca
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = ci.ca,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, variable == "ci.ca"),
            aes(y = 0.95, 
                label = compact), 
            fontface = "bold", size = 5) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0.4, 1),
                     breaks = seq(0.4, 1, 0.2)) +
  labs(x = "Treatment",
       y = expression(bold("C"["i"]~": C"["a"]))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## iWUE
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = iwue,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, variable == "iwue"),
            aes(y = 150, 
                label = compact), 
            fontface = "bold", size = 5, hjust = 0.575) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0, 160),
                     breaks = seq(0, 160, 40)) +
  labs(x = "Treatment",
       y = expression(bold("R"["d"]~": V"["cmax25"]))) +
  guides(fill = "none") +
  pubtheme

##########################################################################
## Vcmax:gs
##########################################################################
ggplot(data = data, aes(x = factor(n.trt, levels = c("LN", "HN")),
                        y = vcmax.gs,
                        fill = n.trt)) +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 2, alpha = 0.5) +
  geom_text(data = subset(comp.letters, variable == "vcmax.gs"),
            aes(y = 1500, 
                label = compact), 
            fontface = "bold", size = 5) +
  facet_grid(.~inoc, labeller = labeller(inoc = facet.labels)) +
  scale_fill_manual(values = cbbPalette,
                    labels = c("NI" = "Not inoculated",
                               "YI" = "Inoculated")) +
  scale_x_discrete(labels = c("Low nitrogen", "High nitrogen")) +
  scale_y_continuous(limits = c(0, 1600),
                     breaks = seq(0, 1600, 400)) +
  labs(x = "Treatment",
       y = expression(bold("V"["cmax25"]~": g"["s400"]~"(units)"))) +
  guides(fill = "none") +
  pubtheme




