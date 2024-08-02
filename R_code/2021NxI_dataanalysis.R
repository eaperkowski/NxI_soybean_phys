## Libraries
library(tidyverse)
library(dplyr)
library(lme4)
library(car)
library(emmeans)
library(multcomp)

## Remove digit bounds in emmeans package
emm_options(opt.digits = FALSE)

## Load data
data <- read.csv("../data/2021NxI_trait_data.csv",
                 na.strings = "NA") %>%
  filter(inoc == "yi" | (inoc == "ni" & nodule.biomass == 0))

## Add  nodule biomass and nod.root.biomass
data$nodule.biomass <- ifelse(is.numeric(data$root.biomass) & 
                                is.na(data$nodule.biomass),
                              0, data$nodule.biomass)

## Check data
head(data)

##########################################################################
## Carbon cost to acquire nitrogen
##########################################################################
n.cost <- lmer(log(n.cost) ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(n.cost)
qqnorm(residuals(n.cost))
qqline(residuals(n.cost))
hist(residuals(n.cost))
shapiro.test(residuals(n.cost))
outlierTest(n.cost)

# Model output
summary(n.cost)
Anova(n.cost)

# Pairwise comparisons
emmeans(n.cost, pairwise~n.trt, type = "response")
emmeans(n.cost, pairwise~inoc, type = "response")
cld(emmeans(n.cost, pairwise~n.trt * inoc, type = "response"))

# Write data frame for compact lettering
ncost.pairwise.full <- data.frame(variable = "ncost",
                                  cld(emmeans(n.cost, ~n.trt*inoc, 
                                              type = "response"),
                                      Letters = letters, reversed = TRUE))

##########################################################################
## Belowground carbon
##########################################################################
bg.carbon <- lmer(log(bg.total.c) ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(bg.carbon)
qqnorm(residuals(bg.carbon))
qqline(residuals(bg.carbon))
hist(residuals(bg.carbon))
shapiro.test(residuals(bg.carbon))
outlierTest(bg.carbon)

# Model output
summary(bg.carbon)
Anova(bg.carbon)

# Pairwise comparisons
emmeans(bg.carbon, pairwise~n.trt)
emmeans(bg.carbon, pairwise~inoc, type = "response")

# Write data frame for compact lettering
bgc.pairwise.full <- data.frame(variable = "bgc",
                                cld(emmeans(bg.carbon, ~n.trt*inoc, 
                                            type = "response"),
                                    Letters = letters))

##########################################################################
## Whole plant nitrogen (denominator of carbon cost to acquire nitrogen)
##########################################################################
wp.nitrogen <- lmer(wp.total.n ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(wp.nitrogen)
qqnorm(residuals(wp.nitrogen))
qqline(residuals(wp.nitrogen))
hist(residuals(wp.nitrogen))
shapiro.test(residuals(wp.nitrogen))
outlierTest(wp.nitrogen)

# Model output
summary(wp.nitrogen)
Anova(wp.nitrogen)

# Pairwise comparisons
emmeans(wp.nitrogen, pairwise~n.trt)
emmeans(wp.nitrogen, pairwise~inoc)
emmeans(wp.nitrogen, pairwise~n.trt * inoc)

# Write data frame for compact lettering
wpn.pairwise.full <- data.frame(variable = "wpn",
                                cld(emmeans(wp.nitrogen, ~n.trt*inoc),
                                    Letters = letters)) %>%
  rename(response = emmean)

##########################################################################
## Total leaf area
##########################################################################
tla <- lmer(total.leaf.area ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(tla)
qqnorm(residuals(tla))
qqline(residuals(tla))
hist(residuals(tla))
shapiro.test(residuals(tla))
outlierTest(tla)

# Model output
summary(tla)
Anova(tla)

# Pairwise comparisons
emmeans(tla, pairwise~n.trt*inoc)
emmeans(tla, pairwise~n.trt)
emmeans(tla, pairwise~inoc)


tla.pairwise.full <- data.frame(variable = "total.leaf.area",
                                cld(emmeans(tla, ~n.trt*inoc),
                                    Letters = letters)) %>%
  rename(response = emmean)

##########################################################################
## Total biomass
##########################################################################
totalbiomass <- lmer(log(total.biomass) ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(totalbiomass)
qqnorm(residuals(totalbiomass))
qqline(residuals(totalbiomass))
hist(residuals(totalbiomass))
shapiro.test(residuals(totalbiomass))
outlierTest(totalbiomass)

# Model output
summary(totalbiomass)
Anova(totalbiomass)

# Pairwise comparisons
emmeans(totalbiomass, pairwise~n.trt*inoc)
emmeans(totalbiomass, pairwise~n.trt, type = "response")
emmeans(totalbiomass, pairwise~inoc)


tbio.pairwise.full <- data.frame(variable = "total.biomass",
                                 cld(emmeans(totalbiomass, ~n.trt*inoc),
                                     Letters = letters)) %>%
  rename(response = emmean)

##########################################################################
## Nodule biomass : root biomass
##########################################################################
nod.root <- lmer(nod.root.biomass ~ n.trt + (1 | block), 
                 data = data)

# Check model assumptions
plot(nod.root)
qqnorm(residuals(nod.root))
qqline(residuals(nod.root))
hist(residuals(nod.root))
shapiro.test(residuals(nod.root))
outlierTest(nod.root)

# Model output
summary(nod.root)
Anova(nod.root)

# Pairwise comparisons
emmeans(nod.root, pairwise~n.trt)

# Write data frame for compact lettering
nodroot.pairwise.full <- data.frame(variable = "nodroot",
                                    cld(emmeans(nod.root, ~n.trt,
                                                type = "response"),
                                        Letters = letters, reversed = TRUE)) %>%
  rename(response = emmean)

##########################################################################
## Nodule biomass
##########################################################################
nod <- lmer(nodule.biomass ~ n.trt + (1 | block), 
            data = subset(data, inoc == "yi"))

# Check model assumptions
plot(nod)
qqnorm(residuals(nod))
qqline(residuals(nod))
hist(residuals(nod))
shapiro.test(residuals(nod))
outlierTest(nod)

# Model output
summary(nod)
Anova(nod)

# Pairwise comparisons
emmeans(nod, pairwise~n.trt)

# Write data frame for compact lettering
nod.pairwise.full <- data.frame(variable = "nod",
                                cld(emmeans(nod, ~n.trt),
                                    Letters = letters, reversed = TRUE)) %>%
  rename(response = emmean)

##########################################################################
## Root biomass
##########################################################################
root <- lmer(log(root.biomass) ~ n.trt * inoc + (1 | block), 
             data = data)

# Check model assumptions
plot(root)
qqnorm(residuals(root))
qqline(residuals(root))
hist(residuals(root))
shapiro.test(residuals(root))
outlierTest(root)

# Model output
summary(root)
Anova(root)

# Pairwise comparisons
emmeans(root, pairwise~inoc, type = "response")

# Write data frame for compact lettering
root.pairwise.full <- data.frame(variable = "root",
                                 cld(emmeans(root, ~n.trt*inoc,
                                             type = "response"),
                                     Letters = letters))

##########################################################################
## BVR
##########################################################################
bvr <- lmer(log(bvr) ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(bvr)
qqnorm(residuals(bvr))
qqline(residuals(bvr))
hist(residuals(bvr))
shapiro.test(residuals(bvr))
outlierTest(bvr)

# Model output
summary(bvr)
Anova(bvr)

# Pairwise comparisons
emmeans(bvr, pairwise~n.trt*inoc, type = "response")

# Write data frame for compact lettering
bvr.pairwise.full <- data.frame(variable = "bvr",
                                cld(emmeans(bvr, ~n.trt*inoc,
                                            type = "response"),
                                    Letters = letters))

##########################################################################
## Make merged emmeans file
##########################################################################
comp.letters <- ncost.pairwise.full %>%
  full_join(bgc.pairwise.full) %>%
  full_join(wpn.pairwise.full) %>%
  full_join(tla.pairwise.full) %>%
  full_join(tbio.pairwise.full) %>%
  full_join(nodroot.pairwise.full) %>%
  full_join(nod.pairwise.full) %>%
  full_join(root.pairwise.full) %>%
  full_join(bvr.pairwise.full) %>%
  unite("treatment", n.trt:inoc, remove = FALSE) %>%
  mutate(treatment = factor(treatment, levels = c("70_ni", "630_ni",
                                                  "70_yi", "630_yi")),
         
         .group = trimws(.group, "both"))

comp.letters

## Write pairwise comparison csv file
write.csv(comp.letters, "../data/2021NxI_compact_letters.csv", 
          row.names = FALSE)
