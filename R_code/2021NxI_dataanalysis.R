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
                 na.strings = "NA")

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
emmeans(n.cost, pairwise~n.trt * inoc, type = "response")

# Write data frame for compact lettering
ncost.pairwise.full <- data.frame(variable = "ncost",
                                  treatment = "full",
                                  cld(emmeans(n.cost, ~n.trt*inoc, 
                                              type = "response"),
                                      Letters = letters))
ncost.pairwise.soiln <- data.frame(variable = "ncost",
                                   treatment = "n.trt",
                                   cld(emmeans(n.cost, ~n.trt,
                                               type = "response"),
                                       Letters = letters))
ncost.pairwise.inoc <- data.frame(variable = "ncost",
                                  treatment = "inoc",
                                  cld(emmeans(n.cost, ~inoc,
                                              type = "response"),
                                      Letters = letters))
ncost.pairwise <- ncost.pairwise.full %>%
  full_join(ncost.pairwise.soiln) %>%
  full_join(ncost.pairwise.inoc) %>%
  dplyr::rename(emmean = response) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

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
                                treatment = "full",
                                cld(emmeans(bg.carbon, ~n.trt*inoc, 
                                            type = "response"),
                                    Letters = letters))
bgc.pairwise.soiln <- data.frame(variable = "bgc",
                                 treatment = "n.trt",
                                 cld(emmeans(bg.carbon, ~n.trt,
                                             type = "response"),
                                     Letters = letters))
bgc.pairwise.inoc <- data.frame(variable = "bgc",
                                treatment = "inoc",
                                cld(emmeans(bg.carbon, ~inoc,
                                            type = "response"),
                                    Letters = letters))
bgc.pairwise <- bgc.pairwise.full %>%
  full_join(bgc.pairwise.soiln) %>%
  full_join(bgc.pairwise.inoc) %>%
  dplyr::rename(emmean = response) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

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
                                treatment = "full",
                                cld(emmeans(wp.nitrogen, ~n.trt*inoc),
                                    Letters = letters))
wpn.pairwise.soiln <- data.frame(variable = "wpn",
                                 treatment = "n.trt",
                                 cld(emmeans(wp.nitrogen, ~n.trt),
                                     Letters = letters))
wpn.pairwise.inoc <- data.frame(variable = "wpn",
                                treatment = "inoc",
                                cld(emmeans(wp.nitrogen, ~inoc),
                                    Letters = letters))
wpn.pairwise <- wpn.pairwise.full %>%
  full_join(wpn.pairwise.soiln) %>%
  full_join(wpn.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)


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
r.squaredGLMM(tla)

# Pairwise comparisons
emmeans(tla, pairwise~n.trt*inoc)
emmeans(tla, pairwise~n.trt)
emmeans(tla, pairwise~inoc)


tla.pairwise.full <- data.frame(variable = "total.leaf.area",
                               treatment = "full",
                               cld(emmeans(tla, ~n.trt*inoc),
                                   Letters = letters))
tla.pairwise.soiln <- data.frame(variable = "total.leaf.area",
                                treatment = "n.trt",
                                cld(emmeans(tla, ~n.trt),
                                    Letters = letters))
tla.pairwise.inoc <- data.frame(variable = "total.leaf.area",
                               treatment = "inoc",
                               cld(emmeans(tla, ~inoc),
                                   Letters = letters))
tla.pairwise <- tla.pairwise.full %>%
  full_join(tla.pairwise.soiln) %>%
  full_join(tla.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

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
r.squaredGLMM(totalbiomass)

# Pairwise comparisons
emmeans(totalbiomass, pairwise~n.trt*inoc)
emmeans(totalbiomass, pairwise~n.trt, type = "response")
emmeans(totalbiomass, pairwise~inoc)


tbio.pairwise.full <- data.frame(variable = "total.biomass",
                                treatment = "full",
                                cld(emmeans(totalbiomass, ~n.trt*inoc),
                                    Letters = letters))
tbio.pairwise.soiln <- data.frame(variable = "total.biomass",
                                 treatment = "n.trt",
                                 cld(emmeans(totalbiomass, ~n.trt),
                                     Letters = letters))
tbio.pairwise.inoc <- data.frame(variable = "total.biomass",
                                treatment = "inoc",
                                cld(emmeans(totalbiomass, ~inoc),
                                    Letters = letters))
tbio.pairwise <- tbio.pairwise.full %>%
  full_join(tbio.pairwise.soiln) %>%
  full_join(tbio.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Nodule biomass : root biomass
##########################################################################
nod.root <- lmer(sqrt(nod.root.biomass) ~ n.trt * inoc + (1 | block), data = data)

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
emmeans(nod.root, pairwise~inoc)
cld(emmeans(nod.root, pairwise~n.trt * inoc, type = "response"))

# Write data frame for compact lettering
nodroot.pairwise.full <- data.frame(variable = "nodroot",
                                treatment = "full",
                                cld(emmeans(nod.root, ~n.trt*inoc,
                                            type = "response"),
                                    Letters = letters))
nodroot.pairwise.soiln <- data.frame(variable = "nodroot",
                                 treatment = "n.trt",
                                 cld(emmeans(nod.root, ~n.trt,
                                             type = "response"),
                                     Letters = letters))
nodroot.pairwise.inoc <- data.frame(variable = "nodroot",
                                treatment = "inoc",
                                cld(emmeans(nod.root, ~inoc,
                                            type = "response"),
                                    Letters = letters))
nodroot.pairwise <- nodroot.pairwise.full %>%
  full_join(nodroot.pairwise.soiln) %>%
  full_join(nodroot.pairwise.inoc) %>%
  dplyr::rename(emmean = response) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Nodule biomass
##########################################################################
nod <- lmer(sqrt(nodule.biomass) ~ n.trt * inoc + (1 | block), data = data)

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
emmeans(nod, pairwise~inoc, type = "response")

# Write data frame for compact lettering
nod.pairwise.full <- data.frame(variable = "nod",
                                    treatment = "full",
                                    cld(emmeans(nod, ~n.trt*inoc,
                                                type = "response"),
                                        Letters = letters))
nod.pairwise.soiln <- data.frame(variable = "nod",
                                     treatment = "n.trt",
                                     cld(emmeans(nod, ~n.trt,
                                                 type = "response"),
                                         Letters = letters))
nod.pairwise.inoc <- data.frame(variable = "nod",
                                    treatment = "inoc",
                                    cld(emmeans(nod, ~inoc,
                                                type = "response"),
                                        Letters = letters))
nod.pairwise <- nod.pairwise.full %>%
  full_join(nod.pairwise.soiln) %>%
  full_join(nod.pairwise.inoc) %>%
  dplyr::rename(emmean = response) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Nodule biomass
##########################################################################
root <- lmer(log(root.biomass) ~ n.trt * inoc + (1 | block), data = data)

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
emmeans(root, pairwise~n.trt)
emmeans(root, pairwise~inoc, type = "response")

# Write data frame for compact lettering
root.pairwise.full <- data.frame(variable = "root",
                                treatment = "full",
                                cld(emmeans(root, ~n.trt*inoc,
                                            type = "response"),
                                    Letters = letters))
root.pairwise.soiln <- data.frame(variable = "root",
                                 treatment = "n.trt",
                                 cld(emmeans(root, ~n.trt,
                                             type = "response"),
                                     Letters = letters))
root.pairwise.inoc <- data.frame(variable = "root",
                                treatment = "inoc",
                                cld(emmeans(root, ~inoc,
                                            type = "response"),
                                    Letters = letters))
root.pairwise <- root.pairwise.full %>%
  full_join(root.pairwise.soiln) %>%
  full_join(root.pairwise.inoc) %>%
  dplyr::rename(emmean = response) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

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
                                 treatment = "full",
                                 cld(emmeans(bvr, ~n.trt*inoc,
                                             type = "response"),
                                     Letters = letters))
bvr.pairwise.soiln <- data.frame(variable = "bvr",
                                  treatment = "n.trt",
                                  cld(emmeans(bvr, ~n.trt,
                                              type = "response"),
                                      Letters = letters))
bvr.pairwise.inoc <- data.frame(variable = "bvr",
                                 treatment = "inoc",
                                 cld(emmeans(bvr, ~inoc,
                                             type = "response"),
                                     Letters = letters))
bvr.pairwise <- bvr.pairwise.full %>%
  full_join(bvr.pairwise.soiln) %>%
  full_join(bvr.pairwise.inoc) %>%
  dplyr::rename(emmean = response) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Make merged emmeans file
##########################################################################
comp.letters <- ncost.pairwise %>%
  full_join(bgc.pairwise) %>%
  full_join(wpn.pairwise) %>%
  full_join(tla.pairwise) %>%
  full_join(tbio.pairwise) %>%
  full_join(nodroot.pairwise) %>%
  full_join(nod.pairwise) %>%
  full_join(root.pairwise)
  full_join(bvr.pairwise) %>%
  dplyr::rename(comparison = treatment) %>%
  unite("treatment", n.trt:inoc, remove = "FALSE") %>%
  mutate(treatment = factor(treatment, levels = c("ln_ni", "hn_ni",
                                                  "ln_yi", "hn_yi")),
         compact = tolower(compact))

comp.letters

## Write pairwise comparison csv file
write.csv(comp.letters, "../data/2021NxI_compact_letters.csv", 
          row.names = FALSE)
