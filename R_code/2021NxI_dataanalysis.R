## Libraries
library(lme4)
library(emmeans)
library(car)
library(tidyverse)
library(dplyr)
library(MuMIn)
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

data$nod.root.biomass <- data$nodule.biomass / data$root.biomass

## Check data
head(data)


##########################################################################
## Nmass
##########################################################################
nmass <- lmer(leaf.n ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(nmass)
qqnorm(residuals(nmass))
qqline(residuals(nmass))
hist(residuals(nmass))
shapiro.test(residuals(nmass))
outlierTest(nmass)

# Model output
summary(nmass)
Anova(nmass)
r.squaredGLMM(nmass)

# Pairwise comparisons
emmeans(nmass, pairwise ~ n.trt)
emmeans(nmass, pairwise ~ n.trt * inoc)

# Emmean for fig making
nmass.pairwise.full <- data.frame(variable = "nmass",
                                treatment = "full",
                                cld(emmeans(nmass, ~n.trt*inoc),
                                    Letters = letters))
nmass.pairwise.soiln <- data.frame(variable = "nmass",
                                 treatment = "n.trt",
                                 cld(emmeans(nmass,
                                             ~n.trt),
                                     Letters = letters))
nmass.pairwise.inoc <- data.frame(variable = "nmass",
                                treatment = "inoc",
                                cld(emmeans(nmass,
                                            ~inoc),
                                    Letters = letters))
nmass.pairwise <- nmass.pairwise.full %>%
  full_join(nmass.pairwise.soiln) %>%
  full_join(nmass.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Specific leaf area
##########################################################################
sla <- lmer(sla ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(sla)
qqnorm(residuals(sla))
qqline(residuals(sla))
hist(residuals(sla))
shapiro.test(residuals(sla))
outlierTest(sla)

# Model output
summary(sla)
Anova(sla)
r.squaredGLMM(sla)

# Pairwise comparisons
emmeans(sla, pairwise~inoc)
emmeans(sla, pairwise~n.trt)

# Emmean for fig making
sla.pairwise.full <- data.frame(variable = "sla",
                                treatment = "full",
                                cld(emmeans(sla, ~n.trt*inoc),
                                    Letters = letters))
sla.pairwise.soiln <- data.frame(variable = "sla",
                                 treatment = "n.trt",
                                 cld(emmeans(sla,
                                             ~n.trt),
                                     Letters = letters))
sla.pairwise.inoc <- data.frame(variable = "sla",
                                treatment = "inoc",
                                cld(emmeans(sla,
                                            ~inoc),
                                    Letters = letters))
sla.pairwise <- sla.pairwise.full %>%
  full_join(sla.pairwise.soiln) %>%
  full_join(sla.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Narea
##########################################################################
data$narea[11] <- NA

narea <- lmer(narea ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(narea)
qqnorm(residuals(narea))
qqline(residuals(narea))
hist(residuals(narea))
shapiro.test(residuals(narea))
outlierTest(narea)

# Model output
summary(narea)
Anova(narea)
r.squaredGLMM(narea)

# Pairwise comparisons
emmeans(narea, pairwise ~ n.trt * inoc)
emmeans(narea, pairwise~n.trt)

# Emmean for fig making
narea.pairwise.full <- data.frame(variable = "narea",
                                  treatment = "full",
                                  cld(emmeans(narea, ~n.trt*inoc),
                                      Letters = letters))
narea.pairwise.soiln <- data.frame(variable = "narea",
                                   treatment = "n.trt",
                                   cld(emmeans(narea,
                                               ~n.trt),
                                       Letters = letters))
narea.pairwise.inoc <- data.frame(variable = "narea",
                                  treatment = "inoc",
                                  cld(emmeans(narea,
                                              ~inoc),
                                      Letters = letters))
narea.pairwise <- narea.pairwise.full %>%
  full_join(narea.pairwise.soiln) %>%
  full_join(narea.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## A400
##########################################################################
a400 <- lmer(anet ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(a400)
qqnorm(residuals(a400))
qqline(residuals(a400))
hist(residuals(a400))
shapiro.test(residuals(a400))
outlierTest(a400)

# Model output
summary(a400)
Anova(a400)
r.squaredGLMM(a400)

# Pairwise comparisons
emmeans(a400, pairwise~n.trt)

# Emmean output for fig making
a400.pairwise.full <- data.frame(variable = "a400",
                                 treatment = "full",
                                 cld(emmeans(a400, ~n.trt*inoc),
                                Letters = letters))
a400.pairwise.soiln <- data.frame(variable = "a400",
                                  treatment = "n.trt",
                                  cld(emmeans(a400,
                                              ~n.trt),
                                      Letters = letters))
a400.pairwise.inoc <- data.frame(variable = "a400",
                                  treatment = "inoc",
                                  cld(emmeans(a400,
                                              ~inoc),
                                      Letters = letters))

a400.pairwise <- a400.pairwise.full %>%
  full_join(a400.pairwise.soiln) %>%
  full_join(a400.pairwise.inoc) %>%
  mutate(.group = trimws(.group, which = "both"),
         compact = c("b", "ab", "a", "a", "b", "a", "a", "a")) %>%
  data.frame()

##########################################################################
## Vcmax25 with TPU
##########################################################################
vcmax25 <- lmer(vcmax25 ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(vcmax25)
qqnorm(residuals(vcmax25))
qqline(residuals(vcmax25))
hist(residuals(vcmax25))
shapiro.test(residuals(vcmax25))
outlierTest(vcmax25)

# Model output
summary(vcmax25)
Anova(vcmax25)
r.squaredGLMM(vcmax25)

# Pairwise comparisons
emmeans(vcmax25, pairwise~n.trt) # Nitrogen addition decreases Vcmax25??

# Emmean output for fig making
vcmax.pairwise.full <- data.frame(variable = "vcmax25",
                                  treatment = "full",
                                  cld(emmeans(vcmax25, ~n.trt*inoc),
                                      Letters = letters))
vcmax.pairwise.soiln <- data.frame(variable = "vcmax25",
                                   treatment = "n.trt",
                                   cld(emmeans(vcmax25, ~n.trt),
                                       Letters = letters))
vcmax.pairwise.inoc <- data.frame(variable = "vcmax25",
                                  treatment = "inoc",
                                  cld(emmeans(vcmax25, ~inoc),
                                      Letters = letters))
vcmax.pairwise <- vcmax.pairwise.full %>%
  full_join(vcmax.pairwise.soiln) %>%
  full_join(vcmax.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = c("a", "a", "a", "a", "b", "a", "a", "a"))

##########################################################################
## Jmax25 with TPU
##########################################################################
jmax25 <- lmer(jmax25 ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(jmax25)
qqnorm(residuals(jmax25))
qqline(residuals(jmax25))
hist(residuals(jmax25))
shapiro.test(residuals(jmax25))
outlierTest(jmax25)

# Model output
summary(jmax25)
Anova(jmax25)
r.squaredGLMM(jmax25)

# Pairwise comparisons
emmeans(jmax25, pairwise~n.trt)
emmeans(jmax25, pairwise~inoc)

# Emmean output for fig making
jmax.pairwise.full <- data.frame(variable = "jmax25",
                                 treatment = "full",
                                 cld(emmeans(jmax25, ~n.trt*inoc),
                                 Letters = letters))
jmax.pairwise.soiln <- data.frame(variable = "jmax25",
                                  treatment = "n.trt",
                                  cld(emmeans(jmax25, ~n.trt),
                                      Letters = letters))
jmax.pairwise.inoc <- data.frame(variable = "jmax25",
                                 treatment = "inoc",
                                 cld(emmeans(jmax25, ~inoc),
                                     Letters = letters))
jmax.pairwise <- jmax.pairwise.full %>%
  full_join(jmax.pairwise.soiln) %>%
  full_join(jmax.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = c("a", "a", "a", "a", "b", "a", "a", "a"))

##########################################################################
## Jmax25:Vcmax25
##########################################################################
data$jmax25.vcmax25[c(17, 62)] <- NA

jmax25.vcmax25 <- lmer(log(jmax25.vcmax25) ~ n.trt * inoc + (1 | block),
                       data = data)


# Check model assumptions
plot(jmax25.vcmax25)
qqnorm(residuals(jmax25.vcmax25))
qqline(residuals(jmax25.vcmax25))
hist(residuals(jmax25.vcmax25))
shapiro.test(residuals(jmax25.vcmax25))
outlierTest(jmax25.vcmax25)

# Model output
summary(jmax25.vcmax25)
Anova(jmax25.vcmax25)
r.squaredGLMM(jmax25.vcmax25)

# Pairwise comparisons
emmeans(jmax25.vcmax25, pairwise~n.trt, type = "response")
emmeans(jmax25.vcmax25, pairwise~inoc, type = "response")

cld(emmeans(jmax25.vcmax25, pairwise~n.trt*inoc, type = "response")) 
# High nitrogen decreases Jmax25:Vcmax25 when plants are not inoculated
# Nitrogen status has no impact on Jmax25:Vcmax25 when plants are inoculated

# Emmean output for fig making
jmax.vcmax.pairwise.full <- data.frame(variable = "jmax25.vcmax25",
                                       treatment = "full",
                                       cld(emmeans(jmax25.vcmax25, ~n.trt*inoc,
                                                   type = "response"),
                                           Letters = letters))
jmax.vcmax.pairwise.soiln <- data.frame(variable = "jmax25.vcmax25",
                                        treatment = "n.trt",
                                        cld(emmeans(jmax25.vcmax25, ~n.trt,
                                                    type = "response"),
                                            Letters = letters))
jmax.vcmax.pairwise.inoc <- data.frame(variable = "jmax25.vcmax25",
                                        treatment = "inoc",
                                        cld(emmeans(jmax25.vcmax25, ~inoc,
                                                    type = "response"),
                                            Letters = letters))
jmax.vcmax.pairwise <- jmax.vcmax.pairwise.full %>%
  full_join(jmax.vcmax.pairwise.soiln) %>%
  full_join(jmax.vcmax.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = c("a", "a", "a", "a", "a", "a", "a", "b"))

##########################################################################
## Rd (standardized to 25 deg C)
##########################################################################
data$rd25[c(7, 31)] <- NA

rd <- lmer(log(rd25) ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(rd)
qqnorm(residuals(rd))
qqline(residuals(rd))
hist(residuals(rd))
shapiro.test(residuals(rd))
outlierTest(rd)

# Model output
summary(rd)
Anova(rd)
r.squaredGLMM(rd)

# Pairwise comparisons
emmeans(rd, pairwise~n.trt, type = "response") 
emmeans(rd, pairwise~n.trt*inoc, type = "response")

# Emmean for fig making
rd.pairwise.full <- data.frame(variable = "rd25",
                               treatment = "full",
                               cld(emmeans(rd, ~n.trt*inoc, 
                                           type = "response"),
                                   Letters = letters))
rd.pairwise.soiln <- data.frame(variable = "rd25",
                                treatment = "n.trt",
                                cld(emmeans(rd, ~n.trt, 
                                            type = "response"),
                                    Letters = letters))
rd.pairwise.inoc <- data.frame(variable = "rd25",
                               treatment = "inoc",
                               cld(emmeans(rd, ~inoc, 
                                           type = "response"),
                                    Letters = letters))
rd.pairwise <- rd.pairwise.full %>%
  full_join(rd.pairwise.soiln) %>%
  full_join(rd.pairwise.inoc) %>%
  dplyr::rename(emmean = response) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Rd:Vcmax (Vcmax is not standardized since Rd is not temp standardized)
##########################################################################
data$rd25.vcmax25[c(7, 31)] <- NA

rd.vcmax <- lmer(rd25.vcmax25 ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(rd.vcmax)
qqnorm(residuals(rd.vcmax))
qqline(residuals(rd.vcmax))
hist(residuals(rd.vcmax))
shapiro.test(residuals(rd.vcmax))
outlierTest(rd.vcmax)

# Model output
summary(rd.vcmax)
Anova(rd.vcmax)
r.squaredGLMM(rd.vcmax)

# Pairwise comparisons
emmeans(rd.vcmax, pairwise~n.trt)
# High nitrogen increases Rd:Vcmax regardless of inoculation status
# Driven by both an increase in Rd and reduction in Vcmax with increasing soil N
emmeans(rd.vcmax, pairwise~n.trt*inoc)
cld(emmeans(rd.vcmax, pairwise~n.trt*inoc))

# Emmean for fig making
rd.vcmax.pairwise.full <- data.frame(variable = "rd.vcmax",
                                     treatment = "full",
                                     cld(emmeans(rd.vcmax, ~n.trt*inoc),
                                         Letters = letters))
rd.vcmax.pairwise.soiln <- data.frame(variable = "rd25.vcmax25",
                                      treatment = "n.trt",
                                      cld(emmeans(rd.vcmax, ~n.trt),
                                          Letters = letters))
rd.vcmax.pairwise.inoc <- data.frame(variable = "rd25.vcmax25",
                                     treatment = "inoc",
                                     cld(emmeans(rd.vcmax, ~inoc),
                                         Letters = letters))
rd.vcmax.pairwise <- rd.vcmax.pairwise.full %>%
  full_join(rd.vcmax.pairwise.soiln) %>%
  full_join(rd.vcmax.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Gs
##########################################################################
data$gsw[64] <- NA

gs <- lmer(gsw ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(gs)
qqnorm(residuals(gs))
qqline(residuals(gs))
hist(residuals(gs))
shapiro.test(residuals(gs))
outlierTest(gs)

# Model output
summary(gs)
Anova(gs)
r.squaredGLMM(gs)

# Pairwise comparisons
emmeans(gs, pairwise~n.trt)
# Increasing nitrogen decreased stomatal conductance

# Emmean for fig making
gs.pairwise.full <- data.frame(variable = "gs",
                               treatment = "full",
                               cld(emmeans(gs, ~n.trt*inoc),
                                   Letters = letters))
gs.pairwise.soiln <- data.frame(variable = "gs",
                                treatment = "n.trt",
                                cld(emmeans(gs, ~n.trt),
                                    Letters = letters))
gs.pairwise.inoc <- data.frame(variable = "gs",
                               treatment = "inoc",
                               cld(emmeans(gs, ~inoc),
                                   Letters = letters))
gs.pairwise <- gs.pairwise.full %>%
  full_join(gs.pairwise.soiln) %>%
  full_join(gs.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = c("b", "b", "a", "a", "b", "a", "a", "a"))

##########################################################################
## Ci:Ca
##########################################################################
cica <- lmer(ci.ca ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(cica)
qqnorm(residuals(cica))
qqline(residuals(cica))
hist(residuals(cica))
shapiro.test(residuals(cica))
outlierTest(cica)

# Model output
summary(cica)
Anova(cica)
r.squaredGLMM(cica)

# Pairwise comparisons
emmeans(cica, pairwise~n.trt)
# Increasing nitrogen decreases Ci:Ca

# Emmean for fig making
cica.pairwise.full <- data.frame(variable = "ci.ca",
                                 treatment = "full",
                            cld(emmeans(cica, ~n.trt*inoc),
                                Letters = letters))
cica.pairwise.soiln <- data.frame(variable = "ci.ca",
                                  treatment = "n.trt",
                                  cld(emmeans(cica, ~n.trt),
                                      Letters = letters))
cica.pairwise.inoc <- data.frame(variable = "ci.ca",
                                 treatment = "inoc",
                                 cld(emmeans(cica, ~inoc),
                                      Letters = letters))
cica.pairwise <- cica.pairwise.full %>%
  full_join(cica.pairwise.soiln) %>%
  full_join(cica.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = c("a", "a", "a", "a", "b", "a", "a", "a"))

##########################################################################
## PNUE
##########################################################################
data$pnue[11] <- NA

pnue <- lmer(pnue ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(pnue)
qqnorm(residuals(pnue))
qqline(residuals(pnue))
hist(residuals(pnue))
shapiro.test(residuals(pnue))
outlierTest(pnue)

# Model output
summary(pnue)
Anova(pnue)
r.squaredGLMM(pnue)

# Pairwise comparisons
emmeans(pnue, pairwise~n.trt*inoc)
emmeans(pnue, pairwise~n.trt)

# Emmean for fig making
pnue.pairwise.full <- data.frame(variable = "pnue",
                                 treatment = "full",
                                 cld(emmeans(pnue, ~n.trt*inoc,
                                             type = "response"),
                                     Letters = letters))
pnue.pairwise.soiln <- data.frame(variable = "pnue",
                                  treatment = "n.trt",
                                  cld(emmeans(pnue, ~n.trt, 
                                              type = "response"),
                                      Letters = letters))
pnue.pairwise.inoc <- data.frame(variable = "pnue",
                                 treatment = "inoc",
                                 cld(emmeans(pnue, ~inoc, 
                                             type = "response"),
                                     Letters = letters))
pnue.pairwise <- pnue.pairwise.full %>%
  full_join(pnue.pairwise.soiln) %>%
  full_join(pnue.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = c("b", "b", "a", "a", "b", "a", "a", "a"))

##########################################################################
## iWUE
##########################################################################
iwue <- lmer(iwue ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(iwue)
qqnorm(residuals(iwue))
qqline(residuals(iwue))
hist(residuals(iwue))
shapiro.test(residuals(iwue))
outlierTest(iwue)

# Model output
summary(iwue)
Anova(iwue)
r.squaredGLMM(iwue)

# Pairwise comparisons
emmeans(iwue, pairwise~n.trt, type = "response")
# Increasing nitrogen increases iWUE

# Emmean for fig making
iwue.pairwise.full <- data.frame(variable = "iwue",
                                 treatment = "full",
                                 cld(emmeans(iwue, ~n.trt*inoc,
                                             type = "response"),
                                     Letters = letters))
iwue.pairwise.soiln <- data.frame(variable = "iwue",
                                  treatment = "n.trt",
                                  cld(emmeans(iwue, ~n.trt, 
                                              type = "response"),
                                      Letters = letters))
iwue.pairwise.inoc <- data.frame(variable = "iwue",
                                 treatment = "inoc",
                                  cld(emmeans(iwue, ~inoc, 
                                              type = "response"),
                                      Letters = letters))
iwue.pairwise <- iwue.pairwise.full %>%
  full_join(iwue.pairwise.soiln) %>%
  full_join(iwue.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Vcmax:gs
##########################################################################
vcmax.gs <- lmer(log(vcmax.gs) ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(vcmax.gs)
qqnorm(residuals(vcmax.gs))
qqline(residuals(vcmax.gs))
hist(residuals(vcmax.gs))
shapiro.test(residuals(vcmax.gs))
outlierTest(vcmax.gs)

# Model output
summary(vcmax.gs)
Anova(vcmax.gs)
r.squaredGLMM(vcmax.gs)

# Pairwise comparisons
emmeans(vcmax.gs, pairwise~n.trt)

# Emmean for fig making
vcmax.gs.pairwise.full <- data.frame(variable = "vcmax.gs",
                                     treatment = "full",
                                     cld(emmeans(vcmax.gs, ~n.trt*inoc, 
                                                 type = "response"),
                                         Letters = letters))
vcmax.gs.pairwise.soiln <- data.frame(variable = "vcmax.gs",
                                      treatment = "n.trt",
                                      cld(emmeans(vcmax.gs, ~n.trt, 
                                                  type = "response"),
                                          Letters = letters))
vcmax.gs.pairwise.inoc <- data.frame(variable = "vcmax.gs",
                                      treatment = "inoc",
                                      cld(emmeans(vcmax.gs, ~inoc, 
                                                  type = "response"),
                                          Letters = letters))
vcmax.gs.pairwise <- vcmax.gs.pairwise.full %>%
  full_join(vcmax.gs.pairwise.soiln) %>%
  full_join(vcmax.gs.pairwise.inoc) %>%
  dplyr::rename(emmean = response) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Narea:gs
##########################################################################
data$narea.gs[11] <- NA

narea.gs <- lmer(log(narea.gs) ~ n.trt * inoc + (1 | block), data = data)

# Check model assumptions
plot(narea.gs)
qqnorm(residuals(narea.gs))
qqline(residuals(narea.gs))
hist(residuals(narea.gs))
shapiro.test(residuals(narea.gs))
outlierTest(narea.gs)

# Model output
summary(narea.gs)
Anova(narea.gs)
r.squaredGLMM(narea.gs)

# Pairwise comparisons
emmeans(narea.gs, pairwise~n.trt)
# Increasing nitrogen increases Narea:gs through and increase in Narea and
# decrease in gs

# Emmean for fig making
narea.gs.pairwise.full <- data.frame(variable = "narea.gs",
                                     treatment = "full",
                                     cld(emmeans(narea.gs, ~n.trt*inoc, 
                                                 type = "response"),
                                         Letters = letters))
narea.gs.pairwise.soiln <- data.frame(variable = "narea.gs",
                                      treatment = "n.trt",
                                      cld(emmeans(narea.gs, ~n.trt, 
                                                  type = "response"),
                                          Letters = letters))
narea.gs.pairwise.inoc <- data.frame(variable = "narea.gs",
                                     treatment = "inoc",
                                     cld(emmeans(narea.gs, ~inoc, 
                                                 type = "response"),
                                         Letters = letters))
narea.gs.pairwise <- narea.gs.pairwise.full %>%
  full_join(narea.gs.pairwise.soiln) %>%
  full_join(narea.gs.pairwise.inoc) %>%
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
## Carbon cost to acquire nitrogen
##########################################################################
## Set up model
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
## Set up model
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
## Set up model
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
## Nodule biomass : root biomass
##########################################################################
## Set up model
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
## Set up model
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
## Root biomass
##########################################################################
## Set up model
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
## Set up model
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
emmeans(bvr, pairwise~n.trt)

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
comp.letters <- nmass.pairwise %>%
  full_join(sla.pairwise) %>%
  full_join(narea.pairwise) %>%
  full_join(tla.pairwise) %>%
  full_join(a400.pairwise) %>%
  full_join(vcmax.pairwise) %>%
  full_join(jmax.pairwise) %>%
  full_join(jmax.vcmax.pairwise) %>%
  full_join(rd.pairwise) %>%
  full_join(rd.vcmax.pairwise) %>%
  full_join(gs.pairwise) %>%
  full_join(cica.pairwise) %>%
  full_join(pnue.pairwise) %>%
  full_join(iwue.pairwise) %>%
  full_join(vcmax.gs.pairwise) %>%
  full_join(narea.gs.pairwise) %>%
  full_join(tbio.pairwise) %>%
  full_join(ncost.pairwise) %>%
  full_join(bgc.pairwise) %>%
  full_join(wpn.pairwise) %>%
  full_join(nodroot.pairwise) %>%
  full_join(nod.pairwise) %>%
  full_join(root.pairwise) %>%
  full_join(bvr.pairwise) %>%
  dplyr::rename(comparison = treatment) %>%
  unite("treatment", n.trt:inoc, remove = "FALSE") %>%
  mutate(treatment = factor(treatment, levels = c("ln_ni", "hn_ni",
                                                  "ln_yi", "hn_yi")),
         compact = tolower(compact)) %>%
  data.frame()
comp.letters

## Write pairwise comparison csv file
write.csv(comp.letters, "../data/2021NxI_compact_letters.csv", 
          row.names = FALSE)
