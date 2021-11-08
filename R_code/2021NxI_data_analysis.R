## Libraries
library(lme4)
library(emmeans)
library(car)
library(tidyverse)
library(dplyr)
library(MuMIn)
library(multcomp)
library(rcompanion)

## Remove digit bounds in emmeans package
emm_options(opt.digits = FALSE)

## Load data
data <- read.csv("../data/2021NxI_trait_data.csv",
                 na.strings = "NA")

## Check data
head(data)

##########################################################################
## A400
##########################################################################
a400 <- lmer(a ~ n.trt * inoc + (1 | block),
             data = data)

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
                                Letters = LETTERS))
a400.pairwise.soiln <- data.frame(variable = "a400",
                                  treatment = "n.trt",
                                  cld(emmeans(a400,
                                              ~n.trt),
                                      Letters = LETTERS))
a400.pairwise.inoc <- data.frame(variable = "a400",
                                  treatment = "inoc",
                                  cld(emmeans(a400,
                                              ~inoc),
                                      Letters = LETTERS))

a400.pairwise <- a400.pairwise.full %>%
  full_join(a400.pairwise.soiln) %>%
  full_join(a400.pairwise.inoc) %>%
  mutate(.group = trimws(.group, which = "both"),
         compact = c("B", "AB", "A", "A", "B", "A", "A", "A")) %>%
  data.frame()

##########################################################################
## Vcmax25 with TPU
##########################################################################
vcmax25 <- lmer(vcmax25 ~ n.trt * inoc + (1 | block),
                data = data)

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
                                      Letters = LETTERS))
vcmax.pairwise.soiln <- data.frame(variable = "vcmax25",
                                   treatment = "n.trt",
                                   cld(emmeans(vcmax25, ~n.trt),
                                       Letters = LETTERS))
vcmax.pairwise.inoc <- data.frame(variable = "vcmax25",
                                  treatment = "inoc",
                                  cld(emmeans(vcmax25, ~inoc),
                                      Letters = LETTERS))
vcmax.pairwise <- vcmax.pairwise.full %>%
  full_join(vcmax.pairwise.soiln) %>%
  full_join(vcmax.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = c("A", "A", "A", "A", "B", "A", "A", "A"))

##########################################################################
## Jmax25 with TPU
##########################################################################
jmax25 <- lmer(jmax25 ~ n.trt * inoc + (1 | block),
               data = data)

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
emmeans(jmax25, pairwise~n.trt) # Nitrogen addition decreases Jmax25?

# Emmean output for fig making
jmax.pairwise.full <- data.frame(variable = "jmax25",
                                 treatment = "full",
                                 cld(emmeans(jmax25, ~n.trt*inoc),
                                 Letters = LETTERS))
jmax.pairwise.soiln <- data.frame(variable = "jmax25",
                                  treatment = "n.trt",
                                  cld(emmeans(jmax25, ~n.trt),
                                      Letters = LETTERS))
jmax.pairwise.inoc <- data.frame(variable = "jmax25",
                                 treatment = "inoc",
                                 cld(emmeans(jmax25, ~inoc),
                                     Letters = LETTERS))
jmax.pairwise <- jmax.pairwise.full %>%
  full_join(jmax.pairwise.soiln) %>%
  full_join(jmax.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = c("B", "AB", "AB", "A", "B", "A", "A", "A"))

##########################################################################
## Jmax:Vcmax25
##########################################################################
data$jmax25.vcmax25[c(46, 49)] <- NA

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
cld(emmeans(jmax25.vcmax25, pairwise~n.trt*inoc)) 
# High nitrogen decreases Jmax25:Vcmax25 when plants are not inoculated
# Nitrogen status has no impact on Jmax25:Vcmax25 when plants are inoculated

# Emmean output for fig making
jmax.vcmax.pairwise.full <- data.frame(variable = "jmax25.vcmax25",
                                       treatment = "full",
                                       cld(emmeans(jmax25.vcmax25, ~n.trt*inoc,
                                                   type = "response"),
                                           Letters = LETTERS))
jmax.vcmax.pairwise.soiln <- data.frame(variable = "jmax25.vcmax25",
                                        treatment = "n.trt",
                                        cld(emmeans(jmax25.vcmax25, ~n.trt,
                                                    type = "response"),
                                            Letters = LETTERS))
jmax.vcmax.pairwise.inoc <- data.frame(variable = "jmax25.vcmax25",
                                        treatment = "inoc",
                                        cld(emmeans(jmax25.vcmax25, ~inoc,
                                                    type = "response"),
                                            Letters = LETTERS))
jmax.vcmax.pairwise <- jmax.vcmax.pairwise.full %>%
  full_join(jmax.vcmax.pairwise.soiln) %>%
  full_join(jmax.vcmax.pairwise.inoc) %>%
  dplyr::rename(emmean = response) %>%
  mutate(.group = trimws(.group, "both"),
         compact = c("B", "AB", "A", "A", "B", "A", "A", "A"))

##########################################################################
## Rd (standardized to 25 deg C)
##########################################################################
data$rd25[c(35, 63)] <- NA

rd <- lmer(log(rd25) ~ n.trt * inoc + (1 | block),
           data = data)

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
emmeans(rd, pairwise~n.trt) 
# Nitrogen addition increases dark respiration
emmeans(rd, pairwise~n.trt*inoc) 
# Low nitrogen, no inoculation has marginally lower Rd than high nitrogen,
# yes inoculation

# Emmean for fig making
rd.pairwise.full <- data.frame(variable = "rd25",
                               treatment = "full",
                               cld(emmeans(rd, ~n.trt*inoc, 
                                           type = "response"),
                                   Letters = LETTERS))
rd.pairwise.soiln <- data.frame(variable = "rd25",
                                treatment = "n.trt",
                                cld(emmeans(rd, ~n.trt, 
                                            type = "response"),
                                    Letters = LETTERS))
rd.pairwise.inoc <- data.frame(variable = "rd25",
                               treatment = "inoc",
                               cld(emmeans(rd, ~inoc, 
                                           type = "response"),
                                    Letters = LETTERS))
rd.pairwise <- rd.pairwise.full %>%
  full_join(rd.pairwise.soiln) %>%
  full_join(rd.pairwise.inoc) %>%
  dplyr::rename(emmean = response) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Rd:Vcmax (Vcmax is not standardized since Rd is not temp standardized)
##########################################################################
data$rd25.vcmax25[c(35)] <- NA

rd.vcmax <- lmer(rd25.vcmax25 ~ n.trt * inoc + (1 | block),
                 data = data)

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

# Emmean for fig making
rd.vcmax.pairwise.full <- data.frame(variable = "rd.vcmax",
                                     treatment = "full",
                                     cld(emmeans(rd.vcmax, ~n.trt*inoc),
                                         Letters = LETTERS))
rd.vcmax.pairwise.soiln <- data.frame(variable = "rd25.vcmax25",
                                      treatment = "n.trt",
                                      cld(emmeans(rd.vcmax, ~n.trt),
                                          Letters = LETTERS))
rd.vcmax.pairwise.inoc <- data.frame(variable = "rd25.vcmax25",
                                     treatment = "inoc",
                                     cld(emmeans(rd.vcmax, ~inoc),
                                         Letters = LETTERS))
rd.vcmax.pairwise <- rd.vcmax.pairwise.full %>%
  full_join(rd.vcmax.pairwise.soiln) %>%
  full_join(rd.vcmax.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Gs
##########################################################################
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
                                   Letters = LETTERS))
gs.pairwise.soiln <- data.frame(variable = "gs",
                                treatment = "n.trt",
                                cld(emmeans(gs, ~n.trt),
                                    Letters = LETTERS))
gs.pairwise.inoc <- data.frame(variable = "gs",
                               treatment = "inoc",
                               cld(emmeans(gs, ~inoc),
                                   Letters = LETTERS))
gs.pairwise <- gs.pairwise.full %>%
  full_join(gs.pairwise.soiln) %>%
  full_join(gs.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = c("B", "B", "A", "A", "B", "A", "A", "A"))

##########################################################################
## Ci:Ca
##########################################################################
data$ci.ca[c(23)] <- NA

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
                                Letters = LETTERS))
cica.pairwise.soiln <- data.frame(variable = "ci.ca",
                                  treatment = "n.trt",
                                  cld(emmeans(cica, ~n.trt),
                                      Letters = LETTERS))
cica.pairwise.inoc <- data.frame(variable = "ci.ca",
                                 treatment = "inoc",
                                 cld(emmeans(cica, ~inoc),
                                      Letters = LETTERS))
cica.pairwise <- cica.pairwise.full %>%
  full_join(cica.pairwise.soiln) %>%
  full_join(cica.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = c("A", "A", "A", "A", "B", "A", "A", "A"))

##########################################################################
## iWUE
##########################################################################
iwue <- lmer(log(iwue) ~ n.trt * inoc + (1 | block),
             data = data)

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
emmeans(iwue, pairwise~n.trt)
# Increasing nitrogen increases iWUE

# Emmean for fig making
iwue.pairwise.full <- data.frame(variable = "iwue",
                                 treatment = "full",
                                 cld(emmeans(iwue, ~n.trt*inoc,
                                             type = "response"),
                                     Letters = LETTERS))
iwue.pairwise.soiln <- data.frame(variable = "iwue",
                                  treatment = "n.trt",
                                  cld(emmeans(iwue, ~n.trt, 
                                              type = "response"),
                                      Letters = LETTERS))
iwue.pairwise.inoc <- data.frame(variable = "iwue",
                                 treatment = "inoc",
                                  cld(emmeans(iwue, ~inoc, 
                                              type = "response"),
                                      Letters = LETTERS))
iwue.pairwise <- iwue.pairwise.full %>%
  full_join(iwue.pairwise.soiln) %>%
  full_join(iwue.pairwise.inoc) %>%
  dplyr::rename(emmean = response) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Vcmax:gs
##########################################################################
vcmax.gs <- lmer(log(vcmax.gs) ~ n.trt * inoc + (1 | block), 
                 data = data)

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
# Increasing nitrogen increases Vcmax:gs through a relatively
# larger decrease in gs than Vcmax25

# Emmean for fig making
vcmax.gs.pairwise.full <- data.frame(variable = "vcmax.gs",
                                     treatment = "full",
                                     cld(emmeans(vcmax.gs, ~n.trt*inoc, 
                                                 type = "response"),
                                         Letters = LETTERS))
vcmax.gs.pairwise.soiln <- data.frame(variable = "vcmax.gs",
                                      treatment = "n.trt",
                                      cld(emmeans(vcmax.gs, ~n.trt, 
                                                  type = "response"),
                                          Letters = LETTERS))
vcmax.gs.pairwise.inoc <- data.frame(variable = "vcmax.gs",
                                      treatment = "inoc",
                                      cld(emmeans(vcmax.gs, ~inoc, 
                                                  type = "response"),
                                          Letters = LETTERS))
vcmax.gs.pairwise <- vcmax.gs.pairwise.full %>%
  full_join(vcmax.gs.pairwise.soiln) %>%
  full_join(vcmax.gs.pairwise.inoc) %>%
  dplyr::rename(emmean = response) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)


##########################################################################
## Specific leaf area
##########################################################################
sla <- lmer(sla ~ n.trt * inoc + (1 | block), 
                 data = data)

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
## Inoculated individuals generally have larger SLA than non-inoculated
## individuals

emmeans(sla, pairwise~n.trt)
## Increasing soil nitrogen has a marginally positive effect on SLA

# Emmean for fig making
sla.pairwise.full <- data.frame(variable = "sla",
                                treatment = "full",
                                cld(emmeans(sla, ~n.trt*inoc),
                                    Letters = LETTERS))
sla.pairwise.soiln <- data.frame(variable = "sla",
                                 treatment = "n.trt",
                                 cld(emmeans(sla,
                                             ~n.trt),
                                     Letters = LETTERS))
sla.pairwise.inoc <- data.frame(variable = "sla",
                                 treatment = "inoc",
                                 cld(emmeans(sla,
                                             ~inoc),
                                     Letters = LETTERS))
sla.pairwise <- sla.pairwise.full %>%
  full_join(sla.pairwise.soiln) %>%
  full_join(sla.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Focal leaf area
##########################################################################
fa <- lmer(focal.area ~ n.trt * inoc + (1 | block), 
            data = data)

# Check model assumptions
plot(fa)
qqnorm(residuals(fa))
qqline(residuals(fa))
hist(residuals(fa))
shapiro.test(residuals(fa))
outlierTest(fa)

# Model output
summary(fa)
Anova(fa)
r.squaredGLMM(fa)

# Pairwise comparisons
emmeans(fa, pairwise~inoc)
## Inoculated individuals generally have larger leaves than non-inoculated
## individuals

emmeans(fa, pairwise~n.trt)
## Increasing soil nitrogen generally increases leaf area

## Emmeans for fig making
fa.pairwise.full <- data.frame(variable = "focal.area",
                               treatment = "full",
                               cld(emmeans(fa, ~n.trt*inoc),
                                   Letters = LETTERS))
fa.pairwise.soiln <- data.frame(variable = "focal.area",
                                treatment = "n.trt",
                                cld(emmeans(fa,  ~n.trt),
                                    Letters = LETTERS))
fa.pairwise.inoc <- data.frame(variable = "focal.area",
                               treatment = "inoc",
                               cld(emmeans(fa, ~inoc),
                                    Letters = LETTERS))
fa.pairwise <- fa.pairwise.full %>%
  full_join(fa.pairwise.soiln) %>%
  full_join(fa.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Focal leaf biomass (dry)
##########################################################################
focal.bio <- lmer(dry.biomass ~ n.trt * inoc + (1 | block), 
                  data = data)

# Check model assumptions
plot(focal.bio)
qqnorm(residuals(focal.bio))
qqline(residuals(focal.bio))
hist(residuals(focal.bio))
shapiro.test(residuals(focal.bio))
outlierTest(focal.bio)

# Model output
summary(focal.bio)
Anova(focal.bio)
r.squaredGLMM(focal.bio)

# Pairwise comparisons
emmeans(focal.bio, pairwise~n.trt)
## Increasing soil nitrogen generally increases dry biomass

## NOTE: the marginal impact of n.trt on SLA was driven by a marginally
## larger increase in leaf area than dry biomass (but both leaf area and
## dry biomass increased)

fb.pairwise.full <- data.frame(variable = "focal.biomass",
                               treatment = "full",
                               cld(emmeans(focal.bio, ~n.trt*inoc),
                                   Letters = LETTERS))
fb.pairwise.soiln <- data.frame(variable = "focal.biomass",
                                treatment = "n.trt",
                                cld(emmeans(focal.bio, ~n.trt),
                                    Letters = LETTERS))
fb.pairwise.inoc <- data.frame(variable = "focal.biomass",
                               treatment = "inoc",
                               cld(emmeans(focal.bio, ~inoc),
                                   Letters = LETTERS))
fb.pairwise <- fb.pairwise.full %>%
  full_join(fb.pairwise.soiln) %>%
  full_join(fb.pairwise.inoc) %>%
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
                                   Letters = LETTERS))
tla.pairwise.soiln <- data.frame(variable = "total.leaf.area",
                                treatment = "n.trt",
                                cld(emmeans(tla, ~n.trt),
                                    Letters = LETTERS))
tla.pairwise.inoc <- data.frame(variable = "total.leaf.area",
                               treatment = "inoc",
                               cld(emmeans(tla, ~inoc),
                                   Letters = LETTERS))
tla.pairwise <- tla.pairwise.full %>%
  full_join(tla.pairwise.soiln) %>%
  full_join(tla.pairwise.inoc) %>%
  mutate(.group = trimws(.group, "both"),
         compact = .group)

##########################################################################
## Make merged emmeans file
##########################################################################
comp.letters <- a400.pairwise %>%
  full_join(vcmax.pairwise) %>%
  full_join(jmax.pairwise) %>%
  full_join(jmax.vcmax.pairwise) %>%
  full_join(rd.pairwise) %>%
  full_join(rd.vcmax.pairwise) %>%
  full_join(gs.pairwise) %>%
  full_join(cica.pairwise) %>%
  full_join(iwue.pairwise) %>%
  full_join(vcmax.gs.pairwise) %>%
  full_join(sla.pairwise) %>%
  full_join(fa.pairwise) %>%
  full_join(fb.pairwise) %>%
  full_join(tla.pairwise) %>%
  dplyr::rename(comparison = treatment) %>%
  unite("treatment", n.trt:inoc, remove = "FALSE") %>%
  mutate(treatment = factor(treatment, levels = c("LN_NI", "HN_NI",
                                                  "LN_YI", "HN_YI")),
         compact = tolower(compact)) %>%
  data.frame()
comp.letters

## Write pairwise comparison csv file
write.csv(comp.letters, "../data/2021NxI_compact_letters.csv")
