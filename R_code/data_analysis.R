## Libraries
library(lme4)
library(emmeans)
library(car)
library(tidyverse)
library(dplyr)
library(MuMIn)
library(multcomp)
library(multcompView)

## Remove digit bounds in emmeans package
emm_options(opt.digits = FALSE)

## Load data
data <- read.csv("../data/trait_data.csv",
                 na.strings = "NA")

## Add Rd:Vcmax, Jmax25:Vcmax25, Vcmax:gs column
data$Rd.Vcmax <- data$Rd.TPU / data$Vcmax.TPU
data$Jmax25.Vcmax25 <- data$Jmax25.TPU / data$Vcmax25.TPU
data$Vcmax.gs <-  data$Vcmax.TPU / data$gsw

head(data)

##########################################################################
## A400
##########################################################################
a400 <- lmer(A ~ n.trt * inoc + (1 | block),
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

##########################################################################
## Vcmax25 with TPU
##########################################################################
vcmax25 <- lmer(Vcmax25.TPU ~ n.trt * inoc + (1 | block),
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

##########################################################################
## Jmax25 with TPU
##########################################################################
jmax25 <- lmer(Jmax25.TPU ~ n.trt * inoc + (1 | block),
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

# Pairwise
emmeans(jmax25, pairwise~n.trt) # Nitrogen addition decreases Jmax25?

##########################################################################
## Jmax:Vcmax25 with TPU
##########################################################################
data$Jmax25.Vcmax25[c(17, 62)] <- NA

jmax25.vcmax25 <- lmer(log(Jmax25.Vcmax25) ~ n.trt * inoc + (1 | block),
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

##########################################################################
## Rd 
##########################################################################
data$Rd.TPU[c(2, 7, 31)] <- NA

rd <- lmer(log(Rd.TPU) ~ n.trt * inoc + (1 | block),
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

##########################################################################
## Rd:Vcmax (Vcmax is not standardized since Rd is not temp standardized)
##########################################################################
data$Rd.Vcmax[c(7, 31)] <- NA

rd.vcmax <- lmer(Rd.Vcmax ~ n.trt * inoc + (1 | block),
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

##########################################################################
## Gs
##########################################################################
gs <- lmer(gsw ~ n.trt * inoc + (1 | block),
           data = data)

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

##########################################################################
## Ci:Ca
##########################################################################
cica <- lmer(ci.ca ~ n.trt * inoc + (1 | block),
             data = data)

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

##########################################################################
## Vcmax:gs
##########################################################################
vcmax.gs <- lmer(log(Vcmax.gs) ~ n.trt * inoc + (1 | block), 
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


