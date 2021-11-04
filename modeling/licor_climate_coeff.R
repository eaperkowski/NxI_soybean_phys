## Load libraries
library(dplyr)
library(tidyr)
library(bigleaf)
library(lubridate)
library(ggplot2)
library(car)
library(emmeans)

## Set working directory
setwd("/Users/eaperkowski/git/joseph_greenhouse_phys_2021/modeling/")

## Read LI-COR and Mesowest .csv file
df.licor <- read.csv("par_data_licor.csv")
df.climate <- read.csv("climate_data_mesowest.csv")
head(df.climate)

## Clean LI-COR data to only include date, time and PAR data
df.licor <- df.licor %>%
  select(date, 
         solar.rad.licor = Qamb_out) %>%
  mutate(date = round_date(ymd_hms(strptime(date, 
                                            format = "%Y%m%d %H:%M:%S")),
                           "15 minutes"),
         solar.rad.licor = ifelse(solar.rad.licor > 0.5, solar.rad.licor, 0)) %>%
  slice(-1)
head(df.licor)

## Clean climate data to only include date, time and PAR data
df.climate <- df.climate %>%
  select(date = Date_Time, 
         solar.rad = solar_radiation_set_1) %>%
  mutate(date = round_date(ymd_hms(strptime(date, 
                                            format = "%m/%d/%Y %H:%M")),
                           "15 minutes"),
         solar.rad = Rg.to.PPFD(as.numeric(solar.rad))) %>%
  group_by(date) %>%
  summarize(solar.rad.climate = mean(solar.rad, na.rm = TRUE)) %>%
  filter(date > "2021-09-01" & date <"2021-11-02") %>%
  data.frame()

df <- df.licor %>%
  inner_join(df.climate)

ggplot(data = df, aes(x = date,
                      y = as.numeric(solar.rad.climate))) +
  geom_point() +
  geom_smooth(method = 'loess') +
  scale_y_continuous(limits = c(0, 1800),
                     breaks = c(0, 1800, 600))

test <- as.data.frame(spline(df$date, df$solar.rad.licor))

ggplot(data = df, aes(x = date,
                      y = as.numeric(solar.rad.licor))) +
  geom_point() +
  geom_line() +
  scale_y_continuous(limits = c(0, 1800),
                     breaks = c(0, 1800, 600))



test <- lm(solar.rad.licor ~ solar.rad.climate, data = df)
summary(test)
Anova(test)  

emtrends(test, ~1, var = "solar.rad.climate")
