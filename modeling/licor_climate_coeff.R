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
  dplyr::select(date, 
         solar.rad.licor = Qamb_out) %>%
  mutate(date = round_date(ymd_hms(strptime(date, 
                                            format = "%Y%m%d %H:%M:%S")),
                           "15 minutes"),
         solar.rad.licor = ifelse(solar.rad.licor > 0.5, solar.rad.licor, 0),
         keep.row = "yes") %>%
  slice(-1)

df.licor$keep.row[c(2, 3, 9, 11, 12, 18, 82, 84, 88, 89, 90, 91,
                    93, 98, 99, 105, 106, 107, 108, 178, 180, 185,
                    186, 187, 189, 194, 195, 201, 202, 203, 204,
                    274, 276, 278, 281, 282, 283, 285, 290, 291,
                    297, 298, 299, 300)] <- "no"
df.licor



## Clean climate data to only include date, time and PAR data
df.climate <- df.climate %>%
  dplyr::select(date = Date_Time, 
         solar.rad = solar_radiation_set_1) %>%
  mutate(date = round_date(ymd_hms(strptime(date, 
                                            format = "%m/%d/%Y %H:%M")),
                           "15 minutes"),
         solar.rad = Rg.to.PPFD(as.numeric(solar.rad))) %>%
  group_by(date) %>%
  summarize(solar.rad.climate = mean(solar.rad, na.rm = TRUE)) %>%
  filter(date > "2021-09-16" & date <"2021-11-02") %>%
  data.frame()

df <- df.licor %>%
  inner_join(df.climate) %>%
  mutate(solar.rad.licor = ifelse(keep.row == "no", NA, solar.rad.licor))

## Visualize coefficient plot: log-log transformations seem to give slope
## similar to 1
ggplot(data = subset(df, solar.rad.climate != 0), 
       aes(x = log(as.numeric(solar.rad.climate)),
           y = log(as.numeric(solar.rad.licor)))) +
  geom_point() +
  geom_smooth(method = 'lm', size = 2) +
  scale_x_continuous(limits = c(0, 8),
                     breaks = seq(0, 8, 2)) +
  scale_y_continuous(limits = c(0, 8),
                     breaks = seq(0, 8, 2)) +
  geom_segment(aes(x = 1.278202, xend = 7.4, y = 1.278202,
                   yend = 7.4), size = 2, color = "red", linetype = "dashed") +
  labs(x = expression("log PAR"["clim"]~"(μmol m"^"-2"~"s"^"-1"~")"),
       y = expression("log PAR"["LI-6800"]~"(μmol m"^"-2"~"s"^"-1"~")"))



## Coefficient model
coef.mod <- lm(log(as.numeric(solar.rad.licor)) ~ log(solar.rad.climate),
               data = subset(df, solar.rad.climate != 0))

## Test model assumptions
plot(coef.mod)
qqnorm(residuals(coef.mod))
qqline(residuals(coef.mod))
hist(residuals(coef.mod))
shapiro.test(residuals(coef.mod))
outlierTest(coef.mod)

## Model output to get coefficients
summary(coef.mod)

## Derive helper function to estimate greenhouse PAR
ghpar <- function(clim.par) {
  ghpar = clim.par ^ 1.09949 * exp(-1.40537)
  
  return(ghpar)
}

## Add greenhouse PAR column
df.climate <- df.climate %>%
  mutate(ghpar = ghpar(solar.rad.climate)) %>%
  select(date, ghpar) %>%
  filter(date > "2021-10-01" & date < "2021-10-14 00:00:00")

max(df.climate$ghpar)
df.climate

## Check greenhouse plot for any issues
ggplot(data = df.climate,
       aes(x = date,
           y = ghpar)) +
  geom_point(size = 1.5) +
  scale_y_continuous(limits = c(0, 1200),
                     breaks = seq(0, 1200, 400)) +
  theme_bw()
