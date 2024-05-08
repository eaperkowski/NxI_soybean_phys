library(tidyverse)

block1n <- read.csv("git/NxI_soybean_phys/hobo_temp/Joseph_block1N.csv") %>%
  mutate(date_time = lubridate::mdy_hm(date)) %>%
  select(date_time, temp_block1n = temp)

block1s <- read.csv("git/NxI_soybean_phys/hobo_temp/Joseph_block1S.csv") %>%
  mutate(date_time = lubridate::mdy_hm(date)) %>%
  select(date_time, temp_block1s = temp)

block2n <- read.csv("git/NxI_soybean_phys/hobo_temp/Joseph_block2N.csv") %>%
  mutate(date_time = lubridate::mdy_hm(date)) %>%
  select(date_time, temp_block2n = temp)

block2s <- read.csv("git/NxI_soybean_phys/hobo_temp/Joseph_block2S.csv") %>%
  mutate(date_time = lubridate::mdy_hm(date)) %>%
  select(date_time, temp_block2s = temp)

block3n <- read.csv("git/NxI_soybean_phys/hobo_temp/Joseph_block3N.csv") %>%
  mutate(date_time = lubridate::mdy_hm(date)) %>%
  select(date_time, temp_block3n = temp)

block3s <- read.csv("git/NxI_soybean_phys/hobo_temp/Joseph_block3S.csv") %>%
  mutate(date_time = lubridate::mdy_hm(date)) %>%
  select(date_time, temp_block3s = temp)

block4n <- read.csv("git/NxI_soybean_phys/hobo_temp/Joseph_block4N.csv") %>%
  mutate(date_time = lubridate::mdy_hm(date)) %>%
  select(date_time, temp_block4n = temp)

block4s <- read.csv("git/NxI_soybean_phys/hobo_temp/Joseph_block4S.csv") %>%
  mutate(date_time = lubridate::mdy_hm(date)) %>%
  select(date_time, temp_block4s = temp)




temp_all <- block1n %>%
  full_join(block1s, by = "date_time") %>%
  full_join(block2n, by = "date_time") %>%
  full_join(block2s, by = "date_time") %>%
  full_join(block3n, by = "date_time") %>%
  full_join(block3s, by = "date_time") %>%
  full_join(block4n, by = "date_time") %>%
  full_join(block4s, by = "date_time")

temp_day_summary <- temp_all %>%
  pivot_longer(cols = temp_block1n:temp_block4s,
               names_to = "block", values_to = "temp") %>%
  mutate(date = date(date_time)) %>%
  group_by(date) %>%
  summarize(max_temp = max(temp),
            min_temp = min(temp))


temp_mean_summary <- temp_day_summary %>%
  summarize(max_exp_temp = mean(max_temp, na.rm = TRUE),
            sd_max_exp_temp = sd(max_temp, na.rm = TRUE),
            min_exp_temp = mean(min_temp, na.rm = TRUE),
            sd_min_exp_temp = sd(min_temp, na.rm = TRUE))

