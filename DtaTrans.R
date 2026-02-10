library(tidyverse)
library(tseries)

# Load raw data
Wdta <- read.csv('data/WindSpdDr.csv')

# Clean and transform to wind vectors
df <- Wdta %>% 
  select(1:(ncol(.) - 1)) %>% 
  select(-c(3, 4)) %>% 
  slice(-1) %>% 
  {
    set_names(., .[1,]) %>% 
      slice(-1)
  } %>% 
  rename(
    date = time,
    Wsp = 'wind_speed_10m (km/h)',
    Wdr = 'wind_direction_10m (°)'
  ) %>% 
  mutate(
    date = as.POSIXct(date, format = '%Y-%m-%dT%H:%M', tz = 'UTC'),
    Wsp = as.numeric(Wsp),
    Wdr = as.numeric(Wdr)
  ) %>% 
  mutate(
    # Convert direction to radians
    WdrRad = (Wdr %% 360) * pi / 180,
    u = Wsp * cos(WdrRad), # east-west component
    v = Wsp * sin(WdrRad), # north-south component
    # Time features for temporal analysis
    hour = hour(date),
    day = day(date),
    week = week(date),
    month = month(date)
  )

# Compute summary statistics
SummaryStats <- df %>% 
  summarise(
    'N observations' = n(),
    'Mean speed (km/h)' = mean(Wsp),
    'Sd speed (km/h)' = sd(Wsp),
    'Min speed (km/h)' = min(Wsp),
    'Max speed (km/h)' = max(Wsp),
    'Mean u (km/h)' = mean(u),
    'Mean v (km/h)' = mean(v),
    'Sd u (km/h)' = sd(u),
    'Sd v (km/h)' = sd(v),
    'Correlation(u,v) (km/h)' = cor(u,v),
  ) %>% 
  pivot_longer(
    everything(),
    names_to = 'Statistic',
    values_to = 'Value'
  ) %>% 
  mutate(
    Value = if_else(
      Statistic == 'N observations',
      Value,
      round(Value, 2)
    )
  )

# Test wind speed stationarity
AdfWsp <- adf.test(df$Wsp)
AdfU <- adf.test(df$u)
AdfV <- adf.test(df$v)

# Test results table
AdfRes <- data.frame(
  Variable = c('Wind speed', 'u component', 'v component'),
  Statistic = c(AdfWsp$statistic, AdfU$statistic, AdfV$statistic),
  pvalue = c(AdfWsp$p.value, AdfU$p.value, AdfV$p.value)
)

# Compute monthly statistics
MonthStats <- df %>% 
  group_by(month) %>% 
  summarise(
    MeanSpd = mean(Wsp),
    SdSpd = sd(Wsp),
    MeanU = mean(u),
    MeanV = mean(v),
    n = n()
  )
