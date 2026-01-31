library(tidyverse)
library(openmeteo)

df <- get_historical_weather(latitude = 44.4,longitude = -65.2,start_date = "1970-01-01",end_date = "2020-01-01",variables = "temperature_2m")

degree_days <- df %>%
  # Extract date components

  mutate(
    year = year(datetime),
    month = month(datetime),
    day = day(datetime)
  ) %>%
  # Filter to Jan 1 – July 15 window
  filter(month < 7 | (month == 7 & day <= 15)) %>%
  # Calculate hourly degree contribution (only when temp > 5°C)
  mutate(
    dd_hourly = pmax(prediction - 5, 0)
  ) %>%
  # Sum by year and convert hourly to daily (divide by 24)
  group_by(year) %>%
  summarise(
    degree_days_above_5 = sum(dd_hourly) / 24,
    .groups = "drop"
  )

# Plot with trend line
ggplot(degree_days, aes(x = year, y = degree_days_above_5)) +
  geom_point(color = "steelblue", size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  labs(
    x = "Year",
    y = "Degree Days (>5°C)",
    title = "Accumulated Degree Days (Jan 1 – July 15)",
    subtitle = "Base temperature: 5°C"
  ) +
  theme_minimal()

# Fit linear model to test for trend
trend_model <- lm(degree_days_above_5 ~ year, data = degree_days)
summary(trend_model)

# For a cleaner summary
broom::tidy(trend_model)
broom::glance(trend_model)