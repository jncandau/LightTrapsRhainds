library(dplyr)
library(tidyr)
library(lubridate)

# Read the data
temp_data <- read.csv2("2m_temperature_at_sites.csv")

# Convert date column to Date type (adjust column name as needed)
temp_data <- temp_data |>
  mutate(date = as.Date(datetime),
  DOY = yday(date)) |> 
  dplyr::filter(DOY < 182)

# Calculate daily mean temperature for each id
daily_means <- temp_data |>
  group_by(id,DOY) |>
  summarize(daily_mean_temp = mean(hourly_temperature_2m, na.rm = TRUE))

# Pivot to wide format with temp_1, temp_2, etc.
wide_temp <- daily_means |>
  pivot_wider(
    id_cols = id,
    names_from = DOY,
    values_from = daily_mean_temp,
    names_prefix = "temp."
  )

# Get the number of temp columns
n_days <- ncol(wide_temp) - 1

# Add gap columns (initialized as NA or 0)
lag_cols <- setNames(
  as.data.frame(matrix(rep(0:(n_days-1),each=nrow(wide_temp)), nrow = nrow(wide_temp), ncol = n_days),byrow = FALSE),
  paste0("lag.", 1:n_days)
)

# Combine temp and gap columns
result <- bind_cols(wide_temp, lag_cols)

result


