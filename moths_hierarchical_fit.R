library(readr)
library(cmdstanr)
library(posterior)
library(tidyverse)

# ── 0. Load data ─────────────────────────────────────────────────────────────

Maine.clean  <- read_rds("Maine_clean.rds")
flight_cumul <- read_rds("flight_cumul.rds")

# Keep only site-years with a well-defined flight period (within 1 SD of mean)
Maine.clean <- Maine.clean |>
  semi_join(
    flight_cumul |> filter(duration > 6, duration < 19),
    by = "id"
  )

N <- nrow(Maine.clean)
T <- 365L
K <- 30L   # max kernel lag (days); covers >99% of any reasonable gamma kernel

# ── 1. Helper ─────────────────────────────────────────────────────────────────

# Expand a sparse DOY tibble to a full 365-element vector, filling gaps with 0
expand_to_365 <- function(tb, col) {
  tibble(DOY = 1:365) |>
    left_join(tb, by = "DOY") |>
    mutate({{ col }} := replace_na(pull(pick({{ col }}), {{ col }}), 0)) |>
    pull({{ col }})
}

# ── 2. Site index ─────────────────────────────────────────────────────────────

# Integer site ID: multiple site-years share the same site index,
# enabling partial pooling of kernel parameters across years at each location
site_levels <- unique(Maine.clean$Location)
site_id     <- match(Maine.clean$Location, site_levels)
S           <- length(site_levels)

# ── 3. Pupal PMF matrix [N x 365] ────────────────────────────────────────────

# Each row is the pupal count series for one site-year, expanded to 365 days
# and normalised to a PMF (the input signal to be convolved with the kernel)
pup_mat <- Maine.clean |>
  rowwise() |>
  summarise(pup = list(expand_to_365(Pupae, Pupae)), .groups = "drop") |>
  pull(pup) |>
  (\(lst) {
    m <- matrix(0, nrow = N, ncol = T)
    for (i in seq_along(lst)) m[i, ] <- lst[[i]]
    m
  })()

row_sums  <- rowSums(pup_mat)
zero_rows <- row_sums == 0
pup_pmf   <- pup_mat
pup_pmf[!zero_rows, ] <- pup_mat[!zero_rows, ] / row_sums[!zero_rows]
pup_pmf[ zero_rows, ] <- 1 / T   # flat fallback; scale will be near zero

# ── 4. Observed moth count matrix [N x 365] ──────────────────────────────────

# Raw integer counts are passed directly to Stan. The Dirichlet-multinomial
# likelihood operates on counts, so no normalisation is needed here.
# Stan will compute row totals and the implied PMF internally.
obs_mat <- Maine.clean |>
  rowwise() |>
  summarise(obs = list(expand_to_365(Moths, Moths)), .groups = "drop") |>
  pull(obs) |>
  (\(lst) {
    m <- matrix(0L, nrow = N, ncol = T)
    for (i in seq_along(lst)) m[i, ] <- as.integer(lst[[i]])
    m
  })()

# ── 5. Temperature matrix [N x 365] ──────────────────────────────────────────

# Daily maximum temperature per site-year, used to gate flight activity:
# moths do not fly on days below an estimated threshold temperature
tavg_mat <- Maine.clean |>
  rowwise() |>
  summarise(
    tavg = list({
      tb <- Weather |>
        mutate(DOY = as.integer(format(Date, "%j"))) |>
        group_by(DOY) |>
        summarise(tavg = mean(tavg, na.rm = TRUE), .groups = "drop")
      expand_to_365(tb, tavg)
    }),
    .groups = "drop"
  ) |>
  pull(tavg) |>
  (\(lst) {
    m <- matrix(0, nrow = N, ncol = T)
    for (i in seq_along(lst)) m[i, ] <- lst[[i]]
    m
  })()

# ── 6. Assemble Stan data ─────────────────────────────────────────────────────

stan_data <- list(
  N       = N,
  S       = S,
  T       = T,
  K       = K,
  site_id = site_id,
  pup_pmf = pup_pmf,
  y       = obs_mat,
  tavg    = tavg_mat
)

# ── 7. Compile and sample ─────────────────────────────────────────────────────

mod <- cmdstan_model("moths_hierarchical.stan")

fit <- mod$sample(
  data            = stan_data,
  seed            = 1965,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  adapt_delta     = 0.95,
  max_treedepth   = 12,
  refresh         = 200,
  output_dir      = "chain_outputs"
)

fit$save_object("fit_moths.rds")

fit <- read_rds("fit_moths.rds")

# ── 8. Diagnostics ────────────────────────────────────────────────────────────

fit$print(variables = c(
  "mu_log_shape", "sigma_log_shape",
  "mu_log_rate",  "sigma_log_rate",
  "T_thresh",     "T_steep",
  "concentration"
))

fit$cmdstan_diagnose()

# ── 9. Extract timing summaries ───────────────────────────────────────────────

draws <- fit$draws(format = "df")

# Fix 1: dur_summary — extract the number inside brackets, not the first digits
dur_summary <- draws |>
  select(starts_with("duration_80")) |>
  pivot_longer(everything(), names_to = "par") |>
  mutate(i = as.integer(str_extract(par, "(?<=\\[)\\d+(?=\\])"))) |>
  group_by(i) |>
  summarise(
    dur_mean = mean(value),
    dur_q05  = quantile(value, 0.05),
    dur_q95  = quantile(value, 0.95),
    .groups  = "drop"
  )

# Fix 2: peak_summary — same fix for safety, and ungroup Maine.clean before mutating
peak_summary <- draws |>
  select(starts_with("peak_doy")) |>
  pivot_longer(everything(), names_to = "par") |>
  mutate(i = as.integer(str_extract(par, "(?<=\\[)\\d+(?=\\])"))) |>
  group_by(i) |>
  summarise(
    peak_mean = mean(value),
    peak_q05  = quantile(value, 0.05),
    peak_q95  = quantile(value, 0.95),
    .groups   = "drop"
  )

# Fix 3: ungroup before row_number() to get a global 1:N index
results <- Maine.clean |>
  ungroup() |>
  mutate(i = row_number()) |>
  left_join(peak_summary, by = "i") |>
  left_join(dur_summary,  by = "i") |>
  select(Location, Year, peak_mean, peak_q05, peak_q95,
         dur_mean, dur_q05, dur_q95)

print(results)

# ── 10. Diagnostic plots ──────────────────────────────────────────────────────

# Temperature threshold sigmoid curve
tibble(temp = seq(-5, 30, 0.1)) |>
  mutate(
    flight_prob = plogis(mean(draws$T_steep) * (temp - mean(draws$T_thresh)))
  ) |>
  ggplot(aes(temp, flight_prob)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = mean(draws$T_thresh), linetype = "dashed") +
  labs(x = "Daily max temperature (°C)",
       y = "Relative flight probability",
       title = "Estimated temperature threshold for moth flight")

# Peak timing by site and year
results |>
  ggplot(aes(Year, peak_mean, colour = Location)) +
  geom_point() +
  geom_errorbar(aes(ymin = peak_q05, ymax = peak_q95), width = 0, alpha = 0.4) +
  labs(x = "Year", y = "Posterior mean peak DOY",
       title = "Emergence peak timing by site and year")

