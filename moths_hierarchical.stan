// Hierarchical convolution model for spruce budworm moth trapping data
//
// Goal: infer the TIMING and DURATION of adult moth emergence, not abundance.
//
// Model structure:
//   - The pupal PMF (pup_pmf) is convolved with a site-level gamma kernel
//     to produce a predicted emergence PMF.
//   - A sigmoid temperature gate suppresses predicted flight on cold days,
//     after which the result is renormalised to a valid PMF (pred_pmf).
//   - The observed moth PMF (moth_pmf) is modelled as a Dirichlet-multinomial
//     draw from pred_pmf, weighted by total observed moths (total_obs).
//
// Hierarchy:
//   Population
//     └── Site (log_shape_s, log_rate_s): partial pooling over Locations
//           └── Site-year: inherits site kernel; only timing varies via
//                          the temperature gate (tmax differs by year)
//
// Parameters NOT in this model:
//   - log_scale  (abundance; irrelevant when modelling shape only)
//   - size       (NegBin overdispersion; replaced by Dirichlet concentration)
//   - beta_shape / beta_rate  (GDD regression; removed)

data {
  int<lower=1> N;                            // number of site-years
  int<lower=1> S;                            // number of unique Locations
  int<lower=1> T;                            // days per year (365)
  int<lower=1> K;                            // max kernel lag (days)

  array[N] int<lower=1, upper=S> site_id;   // site index for each site-year

  matrix[N, T] pup_pmf;                     // normalised pupal PMF [N x T]
  array[N, T] int<lower=0> y;               // observed moth counts [N x T]
  matrix[N, T] tmax;                        // daily max temperature [N x T], degrees C
}

transformed data {
  // Row totals: total moths observed per site-year.
  // Used to scale the Dirichlet concentration so that site-years with more
  // data exert proportionally more influence on the likelihood.
  array[N] int total_obs;
  for (i in 1:N) {
    int s = 0;
    for (t in 1:T) s += y[i, t];
    total_obs[i] = s;
  }
}

parameters {
  // ── Population hyperparameters ──────────────────────────────────────────────
  real mu_log_shape;                 // population mean of log(shape)
  real<lower=0> sigma_log_shape;     // population SD   of log(shape)
  real mu_log_rate;                  // population mean of log(rate)
  real<lower=0> sigma_log_rate;      // population SD   of log(rate)

  // ── Site-level random effects (non-centred parameterisation) ────────────────
  vector[S] z_shape;                 // N(0,1) offsets -> site log(shape)
  vector[S] z_rate;                  // N(0,1) offsets -> site log(rate)

  // ── Temperature threshold for flight ────────────────────────────────────────
  real T_thresh;                     // threshold temperature (degrees C)
  real<lower=0> T_steep;             // sigmoid steepness; large -> near-binary gate

  // ── Dirichlet-multinomial concentration ─────────────────────────────────────
  real<lower=0> concentration;       // global precision scalar
}

transformed parameters {
  // ── Site-level kernel parameters (partial pooling over locations) ────────────
  vector[S] log_shape_s = mu_log_shape + sigma_log_shape * z_shape;
  vector[S] log_rate_s  = mu_log_rate  + sigma_log_rate  * z_rate;

  // ── Predicted emergence PMF via convolution + temperature gating ─────────────
  //
  // For each site-year i and day t:
  //   conv[t]       = sum_{k=1}^{K} pup_pmf[i, t-k] * kernel[s, k]
  //   gated[t]      = conv[t] * sigmoid(T_steep * (tmax[i,t] - T_thresh))
  //   pred_pmf[i,t] = gated[t] / sum(gated)
  //
  // Site-years at the same location share the same kernel; they differ only
  // through their year-specific tmax values, which shift the temperature gate.
  //
  matrix[N, T] pred_pmf;
  for (i in 1:N) {
    int s = site_id[i];

    // Discrete gamma PMF kernel of length K, normalised
    vector[K] kv;
    for (k in 1:K)
      kv[k] = exp(gamma_lpdf(k | exp(log_shape_s[s]), exp(log_rate_s[s])));
    vector[K] kernel = kv / sum(kv);

    // Convolution with temperature gate applied simultaneously
    vector[T] gated;
    for (t in 1:T) {
      real cv = 0;
      for (k in 1:K) {
        int src = t - k;
        if (src >= 1) cv += pup_pmf[i, src] * kernel[k];
      }
      real gate = inv_logit(T_steep * (tmax[i, t] - T_thresh));
      gated[t]  = cv * gate;
    }

    // Add a small uniform floor before renormalising so that no day has
    // exactly zero probability. This prevents zero alpha values in the
    // Dirichlet-multinomial even when the temperature gate suppresses
    // entire regions of the season or the pupal signal is absent.
    vector[T] floored = gated + 1e-6;
    pred_pmf[i] = to_row_vector(floored / sum(floored));
  }
}

model {
  // ── Hyperpriors ──────────────────────────────────────────────────────────────
  // shape ~ exp(1.5) ≈ 4.5; rate ~ exp(-0.5) ≈ 0.6 -> mean lag ≈ 7.5 days
  mu_log_shape    ~ normal(1.5, 1);
  sigma_log_shape ~ normal(0, 0.5);
  mu_log_rate     ~ normal(-0.5, 1);
  sigma_log_rate  ~ normal(0, 0.5);

  // ── Site random effects ──────────────────────────────────────────────────────
  z_shape ~ std_normal();
  z_rate  ~ std_normal();

  // ── Temperature threshold ────────────────────────────────────────────────────
  T_thresh ~ normal(10, 5);    // prior centred on 10 degrees C
  T_steep  ~ normal(1,  0.5);  // moderate sharpness

  // ── Dirichlet-multinomial concentration ──────────────────────────────────────
  // Larger concentration = tighter fit of pred_pmf to moth_pmf.
  // Scaled by total_obs so site-years with more data contribute more.
  concentration ~ gamma(2, 0.01);

  // ── Likelihood ───────────────────────────────────────────────────────────────
  // Dirichlet-multinomial: observed counts y[i] are a draw from pred_pmf[i],
  // with precision scaled by total_obs[i] so data-rich site-years dominate.
  // The floor on alpha ensures all elements are strictly positive, as required.
  for (i in 1:N) {
    vector[T] alpha = concentration * total_obs[i] * to_vector(pred_pmf[i]);
    target += dirichlet_multinomial_lpmf(y[i] | fmax(alpha, 1e-6));
  }
}

generated quantities {
  // ── Posterior predictive draws ────────────────────────────────────────────────
  array[N, T] int y_rep;
  for (i in 1:N) {
    vector[T] alpha = concentration * total_obs[i] * to_vector(pred_pmf[i]);
    array[T] int counts = dirichlet_multinomial_rng(fmax(alpha, 1e-6), total_obs[i]);
    for (t in 1:T) y_rep[i, t] = counts[t];
  }

  // ── Peak emergence DOY ────────────────────────────────────────────────────────
  // Day of year with the highest predicted emergence probability
  array[N] int peak_doy;
  for (i in 1:N) {
    real best_p = -1;
    int  best_t = 1;
    for (t in 1:T) {
      if (pred_pmf[i, t] > best_p) {
        best_p = pred_pmf[i, t];
        best_t = t;
      }
    }
    peak_doy[i] = best_t;
  }

  // ── Duration: shortest interval containing 80% of emergence mass ─────────────
  // Computed as the number of days needed when accumulating from the highest
  // probability day downward (highest density interval approximation)
  vector[N] duration_80;
  for (i in 1:N) {
    vector[T] p_sort = sort_desc(to_vector(pred_pmf[i]));
    real cummass = 0;
    int  width   = 0;
    while (cummass < 0.8 && width < T) {
      width   += 1;
      cummass += p_sort[width];
    }
    duration_80[i] = width;
  }
}
