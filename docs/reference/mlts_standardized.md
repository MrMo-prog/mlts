# Get Standardized Estimates for an mlts Model

Get Standardized Estimates for an mlts Model

## Usage

``` r
mlts_standardized(
  object,
  what = c("between", "within", "both"),
  digits = 3,
  prob = 0.95,
  add_cluster_std = FALSE,
  get_samples = FALSE
)
```

## Arguments

- object:

  `mltsfit`. Output of [`mlts_model`](mlts_model.md) and related
  functions.

- what:

  character. Get between-level standardized estimates
  (`what = "between"`, the default), within-level standardized estimates
  averaged over clusters (`what = "within"`), or both (`what = "both"`).

- digits:

  Number of digits. Default is 3.

- prob:

  A value between 0 and 1 to indicate the width of the credible
  interval. Default is .95.

- add_cluster_std:

  logical. If `what = "within"`, within-level standardized effects for
  each cluster are included in the output (defaults to `FALSE`).

- get_samples:

  logical. For internal use.

## Value

A `list` containing between- and within-level standardized parameters.

## Examples

``` r
# \donttest{
# build simple vector-autoregressive mlts model for two time-series variables
var_model <- mlts_model(q = 2)

# fit model with (artificial) dataset ts_data
fit <- mlts_fit(
  model = var_model,
  data = ts_data,
  ts = c("Y1", "Y2"), # time-series variables
  id = "ID", # identifier variable
  time = "time", # time variable
  tinterval = 1, # interval for approximation of continuous-time dynamic model,
  monitor_person_pars = TRUE # person parameters need to be sampled for standardization
)

# inspect standardized parameter estimates
mlts_standardized(fit)
# }
```
