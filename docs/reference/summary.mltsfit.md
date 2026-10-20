# Create a summary of a fitted model with class `mltsfit`

Create a summary of a fitted model with class `mltsfit`

## Usage

``` r
# S3 method for class 'mltsfit'
summary(
  object,
  priors = FALSE,
  se = FALSE,
  prob = 0.95,
  bpe = c("mean"),
  digits = 3,
  flag_signif = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class `mltsfit`.

- priors:

  Add prior information (default = FALSE).

- se:

  Logical. Should the Monte Carlo Standard Error be included in the
  summary? Defaults to `FALSE`.

- prob:

  A value between 0 and 1 to indicate the width of the credible
  interval. Default is .95.

- bpe:

  Bayesian posterior estimate can be either "mean" (the default) or the
  "median" of the posterior distribution.

- digits:

  Number of digits.

- flag_signif:

  Add significance flags based on `prob` (default = FALSE).

- ...:

  Additional arguments affecting the summary produced.

## Value

A summary of model parameters.

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
  time = "time",
  tinterval = 1 # interval for approximation of continuous-time dynamic model,
)

# inspect model summary
summary(fit)
# }
```
