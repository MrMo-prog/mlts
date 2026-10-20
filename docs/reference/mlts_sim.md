# Simulate data from mlts model

Simulate data from mlts model

## Usage

``` r
mlts_sim(
  model,
  default = FALSE,
  N = NULL,
  N_G = NULL,
  TP,
  burn.in = 50,
  seed = NULL,
  seed.true = 1,
  btw.var.sds = NULL,
  exogenous = NULL
)
```

## Arguments

- model:

  `data.frame`. Output of [`mlts_model`](mlts_model.md).

- default:

  logical. If set to `TRUE`, default prior specifications are added.

- N:

  integer Number of observational units.

- N_G:

  vector of integers. Number of observational units per group.

- TP:

  integer. Number of measurements per observational unit.

- burn.in:

  integer. Length of ‘burn-in’ period.

- seed:

  integer. Seed used for data generation.

- seed.true:

  integer. Separate seed used for sampling of true population parameters
  values from plausible ranges for stationary time series.

- btw.var.sds:

  named numeric vector. Provide standard deviation(s) for all exogenous
  between-level variable(s) specified in `model`, e.g.
  (`btw.var.sds = c("covariate1" = 1)`, to set the SD of the variable
  "covariate1" to 1). Mean values of the respective variable(s) will be
  set to 0 per default.

- exogenous:

  Matrix of numeric values of exogenous variables with
  `N`\*(`TP`+`burn.in`) rows and separate columns for each variable.

## Value

An object of class `"mlts_simdata"`. The object is a list containing the
following components:

- model:

  the model object passed to `mlts_sim` with true parameter values used
  in the data generation added in the column `true.val`

- data:

  a long format `data.frame` of the generated time series data

- RE.pars:

  a `matrix` of cluster-specific true values used in the data generation

## Details

A function to generate data from an output of
[`mlts_model`](mlts_model.md).

## Examples

``` r
# \donttest{
# build a simple vector-autoregressive mlts model with two time-series variables
var_model <- mlts_model(q = 2)

# simulate data from this model with default true values
# (true values are randomly drawn from normal distribution)
var_data <- mlts_sim(
  model = var_model,
  N = 50, TP = 30, # number of units and number of measurements per unit
  default = TRUE # use default parameter values
)

# the data set is stored in .$data
head(var_data$data)

# individual parameter values are stored in .$RE.pars
head(var_data$RE.pars)

# if the mltssim-object is used in mlts_fit(), true values
# are added to the fitted object
fit <- mlts_fit(
  model = var_model,
  data = var_data,
  id = "ID", ts = c("Y1", "Y2"), time = "time"
)

# inspect model with true values
head(fit$pop.pars.summary)
# }
```
