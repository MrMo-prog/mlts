# Fit Bayesian Multilevel Manifest or Latent Time-Series Models

Fit Bayesian Multilevel Manifest or Latent Time-Series Models

## Usage

``` r
mlts_fit(
  model,
  data = NULL,
  id,
  group = NULL,
  ts,
  covariates = NULL,
  outcomes = NULL,
  outcome_pred_btw = NULL,
  center_covs = TRUE,
  time = NULL,
  tinterval = NULL,
  beep = NULL,
  days = NULL,
  n_overnight_NAs,
  max_NA_seq = NULL,
  na.rm = FALSE,
  iter = 500,
  chains = 2,
  cores = 2,
  monitor_person_pars = FALSE,
  monitor_all_pars = FALSE,
  get_SD_latent = FALSE,
  fit_model = TRUE,
  print_message = TRUE,
  print_warning = TRUE,
  ...
)
```

## Arguments

- model:

  `data.frame`. Output of [`mlts_model`](mlts_model.md) and related
  functions.

- data:

  An object of class `data.frame` (or one that can be coerced to that
  class) containing data of all variables used in the model.
  Alternatively, a list object with simulated data created by `mlts_sim`
  can be entered directly and allows for comparison of estimates and
  true population paramter values used in the data generation.

- id:

  Character. The variable in `data` that identifies the observational
  cluster unit. Not necessary when `data` is a list object of simulated
  data generated with `mlts_sim`.

- group:

  Character. The variable in `data` that identifies the grouping
  variable across cluster units.

- ts:

  Character. The variable(s) in `data` that contain the time-series
  construct(s) or their indicator variable(s). If multiple constructs
  are provided in the `model`, multiple entries are necessary. Note that
  the order of variable names provided in `ts` has to match the
  specification made in the `model`. E.g., if multiple constructs (e.g.,
  `mlts_model(q = 2)`) are provided the order of variables names
  provided in `ts` determines which construct is referred to as mu_1,
  phi(1)\_11, etc..

- covariates:

  Named character vector. An optional named vector of characters to
  refer to predictors of random effects as specified in the `model`.
  Note that specifying `covariates` is only necessary if the respective
  variable name(s) in `data` differ from the variables names specified
  in `model`.

- outcomes:

  Named character vector. Similar to `covariates`, an optional named
  vector of characters to refer to outcome predicted by random effects
  as specified in the `model`. Note that specifying `outcomes` is only
  necessary if the respective variable name(s) in `data` differ from the
  outcome variable name(s) specified in `model`.

- outcome_pred_btw:

  Named character vector. Similar to `covariates`, an optional named
  vector of characters to refer to additional between-level variables
  entered as outcome predictor(s) as specified in the `model`. Note that
  specifying `outcome_pred_btw` is only necessary if the respective
  variable name(s) in `data` differ from the variable name(s) specified
  in `model`.

- center_covs:

  Logical. Between-level covariates used as predictors of random effects
  will be grand-mean centered before model fitting by default. Set
  `center_covs` to `FALSE` when including categorical predictors into
  the set of `covariates`. Note that in this case, additional continuous
  covariates should be grand-mean centered prior to using `mlts_fit`. If
  `group` is specified, covariates will be centered on their respective
  group mean.

- time:

  Character. The variable in `data` that contains the (continuous) time
  of observation.

- tinterval:

  The step interval for approximating equally spaced observations in
  time by insertion of missing values, to be specified with respect to
  the time stamp variable provided in time. Procedure for inserting
  missing values resembles the procedure for time shift transformation
  as described in Asparouhov, Hamaker, & Muthén (2018).

- beep:

  Character. The variable in `data` that contains the running beep
  number starting with 1 for each person.

- days:

  Optional. If a running beep identifier is provided via the `beep`
  argument and observations are nested within days (or similar grouping
  unit), the variable in `data` that contains the day identifier can be
  added to correct for overnight lags (see Details).

- n_overnight_NAs:

  Optional. The number of `NA` rows to add after the last observation of
  each day (if `days` is provided).

- max_NA_seq:

  Integer. Specify a maximum number of consecutive missing values. Can
  decrease estimation times drastically in the presence of very long
  sequences of missing values (e.g., when setting tinterval to values of
  small time steps).

- na.rm:

  logical. Per default, missing values remain in the data and will be
  imputed during model estimation. Set to `TRUE` to remove all rows with
  missing values in variables given in `ts`.

- iter:

  A positive integer specifying the number of iterations for each chain
  (including 50% used as warmup). The default is 500.

- chains:

  A positive integer specifying the number of Markov chains. The default
  is 2.

- cores:

  The number of cores to use when executing the Markov chains in
  parallel. The default is 2 (see
  [`stan`](https://mc-stan.org/rstan/reference/stan.html)).

- monitor_person_pars:

  Logical. Should person parameters (i.e., values of the latent
  variables) be stored? Default is FALSE.

- monitor_all_pars:

  Logical. Should all parameters be stored? Default is FALSE.

- get_SD_latent:

  Logical. Set to `TRUE` to obtain standardized estimates in
  multiple-indicator models.

- fit_model:

  Logical. Set to FALSE to avoid fitting the model which may be helpful
  to inspect prepared data used for model estimation (default = TRUE).

- print_message:

  Logical. Print messages based on defined inputs (default = TRUE).

- print_warning:

  Logical. Print warnings based on defined inputs (default = TRUE).

- ...:

  Additional arguments passed to
  [`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

## Value

An object of class `mltsfit`. The object is a list containing the
following components:

- model:

  the model object passed to `mlts_fit`

- data:

  the preprocessed data used for fitting the model

- param.labels:

  a `data.frame` that provides the names of parameters used in the stan
  model. These parameter names are necessary when running standard
  post-processing functions using `mlts_fit$stanfit`

- pop.pars.summary:

  a `data.frame` that contains summary statistics for all parameter in
  `model`

- person.pars.summary:

  if `monitor_person_pars = TRUE`, a `data.frame` containing summary
  statistics for cluster-specific parameters is provided

- standata:

  a `list` with the data as passed to
  [`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)

- stanfit:

  an object of class `stanfit` with the raw output created by
  [`sampling`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html)

- posteriors:

  an `array` of the MCMC chain results for all parameters in `model`
  created by
  [`rstan::extract`](https://mc-stan.org/rstan/reference/stanfit-method-extract.html)
  with `dimnames` adapted to match the parameter names provided in
  `model`

## References

Asparouhov, T., Hamaker, E. L., & Muthén, B. (2018). Dynamic Structural
Equation Models. Structural Equation Modeling: *A Multidisciplinary
Journal*, *25*(3), 359–388.
[doi:10.1080/10705511.2017.1406803](https://doi.org/10.1080/10705511.2017.1406803)

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
  id = "ID", # cluster identifier variable
  time = "time", # time variable
  tinterval = 1 # interval for approximation of equidistant measurements,
)

# inspect model summary
summary(fit)
# }
```
