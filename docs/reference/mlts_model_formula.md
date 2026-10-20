# Create TeX Model Formula from mlts model object

Create TeX Model Formula from mlts model object

## Usage

``` r
mlts_model_formula(
  model,
  file = NULL,
  keep_tex = FALSE,
  ts = NULL,
  covariates = NULL,
  outcomes = NULL
)
```

## Arguments

- model:

  A model built with [`mlts_model`](mlts_model.md).

- file:

  An optional string containing the name of the file and file path. Has
  to end with .pdf file format.

- keep_tex:

  Logical. Should the TeX file be kept (additional to the Rmd file)?
  Defaults to `FALSE`.

- ts:

  To be included in future releases. An optional character vector
  containing the names of the time-series variables or indicators.

- covariates:

  To be included in future releases. An optional character vector
  containing the names of the between-level covariates.

- outcomes:

  To be included in future releases. An optional character vector
  containing the names of the between-level outcomes.

## Value

An RMarkdown file that is automatically rendered to a pdf document.

## Examples

``` r
# \donttest{
# build a simple vector-autoregressive mlts model with two time-series variables
var_model <- mlts_model(q = 2)

# create formula from the specified model
mlts_model_formula(model = var_model)
# }
```
