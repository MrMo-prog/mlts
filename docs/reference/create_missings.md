# Create Missings for Approximation of Continuous Time Dynamic Models

Create Missings for Approximation of Continuous Time Dynamic Models

## Usage

``` r
create_missings(data, tinterval, id, time, btw_vars = NULL)
```

## Arguments

- data:

  An object of class `data.frame` (or one that can be coerced to that
  class) containing data of all variables used in the model.

- tinterval:

  The step interval for approximation for a continuous time DSEM. The
  smaller the step interval, the better the approximation.

- id:

  The variable in `data` that identifies the person or observational
  unit (as character).

- time:

  The variable in `data` that contains the (continuous) time (as
  string).

- btw_vars:

  The names of between-level variables in the data to be added in newly
  created rows with NAs.

## Value

A `data.frame` with missings imputed for use in
[`mlts_fit`](mlts_fit.md).

## Examples

``` r
# create some data for example
data <- data.frame(
  id = rep(c(1, 2), each = 4),
  time = c(0, 3, 4, 6,
           1, 4, 5, 7)
)

# create missings to approximate continuous time process
create_missings(
  data = data, id = "id", time = "time",
  tinterval = 1 # use time interval of 1 minute
)
```
