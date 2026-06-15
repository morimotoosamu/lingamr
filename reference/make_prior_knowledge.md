# Create a prior knowledge matrix

Create a prior knowledge matrix

## Usage

``` r
make_prior_knowledge(
  n_variables,
  exogenous_variables = NULL,
  sink_variables = NULL,
  paths = NULL,
  no_paths = NULL,
  labels = NULL
)
```

## Arguments

- n_variables:

  Number of variables

- exogenous_variables:

  Exogenous variables (1-based index or variable name, NULL allowed) The
  specified variables are assumed not to be influenced by any other
  variable

- sink_variables:

  Sink variables (1-based index or variable name, NULL allowed) The
  specified variables are assumed not to influence any other variable

- paths:

  Variable pairs that have a directed path (NULL allowed) Of the form
  list(c(from, to), ...). Specified by index or variable name

- no_paths:

  Variable pairs that have no directed path (NULL allowed) Of the form
  list(c(from, to), ...). Specified by index or variable name

- labels:

  Vector of variable names (NULL allowed) Required when specifying by
  variable name. Pass e.g. colnames() of a data.frame

## Value

Prior knowledge matrix (n_variables x n_variables) -1: unknown, 0: no
path, 1: path exists

## Examples

``` r
# Specify by index
pk <- make_prior_knowledge(6, exogenous_variables = c(4))

# Specify by variable name
pk <- make_prior_knowledge(6,
  exogenous_variables = "x3",
  sink_variables = c("x1", "x4"),
  paths = list(c("x3", "x0"), c("x3", "x2")),
  no_paths = list(c("x5", "x2")),
  labels = c("x0", "x1", "x2", "x3", "x4", "x5")
)
```
