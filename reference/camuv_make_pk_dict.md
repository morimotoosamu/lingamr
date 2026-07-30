# Validate CAM-UV prior knowledge and build the forbidden-cause list

Faithful port of `CAMUV._make_pk_dict()`. Unlike the matrix-based prior
knowledge of the other algorithms in this package, CAM-UV takes a list
of pairs `(i, j)` meaning "variable i cannot be a cause of variable j"
(**1-based** in this port; the Python original is 0-based).

## Usage

``` r
camuv_make_pk_dict(prior_knowledge, d)
```

## Arguments

- prior_knowledge:

  NULL, a 2-column numeric matrix, or a list of length-2 numeric vectors

- d:

  number of variables

## Value

NULL, or a list of length `d`; element `j` holds the variables forbidden
from being a cause of variable `j`
