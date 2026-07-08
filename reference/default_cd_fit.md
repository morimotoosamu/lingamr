# Default causal-discovery fit: joint estimation via lingam_multi_group()

Treats the imputed datasets as "groups" sharing a common causal order.

## Usage

``` r
default_cd_fit(X_list, prior_knowledge, apply_prior_knowledge_softly)
```

## Arguments

- X_list:

  List of imputed datasets (one per repeat)

- prior_knowledge:

  Prior knowledge matrix (NULL allowed)

- apply_prior_knowledge_softly:

  Apply prior knowledge softly (logical)

## Value

`list(causal_order = <integer vector>, adjacency_matrices = <list>)`
