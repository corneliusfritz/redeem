# Preprocess Model Terms across Multiple Event Streams

Preprocess Model Terms across Multiple Event Streams

## Usage

``` r
preprocess_multi_stream(
  preprocessed,
  n_nodes,
  verbose,
  directed,
  simultaneous_interactions,
  build_time = NULL,
  max_time = -1,
  model_type = "dem",
  impute_zero = TRUE,
  omit_na = TRUE
)
```

## Arguments

- preprocessed:

  Standard output from \`formula_preprocess\`.

- n_nodes:

  Number of nodes.

- verbose:

  Logical; if TRUE, print progress.

- directed:

  Logical; if TRUE, the model is directed.

- simultaneous_interactions:

  Logical; if TRUE, multiple interactions are allowed.

- build_time:

  Numeric; time at which to start building the dataset.

- max_time:

  Numeric; if positive, events after this time are excluded. Defaults to
  `-1.0` (no upper limit).

- model_type:

  Either "dem" or "rem".

- impute_zero:

  Logical; if TRUE, replace NAs in covariates with 0.

- omit_na:

  Logical; if TRUE, call na.omit() on the final table.

## Value

A data.table containing the unified preprocessed data.
