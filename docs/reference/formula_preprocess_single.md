# Preprocess a single formula for model terms

This function takes an R formula and extracts the necessary information
to build the model's design matrix. It identifies special terms,
transformations, and associated data.

## Usage

``` r
formula_preprocess_single(
  formula,
  n_nodes,
  model_type = "dem",
  process = "0-1",
  directed = FALSE,
  simulation = FALSE
)
```

## Arguments

- formula:

  An R formula object.

- n_nodes:

  Number of nodes.

- model_type:

  Either "dem" or "rem".

- process:

  Either "0-1" (incidence) or "1-0" (duration).

- directed:

  Logical; if TRUE, the model is directed (defaults to FALSE).

- simulation:

  Logical; if TRUE, the formula is being preprocessed for simulation
  (defaults to FALSE).
