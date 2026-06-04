# Convert right-hand side of a formula to a list of term information

Convert right-hand side of a formula to a list of term information

## Usage

``` r
rhs_terms_as_list(
  formula,
  n_nodes,
  env = NULL,
  evaluate_calls = FALSE,
  model_type = "dem",
  process = "0-1",
  directed = FALSE
)
```

## Arguments

- formula:

  The formula to parse.

- n_nodes:

  Number of nodes.

- env:

  The environment in which to evaluate terms.

- evaluate_calls:

  Logical; if \`TRUE\`, evaluates the full calls.

- model_type:

  Either "dem" or "rem".

- process:

  Either "0-1" or "1-0".

- directed:

  Logical; if TRUE, the model is directed.

## Value

A list of term information, including labels, base names, and data.
