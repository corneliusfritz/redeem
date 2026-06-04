# Initialize redeem Model Terms

This is an internal dispatcher that calls the appropriate term
initialization function.

## Usage

``` r
InitRedeemTerm(
  term_name,
  arglist,
  model_type,
  process,
  n_nodes,
  directed = FALSE
)
```

## Arguments

- term_name:

  The name of the term.

- arglist:

  A list of arguments passed to the term in the formula.

- model_type:

  Either "dem" or "rem".

- process:

  Either "0-1" (incidence) or "1-0" (duration).

- n_nodes:

  Number of nodes in the network.

- directed:

  Logical; if TRUE, the model is directed. Used for term-specific logic.
