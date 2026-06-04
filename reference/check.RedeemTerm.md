# Check Arguments for redeem Model Terms

Internal helper to validate arguments and set defaults.

## Usage

``` r
check.RedeemTerm(
  arglist,
  expected = list(),
  defaults = list(),
  allowed_processes = c("0-1", "1-0"),
  allowed_models = c("dem", "rem"),
  model_type = "dem",
  directed = NULL,
  directed_only = FALSE,
  undirected_only = FALSE
)
```

## Arguments

- arglist:

  List of arguments.

- expected:

  List of expected types/values.

- defaults:

  List of default values.

- allowed_processes:

  Vector of allowed processes ("0-1", "1-0").

- model_type:

  Current model type.
