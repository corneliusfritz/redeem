# Match user-provided coefficients to internal model names

Match user-provided coefficients to internal model names

## Usage

``` r
match_coefficients(user_coefs, internal_names, internal_keys = NULL)
```

## Arguments

- user_coefs:

  Named or unnamed vector of coefficients provided by the user.

- internal_names:

  Vector of internal coefficient names (labels).

- internal_keys:

  Vector of internal coefficient keys (e.g., intercept,
  dyadic_cov_identity).

## Value

A numeric vector of the same length as internal_names.
