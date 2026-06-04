# Update core coefficients using R (for correctness assessment)

Update core coefficients using R (for correctness assessment)

## Usage

``` r
update_core_r(
  data,
  covarites,
  prediction,
  est_core,
  identifiable,
  offset_fixed = NULL
)
```

## Arguments

- data:

  Preprocessed data.table.

- prediction:

  Current predicted intensities.

- est_core:

  Initial core coefficients.

- identifiable:

  Logical vector indicating identifiable coefficients.

- offset_fixed:

  Fixed offset (degree + baseline).
