# Validate the Structure of a Durational Event List

This function checks the validity of a dyadic interaction matrix by
ensuring that each interaction start event has a corresponding end
event, that no interactions overlap within a dyad, and that no missing
values are present.

## Usage

``` r
check_matrix(df, return_matrix = FALSE, start_time = NULL)
```

## Arguments

- df:

  A data frame with at least four columns, representing events, where:

  Column 1

  :   Event timing or ID (e.g., timestamp).

  Column 2

  :   "From" node ID for the dyadic interaction.

  Column 3

  :   "To" node ID for the dyadic interaction.

  Column 4

  :   Event type (1 for start, 0 for end of interaction).

- return_matrix:

  Logical; if TRUE, returns the (potentially repaired) event matrix.
  Defaults to FALSE, in which case the function returns `TRUE` if the
  matrix is valid.

- start_time:

  Numeric; optional reference time for adding missing start events.
  Defaults to NULL, in which case the earliest time in the data is used.

## Value

Logical; `TRUE` if all interactions are valid, `FALSE` otherwise. If the
data contains missing values, the function issues a warning and returns
`FALSE`.

## Details

The function performs the following checks:

- Missing values: If any are found, a warning is issued and `FALSE` is
  returned.

- Interaction pairing: Each start event (1) must have a corresponding
  end event (0) without overlap.

- Non-overlapping intervals: Ensures that no start event occurs while
  another interaction is active.

## Examples

``` r
# Create a valid event matrix with durational events (start=1, end=0)
df <- matrix(c(
  1.0, 1, 2, 1,
  2.0, 1, 2, 0,
  1.5, 3, 4, 1,
  3.0, 3, 4, 0
), ncol = 4, byrow = TRUE)
colnames(df) <- c("time", "from", "to", "type")

# Check if the event matrix is valid
check_matrix(df)
#> [1] TRUE
```
