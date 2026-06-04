# Process event actor columns and automatically identify or validate n_nodes

Process event actor columns and automatically identify or validate
n_nodes

## Usage

``` r
process_event_actors(events, n_nodes = NULL, directed = TRUE)
```

## Arguments

- events:

  A matrix or data frame of events with columns `from` and `to` (or
  columns 2 and 3).

- n_nodes:

  Integer; the total number of actors in the network, or `NULL`.

- directed:

  Logical; whether the interaction events are directed. Defaults to
  TRUE.

## Value

A list containing `events` (potentially modified) and `n_nodes`.
