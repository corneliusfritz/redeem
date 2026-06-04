# redeem: Relational and Durational Event Models

[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

**redeem** is an R package for the estimation of **Durational Event Models (DEM)** and **Relational Event Models (REM)**. It features a scalable block-coordinate ascent algorithm designed to handle high-dimensional network data with thousands of actors and time-varying effects.

## Key Features

- **Theory-Aligned**: Implements the framework described in *"Scalable Durational Event Models: Application to Physical and Digital Interactions"* (arXiv:2504.00049).
- **Decoupled intensities**: Separately models the incidence (formation) and duration (dissolution) of events.
- **Scalable Estimation**: Uses a block-coordinate ascent algorithm (MM-based) to efficiently estimate the model.
- **Flexible Baselines**: Supports constant, piecewise-constant, and semiparametric baselines.
- **Support for Constraints**: Handles simultaneous interaction constraints (e.g., phone calls where actors can only be in one call at a time).

## Installation

You can install the development version of **redeem** from GitHub:

```r
# install.packages("devtools")
devtools::install_github("corneliusfritz/redeem")
```

## Quick Start

```r
library(redeem)

# Example: n=50 nodes, directed events
n_nodes <- 50

# events matrix: time, from, to, type (1=start, 0=end)
events <- matrix(c(
  1.0, 1, 5, 1,
  2.5, 1, 5, 0,
  3.2, 2, 10, 1,
  4.8, 2, 10, 0
), ncol = 4, byrow = TRUE)
colnames(events) <- c("time", "from", "to", "type")

# Fit a Durational Event Model
# Modeling both start (0->1) and end (1->0) transitions
fit <- dem(
  events = events,
  n_nodes = n_nodes,
  formula_0_1 = ~ current_interaction() + inertia(),
  formula_1_0 = ~ duration(),
  control = control.redeem(estimate = "Blockwise")
)

# Summarize results
summary(fit)
```

## References

Fritz, C., Rastelli, R., Fop, M., & Caimo, A. (2026). **Scalable Durational Event Models: Application to Physical and Digital Interactions**. *arXiv preprint arXiv:2504.00049*.

