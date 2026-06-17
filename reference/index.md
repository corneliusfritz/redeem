# Package index

## Model Estimation

Main estimation functions and controls for fitting Durational Event
Models (DEM) and Relational Event Models (REM).

- [`dem()`](https://corneliusfritz.github.io/redeem/reference/dem.md) :
  Durational Event Model (DEM) Estimation
- [`rem()`](https://corneliusfritz.github.io/redeem/reference/rem.md) :
  Relational Event Model (REM) Estimation
- [`dem_object`](https://corneliusfritz.github.io/redeem/reference/dem_object.md)
  [`dem-class`](https://corneliusfritz.github.io/redeem/reference/dem_object.md)
  : The dem Object
- [`rem_object`](https://corneliusfritz.github.io/redeem/reference/rem_object.md)
  [`rem-class`](https://corneliusfritz.github.io/redeem/reference/rem_object.md)
  : The rem Object
- [`control.redeem()`](https://corneliusfritz.github.io/redeem/reference/control.redeem.md)
  : Control Parameters for REDEEM Models
- [`dem.simulate()`](https://corneliusfritz.github.io/redeem/reference/dem.simulate.md)
  : Simulate events based on specified formulas and coefficients
- [`rem.simulate()`](https://corneliusfritz.github.io/redeem/reference/rem.simulate.md)
  : Simulate a Relational Event Model (REM)

## Diagnostics & Post-Estimation

Utilities for residual checking, out-of-sample likelihood evaluation,
baseline intensity plotting, and node ranking.

- [`get_residuals()`](https://corneliusfritz.github.io/redeem/reference/get_residuals.md)
  : Get residuals for model diagnostics (Cox-Snell Residuals)

- [`get_oos_likelihood()`](https://corneliusfritz.github.io/redeem/reference/get_oos_likelihood.md)
  : Out-of-sample Log-Likelihood (Proper Scoring Rule)

- [`get_ranking()`](https://corneliusfritz.github.io/redeem/reference/get_ranking.md)
  : Get ranking for test events (Out-of-Sample Goodness-of-Fit)

- [`plot_baseline()`](https://corneliusfritz.github.io/redeem/reference/plot_baseline.md)
  : Plot the Estimated Baseline Intensity

- [`predict_baseline_trend()`](https://corneliusfritz.github.io/redeem/reference/predict_baseline_trend.md)
  : Predict the baseline intensity trend at one or more future time
  points

- [`summary(`*`<redeem_result>`*`)`](https://corneliusfritz.github.io/redeem/reference/summary.redeem_result.md)
  :

  Summary of a `redeem_result` Model Fit

## Model Terms & Sufficient Statistics

Functions representing standard dyadic, triadic, degree-based, and
exogenous covariate terms to be used inside model formulas.

- [`redeem_terms`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`terms`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`statistics`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`Intercept`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`intercept`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`inertia`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`reciprocity`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_interaction`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`duration`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`number_interaction`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_common_partners`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_common_partner`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_common_partner_OSP`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_common_partner_ISP`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_common_partner_OTP`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_common_partner_ITP`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_common_partners`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_common_partner`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_common_partner_OSP`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_common_partner_ISP`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_common_partner_OTP`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_common_partner_ITP`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_triangle`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_triangle`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`dyadic_cov`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`monadic_cov`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`baseline`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_degree_out_sender`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_degree_out_receiver`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_degree_in_sender`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_degree_in_receiver`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_degree_sum`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_degree_absdiff`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_count_out_sender`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_count_out_receiver`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_count_in_sender`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_count_in_receiver`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_count_sum`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`general_count_absdiff`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_degree_out_sender`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_degree_out_receiver`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_degree_in_sender`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_degree_in_receiver`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_degree_sum`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_degree_absdiff`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_count_out_sender`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_count_out_receiver`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_count_in_sender`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_count_in_receiver`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_count_sum`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`current_count_absdiff`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`degree`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`degrees`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`count`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`triangle`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`common_partner`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`psABBA`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`psABBY`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`psABAY`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`psABXA`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`psABXB`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`psABXY`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  [`ps`](https://corneliusfritz.github.io/redeem/reference/redeem_terms.md)
  : redeem Model Terms

## Data Validation & Helpers

Functions for verifying inputs.

- [`check_matrix()`](https://corneliusfritz.github.io/redeem/reference/check_matrix.md)
  : Validate the Structure of a Durational Event List
