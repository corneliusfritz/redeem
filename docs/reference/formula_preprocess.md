# Preprocess Formulas for Model Terms with Event and Node Information

This function processes two model formulas, each of which can specify
transformation and data arguments, and combines the preprocessed results
with additional event and node information. It identifies unique
coefficient names across both formulas, determining which terms to
include based on uniqueness, and returns structured lists for data,
transformations, and term names.

## Usage

``` r
formula_preprocess(
  formula_0_1 = NULL,
  model_type = "dem",
  formula_1_0 = NULL,
  events = matrix(c(0, 0, 0), nrow = 1),
  n_nodes,
  exo_breaks = NULL,
  directed = FALSE,
  simulation = FALSE
)
```

## Arguments

- formula_0_1:

  Optional; an R formula for the \`0 -\> 1\` terms.

- model_type:

  Either "dem" or "rem".

- formula_1_0:

  Optional; an R formula for the \`1 -\> 0\` terms.

- events:

  A data frame or list representing the events.

- n_nodes:

  An integer specifying the number of nodes.

- exo_breaks:

  Optional; a vector or list specifying external breaks.

- directed:

  Logical; if TRUE, the model is directed (defaults to FALSE).

- simulation:

  Logical; if TRUE, the formula is being preprocessed for simulation
  (defaults to FALSE).

## Value

A list containing the following components:

- events:

  The events input, retained for use in model estimation or evaluation.

- n_nodes:

  The number of nodes specified in the input.

- data_list:

  A combined list of matrices for each term’s data, from both formulas,
  where each matrix corresponds to the data for a specific term.

- transformation_list:

  A combined character vector of transformation types for each term,
  with \`"identity"\` for terms without specified transformations.

- coef_names:

  A character vector of coefficient names for each term, combining terms
  across both formulas and ensuring uniqueness.

- term_names:

  A character vector of term names, ordered to match \`data_list\`.

- preprocess_1_0:

  The output list from \`formula_preprocess_single\` applied to
  \`formula_1_0\`.

- preprocess_0_1:

  The output list from \`formula_preprocess_single\` applied to
  \`formula_0_1\`.

- included_1_0:

  A logical vector indicating whether each term in \`coef_names\` comes
  from \`formula_1_0\`.

- included_0_1:

  A logical vector indicating whether each term in \`coef_names\` comes
  from \`formula_0_1\`.

A list containing the preprocessed information.

## Details

The function first calls \`formula_preprocess_single\` on
\`formula_1_0\` and \`formula_0_1\` separately to obtain individual term
processing details. It then identifies unique terms across both formulas
and combines the term data, transformations, and coefficient names into
a single output list, structured for use in further modeling or
evaluation.

## See also

[`formula_preprocess_single`](https://corneliusfritz.github.io/redeem/reference/formula_preprocess_single.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Define simple event data
event_data <- matrix(c(
  1.2, 1, 5, 1,
  2.5, 1, 5, 0,
  3.1, 2, 8, 1,
  4.4, 2, 8, 0
), ncol = 4, byrow = TRUE)
colnames(event_data) <- c("time", "from", "to", "type")

# Preprocess the formulas
formula_preprocess(
  formula_1_0 = ~ current_interaction() + current_common_partners(),
  formula_0_1 = ~ general_common_partners(),
  events = event_data,
  n_nodes = 10
)
} # }
```
