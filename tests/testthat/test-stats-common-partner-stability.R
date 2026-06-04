library(testthat)
library(redeem)

test_that("general_common_partner is stable on dissolution", {
  # Small network
  n_nodes <- 5
  
  # Event sequence: 
  # 1. (1,2) starts at time 1
  # 2. (2,3) starts at time 2 -> (1,3) gets 1 shared partner (2)
  # 3. (1,2) ends at time 3   -> (1,3) should STILL have 1 shared partner (2)
  # 4. (1,3) starts at time 4
  
  events <- matrix(c(
    1, 1, 2, 1,
    2, 2, 3, 1,
    3, 1, 2, 0,
    4, 1, 3, 1,
    10, 2, 3, 0,
    10, 1, 3, 0
  ), ncol = 4, byrow = TRUE)
  colnames(events) <- c("time", "from", "to", "type")
  
  f_0_1 <- ~ general_common_partner(transformation = "identity")
  f_1_0 <- ~ Intercept
  
  # Run dem() directly
  fit <- dem(
    events = events,
    formula_0_1 = f_0_1,
    formula_1_0 = f_1_0,
    n_nodes = n_nodes,
    directed = FALSE,
    control = control.redeem(return_data = TRUE)
  )
  
  expect_s3_class(fit, "dem")
  expect_true(all(is.finite(fit$model_0_1$coefficients)))
  
  # Verify statistic stability
  data_0_1 <- fit$model_0_1$data
  dyad_1_3 <- data_0_1[data_0_1$from == 1 & data_0_1$to == 3, ]
  dyad_1_3 <- dyad_1_3[order(dyad_1_3$time), ]
  
  find_stat <- function(df) {
    cols <- colnames(df)
    cols[grep("general_common_partner", cols)][1]
  }
  stat_col <- find_stat(data_0_1)
  
  get_val <- function(df, target, col) {
    row <- df[df$time <= target & df$time_new > target, ]
    if(nrow(row) == 0) return(tail(df[[col]], 1))
    row[[col]][1]
  }

  val_t2 <- get_val(dyad_1_3, 2, stat_col)
  val_t3 <- get_val(dyad_1_3, 3.1, stat_col)
  
  expect_equal(val_t2, 1, info = "Dyad (1,3) should have 1 common partner after (2,3) starts at t=2")
  expect_equal(val_t3, 1, info = "Dyad (1,3) should STILL have 1 common partner after (1,2) ends at t=3 (it is a 'general' stat)")
})
