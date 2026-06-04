// Implementation of helper functions
#include "helper_functions.h"

// [[Rcpp::export]]
arma::uvec topk_indices(const arma::vec& x, arma::uword k) {
  if (k == 0) return arma::uvec();                  // empty
  arma::uvec ord = arma::sort_index(x, "descend");
  return (k >= x.n_elem) ? ord : ord.head(k);
}
// [[Rcpp::export]]
arma::mat do_call_full_cpp(Rcpp::List res, double max_time) {
  if (res.size() == 0) return arma::mat();

  // 1. Calculate total rows
  arma::uword total_rows = 0;
  for (int i = 0; i < res.size(); ++i) {
    if (res[i] != R_NilValue) {
      arma::mat m = Rcpp::as<arma::mat>(res[i]);
      total_rows += m.n_rows;
    }
  }
  if (total_rows == 0) return arma::mat();

  // 2. Concatenate all matrices
  arma::uword n_cols = 0;
  for (int i = 0; i < res.size(); ++i) {
    if (res[i] != R_NilValue) {
      arma::mat m = Rcpp::as<arma::mat>(res[i]);
      if (m.n_rows > 0) {
        n_cols = m.n_cols;
        break;
      }
    }
  }
  
  if (n_cols == 0) return arma::mat();

  arma::mat full_mat(total_rows, n_cols);
  arma::uword current_row = 0;
  for (int i = 0; i < res.size(); ++i) {
    if (res[i] != R_NilValue) {
      arma::mat m = Rcpp::as<arma::mat>(res[i]);
      if (m.n_rows > 0) {
        if (m.n_cols != n_cols) {
          Rcpp::stop("Inconsistent column count in do_call_full_cpp: expected " + std::to_string(n_cols) + ", got " + std::to_string(m.n_cols));
        }
        full_mat.rows(current_row, current_row + m.n_rows - 1) = m;
        current_row += m.n_rows;
      }
    }
  }

  // Column mapping (based on Current_Stat.cpp and wrapper.cpp/rem.cpp/dem.cpp):
  // 0: Time, 1: Pair ID, 2: Status, 3: Event, 4: From, 5: To, 6: Avail_i, 7: Avail_j, 8+: Stats

  // 3. Sort by Pair ID (col 1) then Time (col 0)
  if (full_mat.n_cols < 2) {
      Rcpp::stop("full_mat has too few columns for sorting: " + std::to_string(full_mat.n_cols));
  }
  arma::uvec sort_idx = arma::regspace<arma::uvec>(0, full_mat.n_rows - 1);
  std::sort(sort_idx.begin(), sort_idx.end(), [&](arma::uword a, arma::uword b) {
      if (full_mat(a, 1) != full_mat(b, 1)) return full_mat(a, 1) < full_mat(b, 1);
      return full_mat(a, 0) < full_mat(b, 0);
  });
  full_mat = full_mat.rows(sort_idx);

  // 4. Transform to interval data
  // Output columns: time_end, time, pair_id, status, event, from, to, ...
  arma::mat result(total_rows, n_cols + 1);
  
  for (arma::uword i = 0; i < total_rows; ++i) {
    bool same_pair = (i < total_rows - 1) && (full_mat(i, 1) == full_mat(i + 1, 1));
    double t_start = full_mat(i, 0);
    double t_end = same_pair ? full_mat(i + 1, 0) : max_time;
    
    double event_end = 0.0;
    if (same_pair) {
      event_end = full_mat(i + 1, 3); // Event is Col 3 in full_mat
    }
    
    result(i, 0) = t_end;      // New col 0: time_end
    result(i, 1) = t_start;    // Original col 0: time
    result(i, 2) = full_mat(i, 1); // Pair ID
    result(i, 3) = full_mat(i, 2); // Status
    result(i, 4) = event_end;      // Event
    
    // Copy remaining columns (From, To, Avail, Stats)
    if (n_cols > 4) {
        result.row(i).cols(5, n_cols) = full_mat.row(i).cols(4, n_cols - 1);
    }
  }

  // 5. Remove rows where time == time_end
  arma::uvec to_keep = arma::find(result.col(0) != result.col(1));
  if (to_keep.is_empty()) {
      return arma::mat(0, n_cols + 1);
  }
  return result.rows(to_keep);
}
