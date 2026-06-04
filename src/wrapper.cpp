#include <set>
#include <algorithm>
#include <string>
#include "dem.h"
#include "rem.h"
#include "current_stat.h"
#include <limits>
#include "helper_functions.h"


// Helper to get matrices at a specific time for initialization
std::vector<arma::mat> update_object_covariates_get_mats(const Rcpp::List& data_list, double current_time, const std::vector<std::string>& transformations) {
  std::vector<arma::mat> current_mats;
  for (R_xlen_t i = 0; i < data_list.size(); ++i) {
    arma::mat raw_mat;
    Rcpp::RObject item = data_list[i];
    if (item == R_NilValue) {
      // Initialize raw_mat to empty matrix for NULL entries
      raw_mat = arma::mat();
    } else if (Rcpp::is<Rcpp::List>(item)) {
      Rcpp::List tv_list = Rcpp::as<Rcpp::List>(item);
      Rcpp::CharacterVector times = tv_list.names();
      int best_idx = -1;
      double max_t = -std::numeric_limits<double>::infinity();
      for (R_xlen_t k = 0; k < times.size(); ++k) {
        double tk = safe_stod(std::string(times[k]), "parsing change points");
        if (tk <= current_time + 1e-10 && tk > max_t) {
          max_t = tk;
          best_idx = k;
        }
      }
      if (best_idx != -1) {
        raw_mat = Rcpp::as<arma::mat>(tv_list[best_idx]);
      } else {
        Rcpp::stop("No measurement time found at or before the current time (" + std::to_string(current_time) + "). Please ensure each time-varying covariate list contains at least one key <= current_time.");
      }
    } else {
      raw_mat = Rcpp::as<arma::mat>(item);
    }
    double K = (raw_mat.n_elem > 0) ? raw_mat.at(0, 0) : 1.0;
    current_mats.push_back(apply_transformation(raw_mat, transformations[i], K));
  }
  return current_mats;
}

// Internal helper to update covariates in a REM object's state
void update_object_covariates(REM& obj, double time) {
  obj.update_covariates(time);
}

void update_object_covariates(DEM& obj, double time) {
  obj.update_covariates(time);
}

// [[Rcpp::export]]
arma::mat preprocess(arma::mat edgelist, std::vector<std::string> terms, Rcpp::List data_list, std::vector<std::string> transformations,
                    unsigned int n_nodes, bool verbose, bool directed, bool simultaneous_interactions,
                    Rcpp::Nullable<Rcpp::NumericVector> window_map = R_NilValue, double build_time = 0.0, double max_time = 0.0) {


  if (edgelist.n_rows == 0) return arma::mat(0, 9 + terms.size());
  arma::mat sorted_edgelist = edgelist.rows(arma::stable_sort_index(edgelist.col(0)));

  for (int i = 0; i < data_list.size(); ++i) {
    if (Rcpp::is<Rcpp::List>(data_list[i])) {
      Rcpp::List tv_list = Rcpp::as<Rcpp::List>(data_list[i]);
      Rcpp::CharacterVector times = tv_list.names();
      double min_t = std::numeric_limits<double>::infinity();
      for (R_xlen_t k = 0; k < times.size(); ++k) {
        double tk = safe_stod(std::string(times[k]), "parsing change points during validation");
        if (tk < min_t) min_t = tk;
      }
      if (min_t > sorted_edgelist(0, 0) + 1e-10) {
        Rcpp::stop("The first measurement of a time-varying covariate (index " + std::to_string(i+1) + ") must be at or before the first event time (" + std::to_string(sorted_edgelist(0, 0)) + ").");
      }
    }
  }

  // 2. Setup the DEM object
  Rcpp::List window_info;
  if (window_map.isNotNull()) {
    Rcpp::NumericVector wm(window_map);
    if (wm.size() > 0 && !Rf_isNull(wm.names())) {
      Rcpp::CharacterVector window_lengths_names = wm.names();
      for (int i = 0; i < wm.size(); ++i) {
        window_info[std::to_string((int)wm[i])] = safe_stod(std::string(window_lengths_names[i]), "parsing window length");
      }
    }
  }

  std::vector<arma::mat> initial_data = update_object_covariates_get_mats(data_list, 0.0, transformations);
  DEM tmp(sorted_edgelist, terms, n_nodes, directed, verbose, data_list, initial_data, transformations, simultaneous_interactions, window_info);

  unsigned int build_idx = 0;
  if (build_time > 0) {
    bool found = false;
    for (unsigned int k = 0; k < tmp.changepoints.size(); ++k) {
      if (tmp.changepoints[k] > build_time + 1e-10) {
        build_idx = k;
        found = true;
        break;
      }
    }
    if (!found) build_idx = tmp.changepoints.size();
  }

  return tmp.preprocess(true, 0, tmp.changepoints.size(), build_idx, max_time);
}

// [[Rcpp::export]]
arma::mat preprocess_rem(arma::mat edgelist, std::vector<std::string> terms, Rcpp::List data_list, std::vector<std::string> transformations,
                         unsigned int n_nodes, bool verbose, bool directed,
                         Rcpp::Nullable<Rcpp::NumericVector> window_map = R_NilValue, double build_time = 0.0, double max_time = 0.0) {

  if (edgelist.n_rows == 0) return arma::mat(0, 9 + terms.size());
  arma::mat sorted_edgelist = edgelist.rows(arma::stable_sort_index(edgelist.col(0)));

  for (int i = 0; i < data_list.size(); ++i) {
    if (Rcpp::is<Rcpp::List>(data_list[i])) {
      Rcpp::List tv_list = Rcpp::as<Rcpp::List>(data_list[i]);
      Rcpp::CharacterVector times = tv_list.names();
      double min_t = std::numeric_limits<double>::infinity();
      for (R_xlen_t k = 0; k < times.size(); ++k) {
        double tk = safe_stod(std::string(times[k]), "parsing change points during validation");
        if (tk < min_t) min_t = tk;
      }
      if (min_t > sorted_edgelist(0, 0) + 1e-10) {
        Rcpp::stop("The first measurement of a time-varying covariate (index " + std::to_string(i+1) + ") must be at or before the first event time (" + std::to_string(sorted_edgelist(0, 0)) + ").");
      }
    }
  }

  // 2. Setup the REM object
  Rcpp::List window_info;
  if (window_map.isNotNull()) {
    Rcpp::NumericVector wm(window_map);
    if (wm.size() > 0 && !Rf_isNull(wm.names())) {
        Rcpp::CharacterVector window_lengths_names = wm.names();
        for (int i = 0; i < wm.size(); ++i) {
            window_info[std::to_string((int)wm[i])] = safe_stod(std::string(window_lengths_names[i]), "parsing window length");
        }
    }
  }

  std::vector<arma::mat> initial_data = update_object_covariates_get_mats(data_list, 0.0, transformations);
  REM tmp(sorted_edgelist, terms, n_nodes, directed, verbose, data_list, initial_data, transformations, window_info);

  unsigned int build_idx = 0;
  if (build_time > 0) {
    bool found = false;
    for (unsigned int k = 0; k < tmp.changepoints.size(); ++k) {
      if (tmp.changepoints[k] > build_time + 1e-10) {
        build_idx = k;
        found = true;
        break;
      }
    }
    if (!found) build_idx = tmp.changepoints.size();
  }

  return tmp.preprocess(true, 0, tmp.changepoints.size(), build_idx, max_time);
}


// [[Rcpp::export]]
double calc_llh_scaled(Rcpp::NumericVector pred,
                       Rcpp::NumericVector intensity,
                                    Rcpp::IntegerVector delta,
                                    Rcpp::IntegerVector pair_id) {
  int n = pred.size();
  if (n == 0) return 0.0;

  double total_llh = 0.0;
  double current_sum_x = 0.0;

  for(int i = 0; i < n; ++i) {
    current_sum_x += pred[i];
    bool is_event = (delta[i] == 1);
    bool is_boundary = (i == n - 1) || (pair_id[i] != pair_id[i + 1]);

    if (is_event || is_boundary) {
      // If the group ends due to a boundary change but no event occurred,
      // the status is 0 (censored). If an event occurred, status is 1.
      int status = delta[i];

      if (status == 1) {
        if (intensity[i] > 0) {
          total_llh += std::log(intensity[i]) - current_sum_x;
        } else {
          return -R_PosInf;
        }
      } else {
        total_llh -= current_sum_x;
      }

      // Reset for the new interval/subject
      current_sum_x = 0.0;
    }
  }

  return total_llh;
}

// [[Rcpp::export]]
arma::uvec find_to(unsigned int to, bool directed, unsigned int n_nodes){
  arma::vec res;
  if(to > n_nodes){
    Rcpp::stop("The values that you want to find have to be below " + std::to_string(n_nodes) );
  }

  if(directed){
    if(to == 1){
      res = (n_nodes-1)*arma::regspace(1,n_nodes-1) ;
    } else if (to == n_nodes) {
      res = (n_nodes-1)*arma::regspace(1,n_nodes-1) -1;
    } else {
      res = arma::join_cols((n_nodes-1)*arma::regspace(0,to-2) + to -2, (n_nodes-1)*arma::regspace(to,n_nodes-1) + to -1);
    }
  } else {
    if(to == 1){
    } else if(to == 2){
      res = to -2;
    } else {
      res = to -2 + arma::cumsum(n_nodes - arma::regspace(1,to-2) - 1) ;
      res = arma::join_cols(arma::regspace(to -2,to -2), res);
    }

  }
  return(arma::conv_to<arma::uvec>::from(res));
}

// [[Rcpp::export]]
arma::uvec find_from(unsigned int from, bool directed, unsigned int n_nodes){
  if(from > n_nodes){
    Rcpp::stop("The values that you want to find have to be below " + std::to_string(n_nodes) );
  }
  if(directed){
    return(arma::conv_to<arma::uvec>::from(arma::regspace(0,n_nodes-2) + (n_nodes-1)*(from-1)));
  } else {
    if(from == n_nodes){
      arma::uvec res;
      return(res);
    } else {
      return(arma::conv_to<arma::uvec>::from(n_nodes*(from-1) - (from)*(from-1)/2 + arma::regspace(1,n_nodes - from) -1));
    }
  }
}

unsigned int find_from_to(unsigned int from, unsigned int to, bool directed, unsigned int n_nodes){
  if(from == to){
    Rcpp::stop("There are no self-loops allowed "
                                   + std::to_string(from) + " to "+
                                     std::to_string(to)+" is thus not possible" );
  }
  if((from > n_nodes) | (to > n_nodes)){
    Rcpp::stop("The values that you want to find have to be below " + std::to_string(n_nodes) );
  }
  if(directed){
    if(from>to){
      return((n_nodes-1)*(from-1) + to-1);
    } else {
      return((n_nodes-1)*(from-1) + to-2);
    }

  } else  {
    if(from>to){
      // If the order is swapped -> swap it back
      int tmp = to;
      to = from;
      from = tmp;
    }
    return(n_nodes*(from-1)- (from-1)*(from)/2 + (to-from) -1);
  }
}

// [[Rcpp::export]]
Rcpp::List get_A_B_C_D_E_F_exact(const arma::vec& from_v, const arma::vec& to_v, 
                                 const arma::vec& weight_v, const arma::mat& covarites,
                                 const arma::vec& time_slices,
                                 unsigned int n_slices, unsigned int n_nodes, bool directed,
                                 bool full_baseline = false) {
  if (to_v.n_elem != from_v.n_elem || time_slices.n_elem != from_v.n_elem) {
    Rcpp::stop("from_v, to_v, and time_slices must have the same length.");
  }
  if (covarites.n_cols > 0 && covarites.n_rows != from_v.n_elem) {
    Rcpp::stop("Covariates matrix must have the same number of rows as the number of events.");
  }
  bool use_weights = (weight_v.n_elem == from_v.n_elem);

  unsigned int n_cov = covarites.n_cols;
  unsigned int n_degree = directed ? 2 * n_nodes : n_nodes;
  unsigned int n_time_params = full_baseline ? n_slices : (n_slices > 0 ? n_slices - 1 : 0);

  arma::mat res_A(n_cov, n_cov, arma::fill::zeros);
  arma::mat res_B(n_degree, n_degree, arma::fill::zeros);
  arma::vec res_C(n_time_params, arma::fill::zeros);
  arma::mat res_D(n_degree, n_cov, arma::fill::zeros);
  arma::mat res_E(n_cov, n_time_params, arma::fill::zeros);
  arma::mat res_F(n_degree, n_time_params, arma::fill::zeros);

  // For each row in the data
  for(unsigned int i = 0; i < from_v.n_elem; i++){
    Rcpp::checkUserInterrupt();
    double weight = use_weights ? weight_v.at(i) : 1.0;

    int from = (int)from_v.at(i) - 1;
    int to = (int)to_v.at(i) - 1;

    if (from < 0 || from >= (int)n_nodes || to < 0 || to >= (int)n_nodes) {
       Rcpp::stop("Node indices must be between 1 and n_nodes. Found indices: " + std::to_string(from+1) + ", " + std::to_string(to+1));
    }

    arma::rowvec row_i;
    if (n_cov > 0) {
        row_i = covarites.row(i);
        res_A += weight * row_i.t() * row_i;
    }

    if (directed) {
      res_B.at(from, from) += weight;
      res_B.at(n_nodes + to, n_nodes + to) += weight;
      res_B.at(from, n_nodes + to) += weight;
      res_B.at(n_nodes + to, from) += weight;

      if (n_cov > 0) {
        res_D.row(from) += weight * row_i;
        res_D.row(n_nodes + to) += weight * row_i;
      }

      int ts = -1;
      if (full_baseline) {
        ts = time_slices.at(i) - 1;
      } else if (time_slices.at(i) != 1) {
        ts = time_slices.at(i) - 2;
      }

      if (ts >= 0) {
        if (ts >= (int)n_time_params) Rcpp::stop("Time slice index out of bounds: " + std::to_string(ts+1));
        res_F.at(from, ts) += weight;
        res_F.at(n_nodes + to, ts) += weight;
      }
    } else {
      res_B.at(from, from) += weight;
      res_B.at(to, to) += weight;
      res_B.at(from, to) += weight;
      res_B.at(to, from) += weight;

      if (n_cov > 0) {
        res_D.row(from) += weight * row_i;
        res_D.row(to) += weight * row_i;
      }

      int ts = -1;
      if (full_baseline) {
        ts = time_slices.at(i) - 1;
      } else if (time_slices.at(i) != 1) {
        ts = time_slices.at(i) - 2;
      }

      if (ts >= 0) {
        if (ts >= (int)n_time_params) Rcpp::stop("Time slice index out of bounds: " + std::to_string(ts+1));
        res_F.at(from, ts) += weight;
        res_F.at(to, ts) += weight;
      }
    }

    int ts = -1;
    if (full_baseline) {
      ts = time_slices.at(i) - 1;
    } else if (time_slices.at(i) != 1) {
      ts = time_slices.at(i) - 2;
    }

    if (ts >= 0) {
      if (ts >= (int)n_time_params) Rcpp::stop("Time slice index out of bounds: " + std::to_string(ts+1));
      res_C.at(ts) += weight;
      if (n_cov > 0) {
        res_E.col(ts) += weight * row_i.t();
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("A_mat") = res_A,
                            Rcpp::Named("B_mat") = res_B,
                            Rcpp::Named("C_mat") = res_C,
                            Rcpp::Named("D_mat") = res_D,
                            Rcpp::Named("E_mat") = res_E,
                            Rcpp::Named("F_mat") = res_F);
}


// [[Rcpp::export]]
Rcpp::List get_A_B_C_exact(const arma::mat& data, int n_nodes, bool directed) {
  if (data.n_rows == 0) {
    unsigned int n_degree = directed ? 2 * n_nodes : n_nodes;
    unsigned int n_cov = data.n_cols > 3 ? data.n_cols - 3 : 0;
    return Rcpp::List::create(Rcpp::Named("A") = arma::mat(n_degree, n_degree, arma::fill::zeros),
                              Rcpp::Named("B") = arma::mat(n_degree, n_cov, arma::fill::zeros),
                              Rcpp::Named("C") = arma::mat(n_cov, n_cov, arma::fill::zeros));
  }
  unsigned int n_degree = directed ? 2 * n_nodes : n_nodes;
  arma::mat res_A(n_degree,n_degree, arma::fill::zeros);
  arma::mat res_B(n_degree,data.n_cols-3, arma::fill::zeros);
  arma::mat res_C(data.n_cols-3,data.n_cols-3, arma::fill::zeros);

  // For each row in the data
  for(unsigned int i = 0; i < data.n_rows; i++){
    Rcpp::checkUserInterrupt();
    double weight = data.at(i, 2);
    arma::rowvec row_i = data.row(i).cols(3, data.n_cols - 1);
    int from = (int)data.at(i,0) - 1;
    int to = (int)data.at(i,1) - 1;
    
    if (from < 0 || from >= n_nodes || to < 0 || to >= n_nodes) {
       Rcpp::stop("Node indices out of bounds in get_A_B_C_exact.");
    }


    if (directed) {
      res_A.at(from, from) += weight;
      res_A.at(n_nodes + to, n_nodes + to) += weight;
      res_A.at(from, n_nodes + to) += weight;
      res_A.at(n_nodes + to, from) += weight;

      res_B.row(from) += weight * row_i;
      res_B.row(n_nodes + to) += weight * row_i;
    } else {
      res_A.at(from, from) += weight;
      res_A.at(to, to) += weight;
      res_A.at(from, to) += weight;
      res_A.at(to, from) += weight;

      res_B.row(from) += weight * row_i;
      res_B.row(to) += weight * row_i;
    }
    res_C += weight * row_i.t() * row_i;
  }
  return Rcpp::List::create(Rcpp::Named("A") = res_A,
                            Rcpp::Named("B") = res_B,
                            Rcpp::Named("C") = res_C);
}

/*
// [[Rcpp::export]]
Rcpp::List get_A_B_C(const arma::mat& data, int n_nodes, bool directed) {
  unsigned int n_degree = directed ? 2 * n_nodes : n_nodes;
  arma::vec res_A(n_degree, arma::fill::zeros);
  arma::mat res_B(n_degree,data.n_cols-3, arma::fill::zeros);
  arma::mat res_C(data.n_cols-3,data.n_cols-3, arma::fill::zeros);
  // For each row in the data
  for(unsigned int i = 0; i < data.n_rows; i++){
    double weight = data.at(i, 2);
    arma::rowvec row_i = data.row(i).cols(3, data.n_cols - 1);
    int from_idx = data.at(i,0) - 1;
    int to_idx = data.at(i,1) - 1;

    if (directed) {
      to_idx += n_nodes;
    }

    res_A.at(from_idx) += weight;
    res_A.at(to_idx) += weight;

    res_B.row(from_idx) += weight * row_i;
    res_B.row(to_idx) += weight * row_i;

    res_C += weight * row_i.t() * row_i;
  }
  return Rcpp::List::create(Rcpp::Named("A") = res_A,
                            Rcpp::Named("B") = res_B,
                            Rcpp::Named("C") = res_C);
}
*/

// [[Rcpp::export]]
arma::mat get_data_degree(Rcpp::NumericVector from_v, Rcpp::NumericVector to_v, Rcpp::NumericVector obs_v, Rcpp::NumericVector pred_v, int n_nodes, bool directed) {
  arma::vec from(from_v.begin(), from_v.size(), false, true);
  arma::vec to(to_v.begin(), to_v.size(), false, true);
  arma::vec obs(obs_v.begin(), obs_v.size(), false, true);
  arma::vec pred(pred_v.begin(), pred_v.size(), false, true);

  arma::mat pairs;
  if (directed) {
    pairs = combinations_directed(n_nodes);
  } else {
    pairs = combinations_undirected(n_nodes);
  }

  arma::vec res_obs = arma::zeros<arma::vec>(pairs.n_rows);
  arma::vec res_pred = arma::zeros<arma::vec>(pairs.n_rows);

  for(unsigned int i = 0; i < from.n_elem; i++){
    unsigned int f = (unsigned int)from.at(i);
    unsigned int t = (unsigned int)to.at(i);
    unsigned int idx;

    idx = find_from_to(f, t, directed, n_nodes);

    res_obs.at(idx) += obs.at(i);
    res_pred.at(idx) += pred.at(i);
  }

  return arma::join_horiz(pairs, res_obs, res_pred);
}

/*
// [[Rcpp::export]]
arma::mat get_data_time(const arma::mat& data, int n_slices) {
  // arma::mat combined(n_slices);
  arma::vec res = arma::regspace(1, n_slices-1);
  arma::vec res_obs = arma::zeros<arma::vec>(n_slices-1);
  arma::vec res_pred = arma::zeros<arma::vec>(n_slices-1);
  for(unsigned int i = 0; i < data.n_rows; i++){
    Rcpp::checkUserInterrupt();
    if(data.at(i,0) != 1){
      res_obs.at(data.at(i,0)-2) += data.at(i,1);
      res_pred.at(data.at(i,0)-2) += data.at(i,2);
    }
  }
  arma::mat res_final = arma::join_horiz(res,res_obs, res_pred);
  return res_final;
}

// Version of the function that is more efficient
// data has the following format (0: time_slice_from, 1:time_slice_to, 2:from, 3:to, 4:pred, 5:obs) for each slice in the data
// [[Rcpp::export]]
arma::mat get_data_time_alt3(const arma::mat& data, int n_slices, arma::vec & changepoints) {
  arma::vec res = arma::regspace(1, n_slices-1);
  arma::vec res_obs = arma::zeros<arma::vec>(n_slices-1);
  arma::vec res_pred = arma::zeros<arma::vec>(n_slices-1);
  double epsilon = 1e-6;

  for(unsigned int i = 0; i < data.n_rows; i++){
    Rcpp::checkUserInterrupt();
    // If the end is also in the first slice, we do nothing as the step function starts at 0 by definition
    // and the respective parameter is set to zero
    if(data.at(i,1) != 1.0){
      // If both observations are not in the same slice
      if(data.at(i,0)!= data.at(i,1)){
        // If the from slice is not the very beginning, the integral starts at the respective changepoint
        if(data.at(i,0) != 1.0){
          // This is the first slice, if the from slice is not the very beginning (if it is in the beginning the effect is zero)
          res_pred.at(data.at(i,0)-2.0) += std::abs(changepoints.at(data.at(i,0)-1.0)-data.at(i,2))*data.at(i,4);
        }
        // Does the information span over one entire slice?
        if (std::abs(data.at(i,1) - data.at(i,0)) > epsilon) {
          // Go over each slice and add the respective information
          for (int j = data.at(i,0) +1;j <= data.at(i,1)-1; ++j) {
            res_pred.at(j-2) += std::abs(changepoints.at(j-1)-changepoints.at(j-2))*data.at(i,4);
          }
        }
        // This adds the information at the end
        res_pred.at(data.at(i,1)-2.0) += std::abs(data.at(i,3)-changepoints.at(data.at(i,1)-2))*data.at(i,4);
      } else {
        // This is what we have to do if the from and to slice are the same
        res_pred.at(data.at(i,1)-2.0) += std::abs(data.at(i,3)-data.at(i,2))*data.at(i,4);
      }
      // If the to slice is the very beginning do not count (as the respective parameter is set to zero)
      res_obs.at(data.at(i,1)-2.0) += data.at(i,5);
    }
  }
  arma::mat res_final = arma::join_horiz(res,res_obs, res_pred);
  return res_final;
}
*/

// [[Rcpp::export]]
arma::mat get_data_time_alt(Rcpp::NumericVector slice_start_v_r, Rcpp::NumericVector slice_end_v_r, Rcpp::NumericVector time_start_v_r, Rcpp::NumericVector time_end_v_r, Rcpp::NumericVector weight_v_r, Rcpp::NumericVector observation_v_r, int n_slices, Rcpp::NumericVector changepoints_r, bool full_baseline = false) {
  // Zero-copy: borrow memory from R vectors
  arma::vec slice_start_v(slice_start_v_r.begin(), slice_start_v_r.size(), false, true);
  arma::vec slice_end_v(slice_end_v_r.begin(), slice_end_v_r.size(), false, true);
  arma::vec time_start_v(time_start_v_r.begin(), time_start_v_r.size(), false, true);
  arma::vec time_end_v(time_end_v_r.begin(), time_end_v_r.size(), false, true);
  arma::vec weight_v(weight_v_r.begin(), weight_v_r.size(), false, true);
  arma::vec observation_v(observation_v_r.begin(), observation_v_r.size(), false, true);
  arma::vec changepoints(changepoints_r.begin(), changepoints_r.size(), false, true);

  if (n_slices <= 0) {
    return arma::mat();
  }

  unsigned int n_time_params = full_baseline ? n_slices : (n_slices > 0 ? n_slices - 1 : 0);
  arma::vec res = arma::regspace(1, n_time_params);
  arma::vec res_obs = arma::zeros<arma::vec>(n_time_params);
  arma::vec res_pred = arma::zeros<arma::vec>(n_time_params);
  double epsilon = 1e-6;

  for (unsigned int i = 0; i < slice_start_v.n_elem; i++) {
    // Retrieve data values once to avoid repeated indexing
    double slice_start = slice_start_v.at(i);
    double slice_end = slice_end_v.at(i);
    double time_start = time_start_v.at(i);
    double time_end = time_end_v.at(i);
    double weight = weight_v.at(i);
    double observation = observation_v.at(i);

    int ts = -1;
    if (full_baseline) {
      ts = slice_end - 1;
    } else if (slice_end != 1.0) {
      ts = slice_end - 2;
    }

    // Skip if not in an estimated slice
    if (ts >= 0) {
      if (slice_start != slice_end) {
        // If the from slice is not the beginning or if we are in full baseline mode,
        // the integral starts at the respective changepoint
        bool slice_start_included = full_baseline || (slice_start != 1.0);
        int ts_start = full_baseline ? (int)slice_start - 1 : (int)slice_start - 2;

        if (slice_start_included) {
          res_pred.at(ts_start) += (changepoints(slice_start - 1.0) - time_start) * weight;
        }

        // If the observation spans over multiple slices, accumulate results
        if (std::abs(slice_end - slice_start) > epsilon) {
          for (int j = slice_start + 1; j <= slice_end - 1; ++j) {
            int ts_mid = full_baseline ? j - 1 : j - 2;
            res_pred.at(ts_mid) += (changepoints.at(j-1)-changepoints.at(j-2))* weight;
          }
        }

        // Add information at the end slice
        double slice_boundary = (slice_end > 1.5) ? changepoints(slice_end - 2.0) : 0.0;
        res_pred.at(ts) += (time_end - slice_boundary) * weight;
      } else {
        // Handle the case where from and to slice are the same
        res_pred.at(ts) += (time_end - time_start) * weight;
      }

      // Update the observation vector
      res_obs.at(ts) += observation;
    }
  }

  // Combine the results into the final matrix
  arma::mat res_final;
  if (n_time_params == 0) {
    res_final = arma::zeros<arma::mat>(0, 3);
  } else {
    res_final = arma::join_horiz(res, res_obs, res_pred);
  }
  return res_final;
}


/*
// [[Rcpp::export]]
arma::mat get_data_time_alt2(const arma::mat& data, int n_slices, arma::vec & changepoints) {
  // arma::mat combined(n_slices);
  arma::vec res = arma::regspace(1, n_slices-1);
  arma::vec res_obs = arma::zeros<arma::vec>(n_slices-1);
  arma::vec res_pred = arma::zeros<arma::vec>(n_slices-1);
  for(unsigned int i = 0; i < data.n_rows; i++){
    Rcpp::checkUserInterrupt();

    // If the end is also in the first slice, we do nothing as the step function starts at 0 by definition
    // and the respective parameter is set to zero
    // Add the information at the very beginning and very end to the respecitve slices
    if(data.at(i,1) != 1.0){
      // If both observations are in the same slice, we only need to take care of the respective parameter
      if(data.at(i,0)!= data.at(i,1)){

        if((data.at(i,0)) !=2.0){

          // if(res_pred.at(0) > 0 ){
          //   break;
          // }
          res_pred.at(data.at(i,0)-2.0) += (changepoints.at(data.at(i,0)-1)-data.at(i,2))*data.at(i,4);
        }
        if((data.at(i,1))!=2.0){
          res_pred.at(data.at(i,1)-2.0) += (data.at(i,3)-changepoints.at(data.at(i,1)-2))*data.at(i,4);
        }

        // Does the information span over one entire slice?
        if(data.at(i,1)-data.at(i,0) > 1){
          // Go over each slice and add the respective information
          for (int j = data.at(i,0)+1;j <= data.at(i,1)-1; ++j) {

            if(j !=2){
              res_pred.at(j-2) += (changepoints.at(j-1)-changepoints.at(j-2))*data.at(i,4);
            }

          }
        }
      } else {
        if(data.at(i,1) !=2.0){
          res_pred.at(data.at(i,1)-2.0) += (data.at(i,3)-data.at(i,2))*data.at(i,4);
        }

      }
      // If the to slice is the very beginning do not count (as the respective parameter is set to zero)
      res_obs.at(data.at(i,1)-2.0) += data.at(i,5);
    }
  }
  arma::mat res_final = arma::join_horiz(res,res_obs, res_pred);
  return res_final;
}
*/
// [[Rcpp::export]]
double integrate_step_function(const arma::vec& change_points, const arma::vec& values, double from, double to, double from_slice, double to_slice) {
  double integral = 0.0;
  // Integrate between the relevant intervals
  // If the slices are the same
  int s_from = (int)from_slice - 1;
  int s_to = (int)to_slice - 1;

  if (s_to != s_from) {
    integral += (change_points[s_from] - from) * values[s_from];
    if (s_to - s_from > 1) {
      for (int i = s_from + 1; i < s_to; ++i) {
        integral += (change_points[i] - change_points[i - 1]) * values[i];
      }
    }
    integral += (to - change_points[s_to - 1]) * values[s_to];
  } else {
    integral += (to - from) * values[s_from];
  }
  return integral;
}


// [[Rcpp::export]]
arma::vec get_time_offset(Rcpp::NumericVector from_slice_r, Rcpp::NumericVector to_slice_r, Rcpp::NumericVector from_time_r, Rcpp::NumericVector to_time_r, Rcpp::NumericVector est_time_r, Rcpp::NumericVector changepoints_r) {
  // Zero-copy: borrow memory from R vectors
  arma::vec from_slice(from_slice_r.begin(), from_slice_r.size(), false, true);
  arma::vec to_slice(to_slice_r.begin(), to_slice_r.size(), false, true);
  arma::vec from_time(from_time_r.begin(), from_time_r.size(), false, true);
  arma::vec to_time(to_time_r.begin(), to_time_r.size(), false, true);
  arma::vec est_time(est_time_r.begin(), est_time_r.size(), false, true);
  arma::vec changepoints(changepoints_r.begin(), changepoints_r.size(), false, true);

  arma::vec res = arma::vec(from_slice.n_elem);
  for(unsigned int i = 0; i < from_slice.n_elem; i++){
    Rcpp::checkUserInterrupt();
    res.at(i) = integrate_step_function(changepoints,est_time, from_time.at(i), to_time.at(i), from_slice.at(i), to_slice.at(i));
  }
  return res;
}



// Memory-efficient MM update for degree effects
// Returns updated degrees in O(M) time where M is number of intervals
// [[Rcpp::export]]
Rcpp::List update_degree_fast(const arma::vec& from_v, const arma::vec& to_v, 
                             const arma::vec& event_v, const arma::vec& prediction_v,
                             const arma::vec& weights,
                             arma::vec est_degree, unsigned int n_nodes, 
                             bool directed, bool update_sender = true) {
  
  if (from_v.n_elem != to_v.n_elem || from_v.n_elem != event_v.n_elem || from_v.n_elem != prediction_v.n_elem) {
    Rcpp::stop("from_v, to_v, event_v, and prediction_v must have the same length.");
  }
  if (weights.n_elem > 0 && weights.n_elem != from_v.n_elem) {
    Rcpp::stop("weights vector must either be empty or match the length of other inputs.");
  }
  unsigned int required_size = directed ? 2 * n_nodes : n_nodes;
  if (est_degree.n_elem != required_size) {
    Rcpp::stop("est_degree has incorrect size. Expected " + std::to_string(required_size) + " but got " + std::to_string(est_degree.n_elem));
  }

  arma::vec obs_sum(n_nodes, arma::fill::zeros);
  arma::vec pred_sum(n_nodes, arma::fill::zeros);
  
  bool use_weights = (weights.n_elem == from_v.n_elem);

  if (directed) {
    if (update_sender) {
      // Update sender effects (mu)
      for (unsigned int i = 0; i < from_v.n_elem; ++i) {
        unsigned int u = (unsigned int)from_v.at(i) - 1;
        if (u >= n_nodes) Rcpp::stop("Sender index out of bounds in update_degree_fast.");
        double w = use_weights ? weights.at(i) : 1.0;
        obs_sum(u) += w * event_v.at(i);
        pred_sum(u) += w * prediction_v.at(i);
      }
      for (unsigned int i = 0; i < n_nodes; ++i) {
        double o = obs_sum(i);
        double p = pred_sum(i);
        if (o < 1e-15) o = 1e-15;
        if (p < 1e-15) p = 1e-15;
        est_degree(i) += std::log(o / p);
      }

    } else {
      // Update receiver effects (nu)
      for (unsigned int i = 0; i < to_v.n_elem; ++i) {
        unsigned int v = (unsigned int)to_v.at(i) - 1;
        if (v >= n_nodes) Rcpp::stop("Receiver index out of bounds in update_degree_fast.");
        double w = use_weights ? weights.at(i) : 1.0;
        obs_sum(v) += w * event_v.at(i);
        pred_sum(v) += w * prediction_v.at(i);
      }
      for (unsigned int i = 0; i < n_nodes; ++i) {
        double o = obs_sum(i);
        double p = pred_sum(i);
        if (o < 1e-15) o = 1e-15;
        if (p < 1e-15) p = 1e-15;
        est_degree(n_nodes + i) += std::log(o / p);
      }

    }
  } else {
    // Undirected: each interval counts for both sender and receiver
    for (unsigned int i = 0; i < from_v.n_elem; ++i) {
      unsigned int u = (unsigned int)from_v.at(i) - 1;
      unsigned int v = (unsigned int)to_v.at(i) - 1;
      if (u >= n_nodes || v >= n_nodes) Rcpp::stop("Node index out of bounds in update_degree_fast (undirected).");
      double w = use_weights ? weights.at(i) : 1.0;
      
      obs_sum(u) += w * event_v.at(i);
      pred_sum(u) += w * prediction_v.at(i);
      
      obs_sum(v) += w * event_v.at(i);
      pred_sum(v) += w * prediction_v.at(i);
    }
    for (unsigned int i = 0; i < n_nodes; ++i) {
      // Matching the original MM implementation's epsilon flooring:
      double o = obs_sum(i);
      double p = pred_sum(i);
      if (o < 1e-15) o = 1e-15;
      if (p < 1e-15) p = 1e-15;
      
      double current_exp_mu = std::exp(est_degree(i));
      double pred_sum_other = p / current_exp_mu;
      est_degree(i) = std::log(std::sqrt(current_exp_mu * o / pred_sum_other));
    }
  }

  
  return Rcpp::List::create(Rcpp::Named("est_degree") = est_degree);
}

// [[Rcpp::export]]
arma::vec update_degree_directed(const arma::mat& update_data, const arma::vec& est_mu, const arma::vec& est_nu, int n_nodes) {
  arma::vec res_mu = arma::zeros<arma::vec>(n_nodes);
  arma::vec res_nu = arma::zeros<arma::vec>(n_nodes);
  arma::vec obs_sum_s = arma::zeros<arma::vec>(n_nodes);
  arma::vec pred_sum_s = arma::zeros<arma::vec>(n_nodes);
  arma::vec obs_sum_r = arma::zeros<arma::vec>(n_nodes);
  arma::vec pred_sum_r = arma::zeros<arma::vec>(n_nodes);

  for(unsigned int k = 0; k < update_data.n_rows; ++k) {
    unsigned int u_idx = (unsigned int)update_data.at(k, 0) - 1;
    unsigned int v_idx = (unsigned int)update_data.at(k, 1) - 1;
    double obs = update_data.at(k, 2);
    double offset_exp = update_data.at(k, 3);

    obs_sum_s(u_idx) += obs;
    pred_sum_s(u_idx) += offset_exp * std::exp(est_nu(v_idx));

    obs_sum_r(v_idx) += obs;
    pred_sum_r(v_idx) += offset_exp * std::exp(est_mu(u_idx));
  }

  for(unsigned int i = 0; i < (unsigned int)n_nodes; ++i) {
    if (pred_sum_s(i) <= 1e-15) pred_sum_s(i) = 1e-15;
    if (pred_sum_r(i) <= 1e-15) pred_sum_r(i) = 1e-15;

    res_mu(i) = std::log(obs_sum_s(i) / pred_sum_s(i));
    res_nu(i) = std::log(obs_sum_r(i) / pred_sum_r(i));
  }

  return join_vert(res_mu, res_nu);
}

// [[Rcpp::export]]
arma::vec update_degree(const arma::mat& update_data, arma::vec& est_degree, int n_nodes, bool directed) {
  if (directed) {
    arma::vec est_mu = est_degree.head(n_nodes);
    arma::vec est_nu = est_degree.tail(n_nodes);
    return update_degree_directed(update_data, est_mu, est_nu, n_nodes);
  }

  arma::vec res = arma::zeros<arma::vec>(n_nodes);
  arma::vec obs_sum = arma::zeros<arma::vec>(n_nodes);
  arma::vec pred_sum = arma::zeros<arma::vec>(n_nodes);

  for(unsigned int k = 0; k < update_data.n_rows; ++k) {
    unsigned int u_idx = (unsigned int)update_data.at(k, 0) - 1;
    unsigned int v_idx = (unsigned int)update_data.at(k, 1) - 1;
    double obs = update_data.at(k, 2);
    double offset_exp = update_data.at(k, 3);

    obs_sum(u_idx) += obs;
    obs_sum(v_idx) += obs;
    pred_sum(u_idx) += offset_exp * std::exp(est_degree(v_idx));
    pred_sum(v_idx) += offset_exp * std::exp(est_degree(u_idx));
  }

  for(unsigned int i = 0; i < (unsigned int)n_nodes; ++i) {
    if (obs_sum(i) <= 1e-15) obs_sum(i) = 1e-15;
    if (pred_sum(i) <= 1e-15) pred_sum(i) = 1e-15;
    res(i) = std::log(std::sqrt(std::exp(est_degree(i)) * obs_sum(i) / (pred_sum(i))));
  }

  return res;
}




// [[Rcpp::export]]
double eval_llh_pois(const arma::vec& outcome, const arma::vec& mean, const arma::vec& weights){
  if (outcome.n_elem == 0) return 0.0;
  bool use_weights = (weights.n_elem == outcome.n_elem);
  
  arma::uvec pos = arma::find(outcome > 0.0);
  if (pos.n_elem > 0) {
      if (arma::any(mean.elem(pos) <= 0.0)) return -R_PosInf;
      double log_term = use_weights ? 
          arma::dot(weights.elem(pos) % outcome.elem(pos), arma::log(mean.elem(pos))) :
          arma::dot(outcome.elem(pos), arma::log(mean.elem(pos)));
      double mean_term = use_weights ?
          arma::dot(weights, mean) :
          arma::accu(mean);
      return log_term - mean_term;
  } else {
      return use_weights ? -arma::dot(weights, mean) : -arma::accu(mean);
  }
}

// [[Rcpp::export]]
double eval_llh_pois_log(const arma::vec& outcome, const arma::vec& log_mean, const arma::vec& weights){
  if (outcome.n_elem == 0) return 0.0;
  bool use_weights = (weights.n_elem == outcome.n_elem);
  
  arma::uvec pos = arma::find(outcome > 0.0);
  if (pos.n_elem > 0) {
      // Check for -Inf in log_mean where outcome > 0
      if (arma::any(log_mean.elem(pos) == -arma::datum::inf)) return -R_PosInf;
      
      double log_term = use_weights ?
          arma::dot(weights.elem(pos) % outcome.elem(pos), log_mean.elem(pos)) :
          arma::dot(outcome.elem(pos), log_mean.elem(pos));
      double exp_term = use_weights ?
          arma::dot(weights, arma::exp(log_mean)) :
          arma::accu(arma::exp(log_mean));
      return log_term - exp_term;
  } else {
      return use_weights ? -arma::dot(weights, arma::exp(log_mean)) : -arma::accu(arma::exp(log_mean));
  }
}


// Gradient computation
// [[Rcpp::export]]
arma::vec gd_step_halfing(const arma::mat& X,
                                 const arma::vec& y,
                                 const arma::vec& beta,
                                 const arma::vec&  offset,
                                 double c) {

  arma::vec mu = exp(X * beta + offset);
  arma::vec gradient = X.t() * (y - mu);
  arma::vec X_beta_offset = X * beta + offset;
  arma::vec X_grad = X * gradient;

  double llh_old = eval_llh_pois(y, exp(X_beta_offset), arma::vec());
  double alpha, llh_new;
  double grad_norm_sq = arma::dot(gradient, gradient);

  for(int i = 0; i < 1000; i++){
    alpha = 2.0 / std::pow(2.0, i);
    // Reuse X_beta_offset and X_grad to compute new linear predictor
    // linear_predictor_new = X * (beta + alpha * gradient) + offset
    //                      = (X * beta + offset) + alpha * (X * gradient)
    llh_new = eval_llh_pois(y, exp(X_beta_offset + alpha * X_grad), arma::vec());

    double threshold = llh_old + c * alpha * grad_norm_sq * alpha;
    if(llh_new > threshold){
      return beta + alpha * gradient;
    }
  }
  return beta;
}

arma::vec gd_stepsize(const arma::mat& X,
                          const arma::vec& y,
                          const arma::vec& beta,
                          const arma::vec&  offset,
                           double stepsize) {
  arma::vec mu = exp(X * beta + offset);
  arma::vec gradient = X.t() * (y - mu)/y.size();
  return(beta + stepsize*gradient);
}


// [[Rcpp::export]]
Rcpp::List gd_estimation(const arma::mat& X,
                          const arma::vec& y,
                          const arma::vec& coef,
                          const arma::vec&  offset,
                          const int & max_iter = 1000,
                          double stepsize = 10.0,
                          const double & tol = 0.001,
                          double c = 0.001) {
  arma::vec beta_new = coef;
  arma::vec beta_old = coef;
  arma::vec llh(max_iter);
  arma::mat beta_hist(max_iter, coef.n_elem);
  int k = 0;
  for(int i = 0; i < max_iter; i++){
    Rcpp::checkUserInterrupt();
    if(stepsize == 10.0){
      beta_new =  gd_step_halfing(X,y,beta_new,offset, c);
    } else {
      beta_new =  gd_stepsize(X,y,beta_new,offset,stepsize);
    }
    llh.at(i) = eval_llh_pois(y, exp(X * beta_new + offset), arma::vec());
    beta_hist.row(i) = beta_new.t();
    if(arma::norm(beta_new - beta_old, 2) < tol){
      k++; // increment k for the final successful iteration
      break;
    } else {
      beta_old = beta_new;
      k ++;
    }
  }
  return Rcpp::List::create(Rcpp::Named("coef") = beta_new,
                            Rcpp::Named("llh") = llh.head(k),
                            Rcpp::Named("coef_hist") = beta_hist.rows(0, k-1));
}


// [[Rcpp::export]]
Rcpp::List rem_simulate_from_empty_timevarying(unsigned int n_nodes,
                                              arma::vec coef,
                                              arma::vec degree_coef,
                                              arma::vec baseline,
                                              unsigned int n_events,
                                              unsigned int max_events,
                                              double time,
                                              arma::vec time_changepoints,
                                              std::vector<std::string> transformations,
                                              std::vector<std::string> terms,
                                              bool verbose,
                                              bool directed,
                                              Rcpp::List data_list,
                                              unsigned int seed,
                                              arma::vec block,
                                              Rcpp::List window_info) {

  // 1. Identify all change points
  std::set<double> change_points;
  change_points.insert(0.0);
  for (unsigned int i = 0; i < time_changepoints.n_elem; ++i) {
      if (time_changepoints[i] < (time > 0 ? time : R_PosInf)) change_points.insert(time_changepoints[i]);
  }
  // Go through data_list to find any additional change points from time-varying covariates
  for (int i = 0; i < data_list.size(); ++i) {
    if (Rcpp::is<Rcpp::List>(data_list[i])) {
      Rcpp::List tv_list = Rcpp::as<Rcpp::List>(data_list[i]);
      Rcpp::CharacterVector times = tv_list.names();
      for (R_xlen_t k = 0; k < times.size(); ++k) {
        double tk = safe_stod(std::string(times[k]), "parsing change points");
        if (tk > 1e-10 && tk < (time > 0 ? time : R_PosInf)) change_points.insert(tk);
      }
    }
  }
  if (time > 0) change_points.insert(time);

  // Ensure baseline and time_changepoints have 0 at start for indexed access
  arma::vec extended_cps = time_changepoints;
  if (extended_cps.n_elem == 0 || extended_cps[0] > 1e-10) {
      extended_cps = arma::join_vert(arma::vec({0.0}), time_changepoints);
  }
  arma::vec extended_baseline = baseline;
  if (extended_baseline.n_elem < extended_cps.n_elem) {
      extended_baseline = arma::join_vert(arma::vec({0.0}), baseline);
  }

  // 2. Setup the REM object
  arma::mat edgelist(0, 4);
  std::vector<arma::mat> initial_data = update_object_covariates_get_mats(data_list, 0.0, transformations);
  REM tmp(edgelist, terms, n_nodes, directed, verbose, data_list, initial_data, transformations, block.n_elem > 1 ? block : arma::vec(), window_info);

  Rcpp::List res_list(extended_cps.n_elem);
  for (int i = 0; i < res_list.size(); ++i) res_list[i] = arma::mat(0, 3);

  unsigned int events_simulated = 0;
  auto cp_it = change_points.begin();
  if (seed > 0) {
      Rcpp::Environment base("package:base");
      Rcpp::Function set_seed = base["set.seed"];
      set_seed(seed);
  }
  unsigned int current_baseline_idx = 0;

    unsigned int target_events = (max_events > 0) ? max_events : (n_events > 0 ? n_events : 400000);
  if (n_events > 0) target_events = n_events;

  while (events_simulated < target_events) {
    double cp_start = *cp_it;
    auto next_cp_it = std::next(cp_it);

    double cp_end;
    if (next_cp_it == change_points.end()) {
        if (n_events > 0 && time == 0) {
            cp_end = 0; // Infinite for sample
        } else {
            break;
        }
    } else {
        cp_end = *next_cp_it;
    }

    update_object_covariates(tmp, cp_start);

    double current_baseline = 0;
    if (extended_baseline.n_elem > 0) {
        for (int i = extended_cps.n_elem - 1; i >= 0; --i) {
            if (extended_cps[i] <= cp_start + 1e-10) {
                unsigned int b_idx = std::min((unsigned int)i, (unsigned int)extended_baseline.n_elem - 1);
                current_baseline = extended_baseline[b_idx];
                current_baseline_idx = i;
                break;
            }
        }
    }

    arma::vec current_coef = coef;
    if (current_coef.n_elem == 0) {
        current_coef = arma::vec({current_baseline});
    } else {
        current_coef.at(0) += current_baseline;
    }

    double interval_duration = (cp_end > 0) ? (cp_end - cp_start) : 0;
    unsigned int remaining = target_events - events_simulated;
    arma::mat interval_events = tmp.sample(n_events > 0 ? remaining : 0, interval_duration, current_coef, degree_coef, max_events);

    if (interval_events.n_rows > 0) {
      // tmp.sample already returns absolute times since tmp.time_old is updated. No shifting needed.


      arma::mat existing = Rcpp::as<arma::mat>(res_list[current_baseline_idx]);
      res_list[current_baseline_idx] = arma::join_vert(existing, interval_events);
      events_simulated += interval_events.n_rows;
    }

    cp_it = next_cp_it;
  }

  // The final status matrix is appended as the LAST element, so the list total size is n_intervals + 2
  Rcpp::List final_res(res_list.size() + 1);
  for(int i=0; i < res_list.size(); ++i) final_res[i] = res_list[i];
  final_res[res_list.size()] = tmp.data_dem.current_stats.return_data();

  return final_res;
}

void print_terms_newline(std::vector<std::string> terms) {
  for (const std::string& term : terms) {
    Rcpp::Rcout << term << "\n";
  }
}


// [[Rcpp::export]]
Rcpp::List simulate_from_empty_timevarying(std::vector<std::string> terms,
                                           Rcpp::List data_list,
                               std::vector<std::string> transformations,
                               unsigned int n_nodes, bool verbose, bool directed,
                               arma::vec coef_0_1, arma::vec coef_1_0,
                               arma::vec degree_coef_0_1, arma::vec degree_coef_1_0,
                               arma::vec time_changepoints,
                               arma::vec baseline_0_1,
                               arma::vec baseline_1_0,
                               bool simultaneous_interactions,
                                unsigned int n_events,
                               unsigned int max_events,
                               double time,
                               Rcpp::List window_info,
                               unsigned int seed = 0) {
  if (seed > 0) {
      Rcpp::Environment base("package:base");
      Rcpp::Function set_seed = base["set.seed"];
      set_seed(seed);
  }

  // Identify all change points of the baseline
  std::set<double> change_points;
  change_points.insert(0.0);
  for (unsigned int i = 0; i < time_changepoints.n_elem; ++i) {
      if (time_changepoints[i] < (time > 0 ? time : R_PosInf)) change_points.insert(time_changepoints[i]);
  }
  // Read in the timevarying covariates and extract change points from there as well
  for (int i = 0; i < data_list.size(); ++i) {
    if (Rcpp::is<Rcpp::List>(data_list[i])) {
      Rcpp::List tv_list = Rcpp::as<Rcpp::List>(data_list[i]);
      Rcpp::CharacterVector times = tv_list.names();
      for (R_xlen_t k = 0; k < times.size(); ++k) {
        double tk = safe_stod(std::string(times[k]), "parsing change points");
        if (tk > 1e-10 && tk < (time > 0 ? time : R_PosInf)) change_points.insert(tk);
      }
    }
  }
  if (time > 0) change_points.insert(time);

  // 2. Setup the DEM object
  arma::mat edgelist(0, 4);
  std::vector<arma::mat> initial_data = update_object_covariates_get_mats(data_list, 0.0, transformations);
  DEM tmp(edgelist, terms, n_nodes, directed, verbose, data_list, initial_data, transformations, simultaneous_interactions, window_info);

  Rcpp::List res_list(time == 0 ? change_points.size() : change_points.size() - 1);
  for (int i = 0; i < res_list.size(); ++i) res_list[i] = arma::mat(0, 4);

  unsigned int events_simulated = 0;

  arma::vec extended_cps = time_changepoints;
  if (extended_cps.n_elem == 0 || extended_cps[0] > 1e-10) {
      extended_cps = arma::join_vert(arma::vec({0.0}), time_changepoints);
  }
  arma::vec extended_baseline_0_1 = baseline_0_1;
  if (extended_baseline_0_1.n_elem < extended_cps.n_elem) {
      extended_baseline_0_1 = arma::join_vert(arma::vec({0.0}), baseline_0_1);
  }
  arma::vec extended_baseline_1_0 = baseline_1_0;
  if (extended_baseline_1_0.n_elem < extended_cps.n_elem) {
      extended_baseline_1_0 = arma::join_vert(arma::vec({0.0}), baseline_1_0);
  }

  auto cp_it = change_points.begin();
  unsigned int current_baseline_idx = 0;

    unsigned int target_events = (max_events > 0) ? max_events : (n_events > 0 ? n_events : 400000);
  if (n_events > 0) target_events = n_events;

  while (events_simulated < target_events) {
    double cp_start = *cp_it;
    auto next_cp_it = std::next(cp_it);

    double cp_end;
    if (next_cp_it == change_points.end()) {
        if (n_events > 0 && time == 0) {
            cp_end = 0; // Infinite for sample
        } else {
            break;
        }
    } else {
        cp_end = *next_cp_it;
    }

    // 1. Refresh state for this sub-interval
    update_object_covariates(tmp, cp_start);

    // 2. Find correct coefficients and baseline for this sub-interval
    double current_b_0_1 = 0;
    double current_b_1_0 = 0;
    if (extended_cps.n_elem > 0) {
        for (int i = extended_cps.n_elem - 1; i >= 0; --i) {
            if (extended_cps[i] <= cp_start + 1e-10) {
                unsigned int idx_0_1 = std::min((unsigned int)i, (unsigned int)extended_baseline_0_1.n_elem - 1);
                unsigned int idx_1_0 = std::min((unsigned int)i, (unsigned int)extended_baseline_1_0.n_elem - 1);
                current_b_0_1 = extended_baseline_0_1[idx_0_1];
                current_b_1_0 = extended_baseline_1_0[idx_1_0];
                current_baseline_idx = i;
                break;
            }
        }
    }

    arma::vec current_coef_0_1 = coef_0_1;
    arma::vec current_coef_1_0 = coef_1_0;
    if (current_coef_0_1.n_elem == 0) {
        current_coef_0_1 = arma::vec({current_b_0_1});
    } else {
        current_coef_0_1.at(0) += current_b_0_1;
    }
    if (current_coef_1_0.n_elem == 0) {
        current_coef_1_0 = arma::vec({current_b_1_0});
    } else {
        current_coef_1_0.at(0) += current_b_1_0;
    }
    // Rcpp::Rcout << "Simulating for interval [" << cp_start << ", " << cp_end << ") with baseline_0_1 = " << current_coef_0_1.t()
    //             << " and baseline_1_0 = " << current_coef_1_0.t() << std::endl;
    // print_terms_newline(terms);
    // 3. Sample events for this sub-interval
    double interval_duration = (cp_end > 0) ? (cp_end - cp_start) : 0;
    unsigned int remaining = target_events - events_simulated;
    arma::mat interval_events = tmp.sample(n_events > 0 ? remaining : 0, interval_duration, current_coef_0_1, current_coef_1_0, degree_coef_0_1, degree_coef_1_0, max_events);

    if (interval_events.n_rows > 0) {
      // tmp.sample already returns absolute times since tmp.time_old is updated. No shifting needed.
      arma::mat existing = Rcpp::as<arma::mat>(res_list[current_baseline_idx]);
      res_list[current_baseline_idx] = arma::join_vert(existing, interval_events);
      events_simulated += interval_events.n_rows;
    }

    cp_it = next_cp_it;
  }

  // The final status matrix is appended as the LAST element, so the list total size is n_intervals + 2
  Rcpp::List final_res(res_list.size() + 1);
  for(int i=0; i < res_list.size(); ++i) final_res[i] = res_list[i];
  final_res[res_list.size()] = tmp.data_dem.current_stats.return_data();

  return final_res;
}

// [[Rcpp::export]]
arma::mat rem_simulate_from_empty(std::vector<std::string> terms, Rcpp::List data_list,
                              std::vector<std::string> transformations,
                              unsigned int n_nodes, bool verbose, bool directed,
                              arma::vec coef,
                              arma::vec degree_coef,
                              unsigned int n_events, double time,
                              unsigned int max_events,
                              arma::vec block,
                              Rcpp::List window_info,
                              unsigned int seed = 0) {

  if (seed > 0) {
      Rcpp::Environment base("package:base");
      Rcpp::Function set_seed = base["set.seed"];
      set_seed(seed);
  }
  // 1. Identify all change points from time-varying covariates
  std::set<double> change_points;
  change_points.insert(0.0);
  bool is_time_varying = false;
  for (int i = 0; i < data_list.size(); ++i) {
    if (Rcpp::is<Rcpp::List>(data_list[i])) {
      is_time_varying = true;
      Rcpp::List tv_list = Rcpp::as<Rcpp::List>(data_list[i]);
      Rcpp::CharacterVector times = tv_list.names();
      for (R_xlen_t k = 0; k < times.size(); ++k) {
        double tk = safe_stod(std::string(times[k]), "parsing change points");
        if (tk > 1e-10) change_points.insert(tk);
      }
    }
  }

  if (!is_time_varying) {
    std::vector<arma::mat> static_data;
    for (int i = 0; i < data_list.size(); ++i) static_data.push_back(Rcpp::as<arma::mat>(data_list[i]));
    arma::mat edgelist(0, 4);
    if (block.n_elem > 1) {
      REM tmp(edgelist, terms, n_nodes, directed, verbose, data_list, static_data, transformations, block, window_info);
      return tmp.sample(n_events, time, coef, degree_coef, max_events);
    } else {
      REM tmp(edgelist, terms, n_nodes, directed, verbose, data_list, static_data, transformations, window_info);
      return tmp.sample(n_events, time, coef, degree_coef, max_events);
    }
  }

  // 2. Handle time-varying case
  arma::mat edgelist(0, 4);
  std::vector<arma::mat> initial_data = update_object_covariates_get_mats(data_list, 0.0, transformations);
  REM tmp(edgelist, terms, n_nodes, directed, verbose, data_list, initial_data, transformations, block.n_elem > 1 ? block : arma::vec(), window_info);

  std::vector<arma::mat> res_list;
  double current_time = 0;
  unsigned int events_simulated = 0;
  auto cp_it = change_points.begin();

  while (events_simulated < max_events && (time == 0 || current_time < time)) {
    auto next_cp_it = std::next(cp_it);
    double next_cp_time = (next_cp_it == change_points.end()) ? R_PosInf : *next_cp_it;
    if (time > 0 && next_cp_time > time) next_cp_time = time;

    // Simulate for this interval
    double interval_duration = (next_cp_time == R_PosInf) ? 0 : (next_cp_time - current_time);
    unsigned int n_to_sim = (n_events > 0) ? (n_events - events_simulated) : 0;

    // We call sample. Note that sample internally resets time_old to 0.
    arma::mat interval_events = tmp.sample(n_to_sim, interval_duration, coef, degree_coef, max_events);

    if (interval_events.n_rows > 0) {
      // tmp.sample already returns absolute times.

      res_list.push_back(interval_events);
      events_simulated += interval_events.n_rows;
      // Now we need to update the history in tmp for the next interval
      // tmp.sample already updated the history for events it simulated.
    }
    current_time = next_cp_time;

    if ((time > 0 && current_time >= time) || (n_events > 0 && events_simulated >= n_events)) break;
    if (next_cp_it == change_points.end()) {
        if (n_events > 0 && events_simulated < n_events && time == 0) {
            // Keep going for one last infinite interval
        } else {
            break;
        }
    }
    cp_it = next_cp_it;
  }

  // Combine all simulated events
  if (res_list.empty()) return arma::mat(0, 3);
  unsigned int total_rows = 0;
  for (const auto& m : res_list) total_rows += m.n_rows;
  arma::mat final_events(total_rows, 3);
  unsigned int current_row = 0;
  for (const auto& m : res_list) {
    final_events.rows(current_row, current_row + m.n_rows - 1) = m;
    current_row += m.n_rows;
  }

  return final_events;
}

// [[Rcpp::export]]
arma::mat simulate_from_empty(std::vector<std::string> terms, Rcpp::List data_list,
                              std::vector<std::string> transformations,
                              unsigned int n_nodes, bool verbose, bool directed,
                              arma::vec coef_0_1, arma::vec coef_1_0,
                              arma::vec degree_coef_0_1, arma::vec degree_coef_1_0,
                              unsigned int n_events, double time,
                              unsigned int max_events, bool simultaneous_interactions,
                              Rcpp::List window_info,
                              unsigned int seed = 0) {
  if (seed > 0) {
      Rcpp::Environment base("package:base");
      Rcpp::Function set_seed = base["set.seed"];
      set_seed(seed);
  }
  std::set<double> change_points;
  change_points.insert(0.0);
  bool is_time_varying = false;
  for (int i = 0; i < data_list.size(); ++i) {
    if (Rcpp::is<Rcpp::List>(data_list[i])) {
      is_time_varying = true;
      Rcpp::List tv_list = Rcpp::as<Rcpp::List>(data_list[i]);
      Rcpp::CharacterVector times = tv_list.names();
      for (R_xlen_t k = 0; k < times.size(); ++k) {
        double tk = safe_stod(std::string(times[k]), "parsing change points");
        if (tk > 1e-10) change_points.insert(tk);
      }
    }
  }

  if (!is_time_varying) {
    std::vector<arma::mat> static_data;
    for (int i = 0; i < data_list.size(); ++i) static_data.push_back(Rcpp::as<arma::mat>(data_list[i]));
    arma::mat edgelist_empty(0, 4);
    DEM tmp(edgelist_empty, terms, n_nodes, directed, verbose, data_list, static_data, transformations, simultaneous_interactions, window_info);
    return tmp.sample(n_events, time, coef_0_1, coef_1_0, degree_coef_0_1, degree_coef_1_0, max_events);
  }

  // Handle time-varying case
  arma::mat edgelist_empty(0, 4);
  std::vector<arma::mat> initial_data = update_object_covariates_get_mats(data_list, 0.0, transformations);
  DEM tmp(edgelist_empty, terms, n_nodes, directed, verbose, data_list, initial_data, transformations, simultaneous_interactions, window_info);

  std::vector<arma::mat> res_list;
  double current_time = 0;
  unsigned int events_simulated = 0;
  auto cp_it = change_points.begin();

  while (events_simulated < max_events && (time == 0 || current_time < time)) {
    double cp_time = *cp_it;
    auto next_cp_it = std::next(cp_it);
    double next_cp_time = (next_cp_it == change_points.end()) ? R_PosInf : *next_cp_it;
    if (time > 0 && next_cp_time > time) next_cp_time = time;

    update_object_covariates(tmp, cp_time);

    double interval_duration = (next_cp_time == R_PosInf) ? 0 : (next_cp_time - current_time);
    unsigned int n_to_sim = (n_events > 0) ? (n_events - events_simulated) : 0;

    arma::mat interval_events = tmp.sample(n_to_sim, interval_duration, coef_0_1, coef_1_0, degree_coef_0_1, degree_coef_1_0, max_events);

    if (interval_events.n_rows > 0) {
      interval_events.col(0) += current_time;
      res_list.push_back(interval_events);
      events_simulated += interval_events.n_rows;
    }
    current_time = next_cp_time;

    if ((time > 0 && current_time >= time) || (n_events > 0 && events_simulated >= n_events)) break;
    if (next_cp_it == change_points.end()) {
        if (n_events > 0 && events_simulated < n_events && time == 0) {
            // Keep going for one last infinite interval
        } else {
            break;
        }
    }
    cp_it = next_cp_it;
  }

  if (res_list.empty()) return arma::mat(0, 4);
  unsigned int total_rows = 0;
  for (const auto& m : res_list) total_rows += m.n_rows;
  arma::mat final_events(total_rows, 4);
  unsigned int current_row = 0;
  for (const auto& m : res_list) {
    final_events.rows(current_row, current_row + m.n_rows - 1) = m;
    current_row += m.n_rows;
  }
  return final_events;
}

// [[Rcpp::export]]
Rcpp::List get_probabilities_per_test_event(
    std::vector<std::string> terms, Rcpp::List data_list,
                              std::vector<std::string> transformations,
                              unsigned int n_nodes,
                              bool verbose, bool directed,
                              arma::vec coef_0_1, arma::vec coef_1_0,
                              arma::vec degree_coef_0_1, arma::vec degree_coef_1_0,
                              bool simultaneous_interactions,
                              arma::mat edgelist_train,
                              arma::mat edgelist_test, unsigned int k,
                              bool is_rem,
                              Rcpp::List window_info,
                              Rcpp::Nullable<Rcpp::NumericVector> baseline_0_1 = R_NilValue,
                              Rcpp::Nullable<Rcpp::NumericVector> baseline_1_0 = R_NilValue) {

  arma::vec b_0_1;
  if (baseline_0_1.isNotNull()) {
    b_0_1 = Rcpp::as<arma::vec>(baseline_0_1);
  }
  if (b_0_1.is_empty()) {
    b_0_1 = arma::vec(1, arma::fill::zeros);
  }

  arma::vec b_1_0;
  if (baseline_1_0.isNotNull()) {
    b_1_0 = Rcpp::as<arma::vec>(baseline_1_0);
  }
  if (b_1_0.is_empty()) {
    b_1_0 = arma::vec(1, arma::fill::zeros);
  }

  if (is_rem) {
    double last_train_time = 0;
    if (edgelist_train.n_rows > 0) {
      last_train_time = edgelist_train.col(0).max();
    }
    std::vector<arma::mat> initial_data = update_object_covariates_get_mats(data_list, last_train_time, transformations);
    REM tmp(edgelist_train, terms, n_nodes, directed, verbose, data_list, initial_data, transformations, window_info);

    tmp.preprocess_build(0, tmp.changepoints.size());
    Rcpp::List res2 = tmp.get_probabilities_per_test_event(coef_0_1, degree_coef_0_1,
                                                           edgelist_test,
                                                           simultaneous_interactions, k,
                                                           b_0_1);
    return res2;
  } else {
    double last_train_time = 0;
    if (edgelist_train.n_rows > 0) {
      last_train_time = edgelist_train.col(0).max();
    }
    std::vector<arma::mat> initial_data = update_object_covariates_get_mats(data_list, last_train_time, transformations);
    DEM tmp(edgelist_train, terms, n_nodes, directed, verbose, data_list, initial_data, transformations, simultaneous_interactions, window_info);

    tmp.preprocess_build(0, tmp.changepoints.size());
    Rcpp::List res2 = tmp.get_probabilities_per_test_event(coef_0_1, degree_coef_0_1,
                                                           coef_1_0, degree_coef_1_0,
                                                           edgelist_test,
                                                           simultaneous_interactions, k,
                                                           b_0_1, b_1_0);
    return res2;
  }
}



// [[Rcpp::export]]
arma::vec get_oos_likelihood_cpp(
    std::vector<std::string> terms, Rcpp::List data_list,
                              std::vector<std::string> transformations,
                              unsigned int n_nodes,
                              bool verbose, bool directed,
                              arma::vec coef_0_1, arma::vec coef_1_0,
                              arma::vec degree_coef_0_1, arma::vec degree_coef_1_0,
                              bool simultaneous_interactions,
                              arma::mat edgelist_train,
                              arma::mat edgelist_test,
                              bool is_rem,
                              Rcpp::List window_info,
                              Rcpp::Nullable<Rcpp::NumericVector> baseline_0_1 = R_NilValue,
                              Rcpp::Nullable<Rcpp::NumericVector> baseline_1_0 = R_NilValue) {

  arma::vec b_0_1;
  if (baseline_0_1.isNotNull()) {
    b_0_1 = Rcpp::as<arma::vec>(baseline_0_1);
  }
  if (b_0_1.is_empty()) {
    b_0_1 = arma::vec(1, arma::fill::zeros);
  }

  arma::vec b_1_0;
  if (baseline_1_0.isNotNull()) {
    b_1_0 = Rcpp::as<arma::vec>(baseline_1_0);
  }
  if (b_1_0.is_empty()) {
    b_1_0 = arma::vec(1, arma::fill::zeros);
  }

  if (is_rem) {
    double last_train_time = 0;
    if (edgelist_train.n_rows > 0) {
      last_train_time = edgelist_train.col(0).max();
    }
    std::vector<arma::mat> initial_data = update_object_covariates_get_mats(data_list, last_train_time, transformations);
    REM tmp(edgelist_train, terms, n_nodes, directed, verbose, data_list, initial_data, transformations, window_info);

    tmp.preprocess_build(0, tmp.changepoints.size());
    arma::vec res2 = tmp.get_oos_likelihood(coef_0_1, degree_coef_0_1,
                                            edgelist_test,
                                            b_0_1);
    return res2;
  } else {
    double last_train_time = 0;
    if (edgelist_train.n_rows > 0) {
      last_train_time = edgelist_train.col(0).max();
    }
    std::vector<arma::mat> initial_data = update_object_covariates_get_mats(data_list, last_train_time, transformations);
    DEM tmp(edgelist_train, terms, n_nodes, directed, verbose, data_list, initial_data, transformations, simultaneous_interactions, window_info);

    tmp.preprocess_build(0, tmp.changepoints.size());
    arma::vec res2 = tmp.get_oos_likelihood(coef_0_1, degree_coef_0_1,
                                            coef_1_0, degree_coef_1_0,
                                            edgelist_test,
                                            simultaneous_interactions,
                                            b_0_1, b_1_0);
    return res2;
  }
}



// [[Rcpp::export]]
arma::vec update_core_cpp(const arma::mat& X,
                            const arma::vec& y,
                            const arma::vec& prediction,
                            arma::vec est_core,
                            const arma::uvec& identifiable,
                            const arma::vec& offset_fixed,
                            const arma::vec& weights) {
  if (identifiable.n_elem == 0) {
    return est_core;
  }

  double old_llh = eval_llh_pois(y, prediction, weights);
  if (!std::isfinite(old_llh)) old_llh = -1e150;

  arma::mat X_id = X.cols(identifiable);
  arma::vec p = prediction;
  for(unsigned int i=0; i<p.n_elem; ++i) {
    if(!std::isfinite(p[i])) p[i] = 1e150; // Use large value for Inf to force step down
    // Removed biased gradient hack that was interfering with NR link
  }

  arma::vec p_w = p;
  if (weights.n_elem == p.n_elem) p_w %= weights;

  arma::mat H = (X_id.each_col() % p_w).t() * X_id;
  arma::mat H_inv;
  if (!arma::solve(H_inv, H, arma::eye(H.n_rows, H.n_cols))) {
    H_inv = arma::pinv(H);
  }

  arma::vec res_w = (weights.n_elem == p.n_elem) ? weights : arma::ones<arma::vec>(p.n_elem);
  arma::vec grad = X_id.t() * ((y - p) % res_w);
  arma::vec delta_beta = H_inv * grad;

  for (int i = 0; i <= 10; ++i) {
    double step_size = std::pow(0.5, i);
    arma::vec new_beta = est_core;
    for(unsigned int j=0; j<identifiable.n_elem; ++j) {
       new_beta(identifiable(j)) += step_size * delta_beta(j);
    }

    arma::vec new_pred = arma::exp(offset_fixed + X * new_beta);
    double max_val = 0.0;
    bool any_finite = false;
    for(unsigned int k=0; k<new_pred.n_elem; ++k) {
      if(std::isfinite(new_pred[k])) {
        if(!any_finite || new_pred[k] > max_val) max_val = new_pred[k];
        any_finite = true;
      }
    }

    for(unsigned int k=0; k<new_pred.n_elem; ++k) {
      if(!std::isfinite(new_pred[k])) {
        new_pred[k] = any_finite ? max_val : 1e-10;
      }
    }

    double new_llh = eval_llh_pois(y, new_pred, weights);
    if (!std::isfinite(new_llh)) new_llh = -1e150;

    // Strict monotonicity check: we only accept the step if it improves the LLH.
    // We already have a line search loop (pow(0.5, i)).
    if (new_llh > old_llh) return new_beta;

    if (i == 10) {
       // If no improvement found even after 10 halvings, stay at old values.
       return est_core;
    }
  }
  return est_core;
}

struct Event {
  double time_val;
  int type; // +1 for start, -1 for end

  bool operator<(const Event& other) const {
    if (time_val != other.time_val) return time_val < other.time_val;
    return type < other.type;
  }
};
// [[Rcpp::export]]
Rcpp::DataFrame get_union_bounds(Rcpp::IntegerVector pair_id, Rcpp::NumericVector time_start, Rcpp::NumericVector time_end) {
  if (pair_id.size() != time_start.size() || pair_id.size() != time_end.size()) {
    Rcpp::stop("All input vectors to get_union_bounds must have the same length.");
  }
  int n_rows = pair_id.size();

  // Group purely by the dyad topology
  std::unordered_map<int, std::vector<Event>> event_grid;
  for (int i = 0; i < n_rows; ++i) {
    int p_id = pair_id[i];
    event_grid[p_id].push_back({time_start[i], 1});
    event_grid[p_id].push_back({time_end[i], -1});
  }

  std::vector<int> out_pair_id;
  std::vector<double> out_start;
  std::vector<double> out_end;

  for (auto& pair : event_grid) {
    int current_pair = pair.first;
    auto& events = pair.second;
    std::sort(events.begin(), events.end());

    int active_sources = 0;
    double current_start = -1;

    for (size_t i = 0; i < events.size(); ++i) {

      // If we are currently active, and the timeline advances, we record the slice
      if (active_sources > 0 && current_start != -1 && current_start < events[i].time_val) {
        out_pair_id.push_back(current_pair);
        out_start.push_back(current_start);
        out_end.push_back(events[i].time_val);
      }

      // Update the global active state tracker
      if (events[i].type == 1) active_sources++;
      else active_sources--;

      // Move the start boundary forward
      current_start = events[i].time_val;
    }
  }

  return Rcpp::DataFrame::create(Rcpp::_["pair_id"]  = out_pair_id,
                                 Rcpp::_["grid_start"] = out_start,
                                 Rcpp::_["grid_end"]   = out_end);
}

