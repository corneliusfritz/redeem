# pragma once
#include <queue>
#include <random>
#include <set>
#include <unordered_map>
#include <RcppArmadillo.h>
#include "data_dem.h"
#include "sufficient_statistics.h"

// Forward declaration of ValidateFunction type
// Note: ValidateFunction and change_statistics_generate are now defined in sufficient_statistics.h

struct ScheduledEvent {
    double time;
    unsigned int from;
    unsigned int to;
    unsigned int type;

    bool operator>(const ScheduledEvent& other) const {
        return time > other.time;
    }
};


// Helper function to safely parse doubles from strings (e.g., covariate measurement times)
// and throw a descriptive R error if the format is invalid.
inline double safe_stod(const std::string& s, const std::string& context = "") {
  std::string trimmed = s;
  trimmed.erase(0, trimmed.find_first_not_of(" \t\n\r"));
  trimmed.erase(trimmed.find_last_not_of(" \t\n\r") + 1);
  
  if (trimmed.empty()) {
    Rcpp::stop("Empty measurement time found for " + context);
  }

  try {
    std::size_t idx = 0;
    double val = std::stod(trimmed, &idx);
    if (idx != trimmed.size()) {
       throw std::invalid_argument("trailing characters");
    }
    return val;
  } catch (const std::exception& e) {
    std::string msg = "All names of the 'data' list must be numeric";
    if (!context.empty()) msg += " while " + context;
    msg += ": '" + s + "'. (" + e.what() + "). Please ensure all names in your covariate data list reflect valid numeric measurement times.";
    Rcpp::stop(msg);
    return 0.0; // Unreachable
  }
}

// Apply transformation to a covariate matrix (for exogenous covariates)
inline arma::mat apply_transformation(const arma::mat& x, const std::string& transformation, double K) {
  if (transformation == "identity") {
    return x;
  } else if (transformation == "log") {
    arma::mat res = arma::log(1.0 + x);
    if (x.n_elem > 0) res.at(0, 0) = x.at(0, 0);
    return res;
  } else if (transformation == "recip") {
    arma::mat res = 1.0 / (1.0 + x);
    if (x.n_elem > 0) res.at(0, 0) = x.at(0, 0);
    return res;
  } else if (transformation == "bin") {
    arma::mat res = arma::zeros(x.n_rows, x.n_cols);
    res.elem(arma::find(x > 0)).ones();
    if (x.n_elem > 0) res.at(0, 0) = x.at(0, 0);
    return res;
  } else if (transformation == "sig") {
    if (K == 0) K = 1.0;
    arma::mat res = x / (x + K);
    if (x.n_elem > 0) res.at(0, 0) = x.at(0, 0);
    return res;
  } else {
    return x;
  }
}

// Return top-K indices by value (descending), sorted by value
inline arma::uvec topk_heap(const arma::vec& v, arma::uword K) {
  using U = arma::uword;
  using P = std::pair<double,U>;  // (value, index)

  K = std::min<U>(K, v.n_elem);
  auto cmp = [](const P& a, const P& b){ return a.first > b.first; }; // min-heap
  std::priority_queue<P, std::vector<P>, decltype(cmp)> pq(cmp);

  for (U i = 0; i < v.n_elem; ++i) {
    const double val = v[i];
    if (pq.size() < K) pq.emplace(val, i);
    else if (val > pq.top().first) { pq.pop(); pq.emplace(val, i); }
  }
  arma::uvec idx(K);
  for (arma::sword t = K - 1; t >= 0; --t) { idx[t] = pq.top().second; pq.pop(); }
  const double* x = v.memptr();
  std::sort(idx.begin(), idx.end(), [&](U i, U j){ return x[i] > x[j]; });
  return idx;
}

// Robust log-sum-exp helper function for a vector of linear predictors
inline double log_sum_exp(const arma::vec& lp) {
  double max_val = -arma::datum::inf;
  for (arma::uword i = 0; i < lp.n_elem; ++i) {
    if (lp[i] > max_val && std::isfinite(lp[i])) {
      max_val = lp[i];
    }
  }
  if (max_val == -arma::datum::inf) {
    return -arma::datum::inf;
  }
  double sum_exp = 0.0;
  for (arma::uword i = 0; i < lp.n_elem; ++i) {
    if (std::isfinite(lp[i])) {
      sum_exp += std::exp(lp[i] - max_val);
    }
  }
  return max_val + std::log(sum_exp);
}

// Declaration only - implementation and [[Rcpp::export]] in helper_impl.cpp
arma::uvec topk_indices(const arma::vec& x, arma::uword k);

// Advance time-dependent statistics at a changepoint (e.g., for current_interaction)
inline arma::uvec advance_delta_time(Data_DEM &object, 
                                     unsigned int number_event, 
                                     std::vector<arma::mat> &data_list,
                                     std::vector<std::string> &transformations,
                                     std::vector<ValidateFunction> &terms) {
  arma::uvec all_stat_changes;
  if (!terms.empty()) {
    unsigned int from_dummy = 0;
    unsigned int to_dummy = 0;
    for (unsigned int a = 0; a < (unsigned int)terms.size(); ++a) {
      // Passing from = 0 signals to the statistic function that it should only perform time-advancement
      // (as implemented in stat_current_interaction)
      arma::uvec changes = terms.at(a)(object, data_list.at(a), from_dummy, to_dummy, number_event, a + 7, transformations.at(a), true);
      if (!changes.is_empty()) {
        all_stat_changes = arma::join_cols(all_stat_changes, changes);
      }
    }
  }
  return arma::unique(all_stat_changes);
}

inline void advance_delta_time(Data_DEM &object, 
                                      unsigned int number_event, 
                                      std::vector<arma::mat> &data_list,
                                      std::vector<std::string> &transformations,
                                      std::vector<ValidateFunction> &terms,
                                      arma::uvec &changes_out) {
  changes_out = advance_delta_time(object, number_event, data_list, transformations, terms);
}

// Function to call all functions in the vector functions and return the changes
inline arma::uvec compute_changes(unsigned int from,
                           unsigned int to,
                           unsigned int number_event,
                           unsigned int type,
                           Data_DEM &object,
                           std::vector<arma::mat> &data_list,
                           std::vector<std::string> &transformation,
                           std::vector<ValidateFunction> &functions,
                           bool simultaneous_interactions){

  const unsigned int n_rows = object.current_stats.data.n_rows;
  std::vector<uint8_t> changed(n_rows, 0);

  // The pair where the event took place always changes
  unsigned int event_pair_idx = object.current_stats.find_from_to(from, to);
  changed[event_pair_idx] = 1;

  // Change object.current_stats to indicate that this pair experienced an event
  object.current_stats.set_event(from,to);

  
  if (functions.size() > 0) {
    for(unsigned int a = 0; a < functions.size(); ++a) {
      
      // skip NON-windowed statistics if the event is a virtual window expiration (type >= 10)
      if (type >= 10 && object.term_names[a].find("_wt") == std::string::npos) continue;
      // Offset 7 points to the first statistics column
      arma::uvec change_stat = functions.at(a)(object, data_list.at(a), from, to, number_event, a + 7, transformation.at(a), type);
      for(unsigned int i = 0; i < change_stat.n_elem; ++i) {
        changed[change_stat[i]] = 1;
      }
    }
  }
  
  // If simultaneous interactions are not possible, the availability of nodes changes
  if(!simultaneous_interactions){
    arma::uvec from_affected = object.current_stats.find_involved(from);
    arma::uvec to_affected = object.current_stats.find_involved(to);
    
    for(unsigned int i = 0; i < from_affected.n_elem; ++i) changed[from_affected[i]] = 1;
    for(unsigned int i = 0; i < to_affected.n_elem; ++i) changed[to_affected[i]] = 1;
  }
  
  // Collect all unique changed indices
  unsigned int count = 0;
  for(uint8_t val : changed) if(val) count++;
  
  arma::uvec all_changes(count);
  unsigned int idx = 0;
  for(unsigned int i = 0; i < n_rows; ++i) {
    if(changed[i]) all_changes[idx++] = i;
  }
  
  if (type == 1) {
    object.last_sender = from;
    object.last_receiver = to;
  }
  
  return all_changes;
}

inline void compute_changes(unsigned int from,
                            unsigned int to,
                            unsigned int number_event,
                            unsigned int type,
                            Data_DEM &object,
                            std::vector<arma::mat> &data_list,
                            std::vector<std::string> &transformation,
                            std::vector<ValidateFunction> &functions,
                            bool simultaneous_interactions,
                            arma::uvec &changes_out){
  changes_out = compute_changes(from, to, number_event, type, object, data_list, transformation, functions, simultaneous_interactions);
}

// Replaces the R function do_call_full with a native C++ implementation
// Logic: intervals are formed between state changes for each dyad.
arma::mat do_call_full_cpp(Rcpp::List res, double max_time);
