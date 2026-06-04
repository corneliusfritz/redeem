#pragma once
#include <RcppArmadillo.h>
#include <set>
#include <string>
#include <vector>
#include <algorithm>
#include "data_dem.h"
#include "redeem/extension_api.hpp"

using namespace Rcpp;

// Functional pointer for dispatching statistics calculations based on term names.
using ValidateFunction = arma::uvec (*)(Data_DEM &, arma::mat &, unsigned int &, unsigned int &, unsigned int &, unsigned int, std::string, unsigned int);

// --- Helpers ---

// Helper to select the correct history based on term name
inline Hist_Events* select_history(Data_DEM &object, unsigned int col_number, bool &is_windowed, int &window_size) {
  window_size = -1;
  is_windowed = false;
  
  if (object.term_names.empty()) return &(object.general_interactions);
  
  // The internal data matrix in Data_DEM always has 7 fixed columns (PairID, Status, Event, From, To, Avail_i, Avail_j).
  // Statistics start at index 7, so term_names[col_number - 7] gives the correct name for the statistic.
  if (col_number >= 7 && (col_number - 7) < object.term_names.size()) {
    const std::string &tname = object.term_names[col_number - 7];
    size_t pos = tname.find("_wt");
    if (pos != std::string::npos) {
      double wt = std::stod(tname.substr(pos + 3));
      for (auto const& [id, length] : object.window_lengths) {
        if (std::abs(length - wt) < 1e-6) {
          is_windowed = true;
          window_size = id;
          return &(object.windowed_history.at(id));
        }
      }
    }
  }

  return &(object.general_interactions);
}

// Helper to get raw statistic value for a given dyad
inline double get_raw_statistic_value(Data_DEM &object, unsigned int idx, unsigned int col_number) {
  unsigned int u = (unsigned int)object.current_stats.data.at(idx, 3);
  unsigned int v = (unsigned int)object.current_stats.data.at(idx, 4);
  const std::string &tname = object.term_names[col_number - 7];

  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);

  if (tname.find("current_interaction") != std::string::npos || tname.find("duration") != std::string::npos) {
    return object.current_stats.data.at(idx, col_number);
  }

  bool is_current = (tname.find("current_") == 0);
  Hist_Events* hist_to_use = is_current ? &(object.current_interactions) : history;

  if (tname.find("inertia") != std::string::npos || tname.find("number_interaction") != std::string::npos) {
    return hist_to_use->get_count(u, v);
  }
  if (tname.find("reciprocity") != std::string::npos) {
    return hist_to_use->get_count(v, u);
  }
  if (tname.find("degree_out_sender") != std::string::npos) {
    return hist_to_use->get_degree(u, "out", false);
  }
  if (tname.find("degree_out_receiver") != std::string::npos) {
    return hist_to_use->get_degree(v, "out", false);
  }
  if (tname.find("degree_in_sender") != std::string::npos) {
    return hist_to_use->get_degree(u, "in", false);
  }
  if (tname.find("degree_in_receiver") != std::string::npos) {
    return hist_to_use->get_degree(v, "in", false);
  }
  if (tname.find("degree_sum") != std::string::npos) {
    return hist_to_use->get_degree(u, "out", false) + hist_to_use->get_degree(v, "out", false);
  }
  if (tname.find("degree_absdiff") != std::string::npos) {
    return std::abs((double)hist_to_use->get_degree(u, "out", false) - (double)hist_to_use->get_degree(v, "out", false));
  }
  if (tname.find("count_out_sender") != std::string::npos) {
    return hist_to_use->get_degree(u, "out", true);
  }
  if (tname.find("count_out_receiver") != std::string::npos) {
    return hist_to_use->get_degree(v, "out", true);
  }
  if (tname.find("count_in_sender") != std::string::npos) {
    return hist_to_use->get_degree(u, "in", true);
  }
  if (tname.find("count_in_receiver") != std::string::npos) {
    return hist_to_use->get_degree(v, "in", true);
  }
  if (tname.find("count_sum") != std::string::npos) {
    return hist_to_use->get_degree(u, "out", true) + hist_to_use->get_degree(v, "out", true);
  }
  if (tname.find("count_absdiff") != std::string::npos) {
    return std::abs((double)hist_to_use->get_degree(u, "out", true) - (double)hist_to_use->get_degree(v, "out", true));
  }
  if (tname.find("common_partner") != std::string::npos) {
    std::string type = "OSP";
    if (tname.find("_ISP") != std::string::npos) type = "ISP";
    else if (tname.find("_OTP") != std::string::npos) type = "OTP";
    else if (tname.find("_ITP") != std::string::npos) type = "ITP";
    return hist_to_use->get_common_partners(u, v, type).n_elem;
  }
  if (tname.find("triangle") != std::string::npos) {
    std::string type = "OSP";
    if (tname.find("_ISP") != std::string::npos) type = "ISP";
    else if (tname.find("_OTP") != std::string::npos) type = "OTP";
    else if (tname.find("_ITP") != std::string::npos) type = "ITP";
    if (!hist_to_use->have_interacted(u, v)) return 0.0;
    return hist_to_use->get_common_partners(u, v, type).n_elem;
  }

  return object.current_stats.data.at(idx, col_number);
}

// Update helper
inline void apply_update(Data_DEM &object, const arma::uvec &indices, unsigned int col_number, double val, std::string transformation, arma::mat &data) {
  if (indices.is_empty()) return;
  if (transformation == "identity") {
    object.current_stats.add_stats_to_col(indices, col_number, val);
  } else if (transformation == "bin") {
    if (val < 0.0) {
      for (unsigned int i = 0; i < indices.n_elem; ++i) {
        unsigned int idx = indices[i];
        double raw_val = get_raw_statistic_value(object, idx, col_number) + val;
        if (raw_val <= 0.0) {
          arma::uvec idx_vec = {idx};
          object.current_stats.set_stats_in_col(idx_vec, col_number, 0.0);
        }
      }
    } else {
      object.current_stats.set_stats_in_col(indices, col_number, 1.0);
    }
  } else if (transformation == "log") {
    object.current_stats.log_add_stats_to_col(indices, col_number, val);
  } else if (transformation == "recip") {
    object.current_stats.recip_add_stats_to_col(indices, col_number, val);
  } else if (transformation == "sig") {
    double K = (data.n_elem > 0) ? data.at(0, 0) : 1.0;
    if (K == 0) K = 1.0;
    object.current_stats.sigmoid_add_stats_to_col(indices, col_number, val, K);
  }
}


double get_stat_move(unsigned int type, bool is_windowed, int window_size);

// --- Generator ---

inline std::vector<ValidateFunction> change_statistics_generate(std::vector<std::string> terms) {
  std::vector<ValidateFunction> functions(terms.size());
  auto& reg = redeem::Registry::instance();
  std::vector<std::string> reg_names = reg.names();
  // Sort keys by length descending to match more specific terms first
  std::sort(reg_names.begin(), reg_names.end(), [](const std::string& a, const std::string& b) {
    return a.length() > b.length();
  });

  for (unsigned int i = 0; i < terms.size(); i++) {
    const std::string &t = terms[i];
    if (t.empty()) {
      functions[i] = reg.get("Intercept");
      continue;
    }
    ValidateFunction found_func = nullptr;
    for (const std::string& reg_name : reg_names) {
      if (t.find(reg_name) != std::string::npos) {
        found_func = reg.get(reg_name);
        break;
      }
    }
    if (found_func) {
      functions[i] = found_func;
    } else {
      functions[i] = reg.get("dyadic_cov");
    }
  }
  return functions;
}
