#pragma once
#include <RcppArmadillo.h>
#include <iostream>

arma::mat combinations_directed(unsigned int x);
arma::mat combinations_undirected(unsigned int x);

class Current_Stat
{
public:
  // unsigned int n_nodes_, unsigned int n_statistics,
  // bool directed_, std::vector<arma::mat> eval_at_zero_
  Current_Stat(unsigned int, unsigned int, bool, std::vector<arma::mat> );
  arma::mat data;
  unsigned int n_nodes;
  bool directed;

  std::vector<uint8_t> is_interacting;
  std::vector<arma::mat> eval_at_zero;

  // We define member function that only work on current_stat objects
  // these functions let us find specific row/rows where the from/to nodes are specified
  arma::uvec find_from(unsigned int);
  arma::uvec find_involved(unsigned int);
  arma::mat return_data();
  arma::uvec find_to(unsigned int);
  unsigned int find_from_to(unsigned int,unsigned  int);
  arma::uvec find_froms(const arma::uvec&);
  arma::uvec find_tos(const arma::uvec&);
  arma::uvec find_from_tos(const arma::uvec&, const arma::uvec&);
  arma::uvec get_currently_interacting();
  arma::uvec get_currently_noninteracting();
  unsigned int get_status(unsigned int, unsigned int);

  // Set the status of pairs to be either 1 (interacting) or 0 (not interacting)
  void set_status(unsigned int, unsigned int, unsigned int);
  void set_event(unsigned int, unsigned int);
  void reset_event(unsigned int, unsigned int);

  arma::vec update_intensity(const arma::vec&, const arma::vec&);

  // Vectorized helpers for status change tracking
  arma::uvec now_avail(unsigned int);
  arma::uvec not_avail(unsigned int);

  inline void add_stats_to_col(const arma::uvec& indices, unsigned int col_idx, double val) {
    if(!indices.is_empty()) {
      data.submat(indices, arma::uvec{col_idx}) += val;
    }
  }
  inline void set_stats_in_col(const arma::uvec& indices, unsigned int col_idx, double val) {
    if(!indices.is_empty()) {
      data.submat(indices, arma::uvec{col_idx}).fill(val);
    }
  }
  inline void log_add_stats_to_col(const arma::uvec& indices, unsigned int col_idx, double val) {
    if(!indices.is_empty()) {
      for (unsigned int i = 0; i < indices.n_elem; ++i) {
        unsigned int idx = indices[i];
        double v = data.at(idx, col_idx);
        if (v < 0.0) v = 0.0;
        double arg = val * std::exp(-v);
        if (arg <= -1.0) {
          data.at(idx, col_idx) = -std::numeric_limits<double>::infinity();
        } else {
          double next_v = v + std::log1p(arg);
          if (next_v < 0.0 && next_v > -1e-12) next_v = 0.0;
          data.at(idx, col_idx) = next_v;
        }
      }
    }
  }
  inline void recip_add_stats_to_col(const arma::uvec& indices, unsigned int col_idx, double val) {
    if(!indices.is_empty()) {
      for (unsigned int i = 0; i < indices.n_elem; ++i) {
        unsigned int idx = indices[i];
        double v = data.at(idx, col_idx);
        if (v == 0.0) {
           data.at(idx, col_idx) = 1.0 / (1.0 + val);
        } else {
           data.at(idx, col_idx) = 1.0 / (1.0/v + val);
        }
      }
    }
  }
  inline void sigmoid_add_stats_to_col(const arma::uvec& indices, unsigned int col_idx, double val, double K) {
    if(!indices.is_empty()) {
      for (unsigned int i = 0; i < indices.n_elem; ++i) {
        unsigned int idx = indices[i];
        double v = data.at(idx, col_idx);
        data.at(idx, col_idx) = (K * v + val * (1.0 - v)) / (K + val * (1.0 - v));
      }
    }
  }

  void set_stat(unsigned int, unsigned int, double);
  bool update_exogenous(unsigned int term_idx, const arma::mat& new_data);
  void update_exogenous(const std::vector<arma::mat>& new_eval_at_zero);

  void clear();

  // Caching for scalable intensity updates (Simulation)
  arma::vec lp_0_1;
  arma::vec lp_1_0;
  arma::vec exp_lp_0_1;
  arma::vec exp_lp_1_0;
  arma::vec lp_degree_0_1;
  arma::vec lp_degree_1_0;
  
  // Fully stateful intensity tracking
  arma::vec combined_intensity;
  arma::uvec currently_interacting;
  arma::uvec currently_noninteracting;
  bool simultaneous_interactions;
  // Current baseline values (stored so update_baseline can compute deltas)
  double current_baseline_0_1 = 0.0;
  double current_baseline_1_0 = 0.0;

  void initialize_intensities(const arma::vec& coef_0_1, const arma::vec& coef_degree_0_1, 
                              const arma::vec& coef_1_0, const arma::vec& coef_degree_1_0,
                              bool simultaneous_interactions_,
                              double baseline_0_1_ = 0.0, double baseline_1_0_ = 0.0);
  void update_baseline(double new_baseline_0_1, double new_baseline_1_0,
                       const arma::vec& coef_0_1, const arma::vec& coef_1_0);
  void update_intensities_at_indices(const arma::uvec& indices, const arma::vec& coef_0_1, const arma::vec& coef_1_0);
  void refresh_combined_intensity(const arma::uvec& indices);

  // Cached indices for performance
  arma::uvec cache_dyad_u; 
  arma::uvec cache_dyad_v; 
  std::vector<arma::uvec> sender_to_dyads;
  std::vector<arma::uvec> receiver_to_dyads;
  std::vector<arma::uvec> node_to_dyads;

protected:
private:
};
