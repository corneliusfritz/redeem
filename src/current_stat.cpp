#include "current_stat.h"
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <stdexcept>

// [[Rcpp::export]]
arma::mat combinations_directed(unsigned int x) {
  arma::vec from(x*(x-1));
  arma::vec to(x*(x-1));
  int n = 0;
  for(unsigned int a = 1; a <= x; a++) {
    for(unsigned int b = 1; b <= x; b++) {
      if( b!= a) {
        from.at(n)= a;
        to.at(n) = b;
        n ++;
      }
    }
  }
  return arma::join_rows(from,to);
}

// [[Rcpp::export]]
arma::mat combinations_undirected(unsigned int x) {
  arma::vec from(x*(x-1)/2);
  arma::vec to(x*(x-1)/2);
  int n = 0;
  for(unsigned int a = 1; a < x; a++) {
    for(unsigned int b = a + 1; b <= x; b++) {
      from.at(n)= a;
      to.at(n) = b;
      n ++;
    }
  }
  return arma::join_rows(from,to);
}

void check_uvec(const arma::umat& to_test, unsigned int max_row, unsigned int max_column) {
  if(to_test.n_elem > 0) {
    if (to_test.n_cols != 2) {
      Rcpp::stop("There need to be 2 columns in the parameter not " + std::to_string(to_test.n_cols));
    }
    unsigned int max_to_test_row = arma::max(to_test.col(0));
    unsigned int max_to_test_column = arma::max(to_test.col(1));
    if (max_to_test_row >= max_row) {
      Rcpp::stop("One of the rows is too high. Row number " + std::to_string(max_to_test_row) +
                  " but only " + std::to_string(max_row) + " provided.");
    }
    if (max_to_test_column >= max_column) {
      Rcpp::stop("One of the columns is too high. Column number " + std::to_string(max_to_test_column) +
                  " but only " + std::to_string(max_column) + " provided.");
    }
  }
}

Current_Stat::Current_Stat(unsigned int n_nodes_, unsigned int n_statistics,
                           bool directed_, std::vector<arma::mat> eval_at_zero_) {
  n_nodes = n_nodes_;
  directed = directed_;
  eval_at_zero = eval_at_zero_;

  arma::mat combinations;
  if (directed) {
    combinations = combinations_directed(n_nodes);
  } else {
    combinations = combinations_undirected(n_nodes);
  }

  unsigned int n_rows = combinations.n_rows;

  // Data matrix structure (7 fixed columns):
  // 0. Pair Id
  // 1. Status (0: not interacting, 1: interacting)
  // 2. Event (0: no recent event, 1: recent event)
  // 3. Node i
  // 4. Node j
  // 5. Is node i available?
  // 6. Is node j available?
  // 7+ Sufficient Statistics

  data.set_size(n_rows, 7 + n_statistics);
  data.col(0) = arma::regspace(1, n_rows);
  data.col(1).zeros(); // Status
  data.col(2).zeros(); // Event
  data.cols(3, 4) = combinations;
  data.col(5).ones(); // Avail i
  data.col(6).ones(); // Avail j
  if (n_statistics > 0) {
    data.cols(7, 6 + n_statistics).zeros();
    // FILL the statistics columns from eval_at_zero
    for (unsigned int i = 0; i < eval_at_zero.size(); ++i) {
        if (eval_at_zero[i].n_rows == n_nodes && eval_at_zero[i].n_cols == n_nodes) {
            arma::umat positions = arma::conv_to<arma::umat>::from(combinations - 1);
            arma::uvec eids = sub2ind(arma::size(eval_at_zero[i]), positions.t());
            data.col(7 + i) = eval_at_zero[i].elem(eids);
        }
    }
  }

  is_interacting.assign(n_rows, 0);

  // Initialize caches
  cache_dyad_u = arma::conv_to<arma::uvec>::from(data.col(3)) - 1;
  cache_dyad_v = arma::conv_to<arma::uvec>::from(data.col(4)) - 1;

  sender_to_dyads.resize(n_nodes);
  receiver_to_dyads.resize(n_nodes);
  node_to_dyads.resize(n_nodes);

  if (directed) {
    for (unsigned int i = 0; i < n_nodes; ++i) {
      // Senders are contiguous: [(N-1)*i, (N-1)*(i+1) - 1]
      sender_to_dyads[i] = arma::regspace<arma::uvec>(i * (n_nodes - 1), (i + 1) * (n_nodes - 1) - 1);

      // Receivers are scattered: one row index per sender block
      receiver_to_dyads[i].set_size(n_nodes - 1);
      unsigned int count = 0;
      for (unsigned int u = 0; u < n_nodes; ++u) {
        if (u == i) continue;
        if (u > i) receiver_to_dyads[i][count++] = u * (n_nodes - 1) + i;
        else receiver_to_dyads[i][count++] = u * (n_nodes - 1) + i - 1;
      }

      // Node dyads: union of sender and receiver roles (keep sorted)
      node_to_dyads[i] = arma::join_cols(sender_to_dyads[i], receiver_to_dyads[i]);
      node_to_dyads[i] = arma::sort(node_to_dyads[i]);
    }
  } else {
    // Undirected: use temporary std::vectors to build indices without O(N) reallocations
    std::vector<std::vector<unsigned int>> s_temp(n_nodes), r_temp(n_nodes), n_temp(n_nodes);
    for (unsigned int i = 0; i < n_rows; ++i) {
      unsigned int u = cache_dyad_u[i];
      unsigned int v = cache_dyad_v[i];
      s_temp[u].push_back(i);
      r_temp[v].push_back(i);
      n_temp[u].push_back(i);
      n_temp[v].push_back(i);
    }
    for (unsigned int i = 0; i < n_nodes; ++i) {
      sender_to_dyads[i] = arma::conv_to<arma::uvec>::from(s_temp[i]);
      receiver_to_dyads[i] = arma::conv_to<arma::uvec>::from(r_temp[i]);
      node_to_dyads[i] = arma::conv_to<arma::uvec>::from(n_temp[i]);
      // Note: n_temp is already sorted because we iterate i from 0 to n_rows
    }
  }
}

void Current_Stat::update_exogenous(const std::vector<arma::mat>& new_eval_at_zero) {
    this->eval_at_zero = new_eval_at_zero;
    if (eval_at_zero.size() > 0) {
        arma::mat combinations = data.cols(3, 4);
        arma::umat positions = arma::conv_to<arma::umat>::from(combinations - 1);
        for (unsigned int n = 0; n < eval_at_zero.size(); n++) {
            const arma::mat& new_mat = eval_at_zero.at(n);
            if (new_mat.n_rows == n_nodes && new_mat.n_cols == n_nodes) {
                arma::uvec eids = sub2ind(size(new_mat), positions.t());
                data.col(7 + n) = new_mat.elem(eids);
            }
        }
    }
}

bool Current_Stat::update_exogenous(unsigned int term_idx, const arma::mat& new_data) {
    if (term_idx >= eval_at_zero.size()) return false;

    // Check for change
    if (arma::size(this->eval_at_zero[term_idx]) == arma::size(new_data)) {
        if (arma::all(arma::vectorise(this->eval_at_zero[term_idx] == new_data))) {
            return false;
        }
    }

    this->eval_at_zero[term_idx] = new_data;

    if (new_data.n_rows == n_nodes && new_data.n_cols == n_nodes) {
        arma::mat combinations = data.cols(3, 4);
        arma::umat positions = arma::conv_to<arma::umat>::from(combinations - 1);
        if (7 + term_idx >= data.n_cols) {
            return true;
        }
        arma::uvec eids = sub2ind(size(new_data), positions.t());
        data.col(7 + term_idx) = new_data.elem(eids);
        return true;
    }
    return false;
}

arma::uvec Current_Stat::find_from(unsigned int from) {
  if(from > n_nodes || from == 0) Rcpp::stop("Node index must be between 1 and n_nodes");
  return sender_to_dyads[from - 1];
}

arma::uvec Current_Stat::find_to(unsigned int to) {
  if(to > n_nodes || to == 0) Rcpp::stop("Node index must be between 1 and n_nodes");
  return receiver_to_dyads[to - 1];
}

arma::uvec Current_Stat::find_involved(unsigned int node) {
  if(node > n_nodes || node == 0) Rcpp::stop("Node index must be between 1 and n_nodes");
  return node_to_dyads[node - 1];
}

unsigned int Current_Stat::find_from_to(unsigned int from, unsigned int to) {
  if(from == to) Rcpp::stop("Self-loops are not allowed.");
  if(from > n_nodes || to > n_nodes || from == 0 || to == 0) Rcpp::stop("Node index out of bounds.");

  if(directed) {
    if(from > to) return (n_nodes - 1) * (from - 1) + to - 1;
    else return (n_nodes - 1) * (from - 1) + to - 2;
  } else {
    if(from > to) std::swap(from, to);
    return n_nodes * (from - 1) - (from * (from - 1) / 2) + (to - from) - 1;
  }
}

arma::uvec Current_Stat::find_froms(const arma::uvec& froms) {
  unsigned int total_size = 0;
  for(unsigned int i = 0; i < froms.n_elem; ++i) total_size += sender_to_dyads[froms[i]-1].n_elem;
  arma::uvec res(total_size);
  unsigned int current_idx = 0;
  for(unsigned int i = 0; i < froms.n_elem; ++i) {
    const arma::uvec& v = sender_to_dyads[froms[i]-1];
    res.subvec(current_idx, current_idx + v.n_elem - 1) = v;
    current_idx += v.n_elem;
  }
  return res;
}

arma::uvec Current_Stat::find_tos(const arma::uvec& tos) {
  unsigned int total_size = 0;
  for(unsigned int i = 0; i < tos.n_elem; ++i) total_size += receiver_to_dyads[tos[i]-1].n_elem;
  arma::uvec res(total_size);
  unsigned int current_idx = 0;
  for(unsigned int i = 0; i < tos.n_elem; ++i) {
    const arma::uvec& v = receiver_to_dyads[tos[i]-1];
    res.subvec(current_idx, current_idx + v.n_elem - 1) = v;
    current_idx += v.n_elem;
  }
  return res;
}

arma::uvec Current_Stat::find_from_tos(const arma::uvec& froms, const arma::uvec& tos) {
  if(froms.n_elem != tos.n_elem) Rcpp::stop("froms and tos must have same length");
  arma::uvec res(froms.n_elem);
  for(unsigned int i = 0; i < froms.n_elem; ++i) {
    res[i] = find_from_to(froms[i], tos[i]);
  }
  return res;
}

arma::uvec Current_Stat::get_currently_interacting() {
  arma::uvec res(data.n_rows);
  unsigned int count = 0;
  for(unsigned int i = 0; i < data.n_rows; ++i) {
    if(is_interacting[i]) res[count++] = i;
  }
  return res.head(count);
}

arma::uvec Current_Stat::get_currently_noninteracting() {
  arma::uvec res(data.n_rows);
  unsigned int count = 0;
  for(unsigned int i = 0; i < data.n_rows; ++i) {
    if(!is_interacting[i]) res[count++] = i;
  }
  return res.head(count);
}

void Current_Stat::set_status(unsigned int from, unsigned int to, unsigned int status) {
  unsigned int pair = find_from_to(from, to);
  data(pair, 1) = status;
  is_interacting[pair] = (uint8_t)status;

  arma::uvec idx_vec = {pair};
  refresh_combined_intensity(idx_vec);
}

unsigned int Current_Stat::get_status(unsigned int from, unsigned int to) {
  return (unsigned int)is_interacting[find_from_to(from, to)];
}

void Current_Stat::set_event(unsigned int from, unsigned int to) {
  data(find_from_to(from, to), 2) = 1;
}

void Current_Stat::reset_event(unsigned int from, unsigned int to) {
  data(find_from_to(from, to), 2) = 0;
}

arma::uvec Current_Stat::not_avail(unsigned int node) {
  arma::uvec i_indices = find_from(node);
  if(!i_indices.is_empty()) data.submat(i_indices, arma::uvec{5}).zeros();

  arma::uvec j_indices = find_to(node);
  if(!j_indices.is_empty()) data.submat(j_indices, arma::uvec{6}).zeros();

  arma::uvec affected = arma::join_cols(i_indices, j_indices);
  refresh_combined_intensity(affected);
  return affected;
}

arma::uvec Current_Stat::now_avail(unsigned int node) {
  arma::uvec i_indices = find_from(node);
  if(!i_indices.is_empty()) data.submat(i_indices, arma::uvec{5}).ones();

  arma::uvec j_indices = find_to(node);
  if(!j_indices.is_empty()) data.submat(j_indices, arma::uvec{6}).ones();

  arma::uvec affected = arma::join_cols(i_indices, j_indices);
  refresh_combined_intensity(affected);
  return affected;
}

arma::vec Current_Stat::update_intensity(const arma::vec& coef, const arma::vec& coef_degree) {
  unsigned int n_stats = data.n_cols - 7;
  arma::vec lp = arma::zeros(data.n_rows);

  if (n_stats > 0 && coef.n_elem == n_stats) {
    lp = data.cols(7, 6 + n_stats) * coef;
  }

  if (coef_degree.n_elem > 0) {
    if (coef_degree.n_elem == 1) {
      lp += 2.0 * coef_degree.at(0);
    } else {
      if (directed) {
        lp += coef_degree.elem(cache_dyad_u) + coef_degree.elem(n_nodes + cache_dyad_v);
      } else {
        lp += coef_degree.elem(cache_dyad_u) + coef_degree.elem(cache_dyad_v);
      }
    }
  }

  arma::vec avail = data.col(5) % data.col(6);
  return arma::exp(lp) % avail;
}

void Current_Stat::initialize_intensities(const arma::vec& coef_0_1, const arma::vec& coef_degree_0_1,
                                          const arma::vec& coef_1_0, const arma::vec& coef_degree_1_0,
                                          bool simultaneous_interactions_,
                                          double baseline_0_1_, double baseline_1_0_) {
  unsigned int n_rows = data.n_rows;
  simultaneous_interactions = simultaneous_interactions_;

  // 1. Calculate degree parts (static during sampling)
  lp_degree_0_1 = arma::zeros(n_rows);
  if (coef_degree_0_1.n_elem == 1) {
    lp_degree_0_1.fill(2.0 * coef_degree_0_1.at(0));
  } else if (coef_degree_0_1.n_elem > 1) {
    if (directed) lp_degree_0_1 = coef_degree_0_1.elem(cache_dyad_u) + coef_degree_0_1.elem(n_nodes + cache_dyad_v);
    else lp_degree_0_1 = coef_degree_0_1.elem(cache_dyad_u) + coef_degree_0_1.elem(cache_dyad_v);
  }
  lp_degree_0_1 += baseline_0_1_;
  current_baseline_0_1 = baseline_0_1_;

  lp_degree_1_0 = arma::zeros(n_rows);
  if (coef_degree_1_0.n_elem == 1) {
    lp_degree_1_0.fill(2.0 * coef_degree_1_0.at(0));
  } else if (coef_degree_1_0.n_elem > 1) {
    if (directed) lp_degree_1_0 = coef_degree_1_0.elem(cache_dyad_u) + coef_degree_1_0.elem(n_nodes + cache_dyad_v);
    else lp_degree_1_0 = coef_degree_1_0.elem(cache_dyad_u) + coef_degree_1_0.elem(cache_dyad_v);
  }
  lp_degree_1_0 += baseline_1_0_;
  current_baseline_1_0 = baseline_1_0_;

  // 2. Initialize LPs and Exp LPs
  lp_0_1 = lp_degree_0_1;
  unsigned int n_stats = data.n_cols - 7;

  if (coef_0_1.n_elem > 0) {
    if (n_stats > 0) {
      lp_0_1 += data.cols(7, 6 + n_stats) * coef_0_1;
    } else {
      lp_0_1 += coef_0_1.at(0);
    }
  }
  exp_lp_0_1 = arma::exp(lp_0_1);

  lp_1_0 = lp_degree_1_0;
  if (coef_1_0.n_elem > 0) {
    if (n_stats > 0) {
      lp_1_0 += data.cols(7, 6 + n_stats) * coef_1_0;
    } else {
      lp_1_0 += coef_1_0.at(0);
    }
  }
  exp_lp_1_0 = arma::exp(lp_1_0);

  // 3. Initialize combined_intensity
  combined_intensity = arma::zeros(n_rows);
  arma::uvec all_indices = arma::regspace<arma::uvec>(0, n_rows - 1);
  refresh_combined_intensity(all_indices);
}

void Current_Stat::update_intensities_at_indices(const arma::uvec& indices, const arma::vec& coef_0_1, const arma::vec& coef_1_0) {
  if (indices.is_empty()) return;

  unsigned int n_stats = data.n_cols - 7;

  if (coef_0_1.n_elem > 0) {
    if (n_stats > 0) {
      arma::uvec stat_cols = arma::regspace<arma::uvec>(7, 6 + n_stats);
      lp_0_1.elem(indices) = lp_degree_0_1.elem(indices) + data.submat(indices, stat_cols) * coef_0_1;
    } else {
      lp_0_1.elem(indices) = lp_degree_0_1.elem(indices) + coef_0_1.at(0);
    }
    exp_lp_0_1.elem(indices) = arma::exp(lp_0_1.elem(indices));
  } else {
    lp_0_1.elem(indices) = lp_degree_0_1.elem(indices);
    exp_lp_0_1.elem(indices) = arma::exp(lp_0_1.elem(indices));
  }

  if (coef_1_0.n_elem > 0) {
    if (n_stats > 0) {
      arma::uvec stat_cols = arma::regspace<arma::uvec>(7, 6 + n_stats);
      lp_1_0.elem(indices) = lp_degree_1_0.elem(indices) + data.submat(indices, stat_cols) * coef_1_0;
    } else {
      lp_1_0.elem(indices) = lp_degree_1_0.elem(indices) + coef_1_0.at(0);
    }
    exp_lp_1_0.elem(indices) = arma::exp(lp_1_0.elem(indices));
  } else {
    lp_1_0.elem(indices) = lp_degree_1_0.elem(indices);
    exp_lp_1_0.elem(indices) = arma::exp(lp_1_0.elem(indices));
  }

  refresh_combined_intensity(indices);
}

void Current_Stat::update_baseline(double new_baseline_0_1, double new_baseline_1_0,
                                   const arma::vec& coef_0_1, const arma::vec& coef_1_0) {
  double delta_0_1 = new_baseline_0_1 - current_baseline_0_1;
  double delta_1_0 = new_baseline_1_0 - current_baseline_1_0;

  if (std::abs(delta_0_1) < 1e-15 && std::abs(delta_1_0) < 1e-15) return;

  // Shift the degree+baseline cache vectors
  if (std::abs(delta_0_1) >= 1e-15) {
    lp_degree_0_1 += delta_0_1;
    current_baseline_0_1 = new_baseline_0_1;
  }
  if (std::abs(delta_1_0) >= 1e-15) {
    lp_degree_1_0 += delta_1_0;
    current_baseline_1_0 = new_baseline_1_0;
  }

  // Re-compute lp, exp_lp for all rows
  arma::uvec all_idx = arma::regspace<arma::uvec>(0, data.n_rows - 1);
  update_intensities_at_indices(all_idx, coef_0_1, coef_1_0);
}

void Current_Stat::refresh_combined_intensity(const arma::uvec& indices) {
    if (indices.is_empty()) return;
    if (combined_intensity.is_empty()) return;

    for (unsigned int i = 0; i < indices.n_elem; ++i) {
        unsigned int idx = indices[i];

        if (is_interacting[idx]) {
            combined_intensity[idx] = exp_lp_1_0[idx];
        } else {
            double intensity = exp_lp_0_1[idx];
            if (!simultaneous_interactions) {
                intensity *= (data(idx, 5) * data(idx, 6));
            }
            combined_intensity[idx] = intensity;
        }
    }
}

arma::mat Current_Stat::return_data(){
  return data;
}

void Current_Stat::clear() {
  unsigned int n_statistics = data.n_cols - 7;
  unsigned int n_rows = data.n_rows;

  data.col(1).zeros();
  data.col(2).zeros();
  data.col(5).ones();
  data.col(6).ones();
  if (n_statistics > 0) {
    this->update_exogenous(this->eval_at_zero);
  }
  is_interacting.assign(n_rows, 0);
}

void Current_Stat::set_stat(unsigned int row, unsigned int col, double val) { data.at(row, col) = val; }
