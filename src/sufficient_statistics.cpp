#include <RcppArmadillo.h>
#include "redeem/sufficient_statistics.h"
#include "redeem/extension_api.hpp"

// Determine the change in statistics based on event type
double get_stat_move(unsigned int type, bool is_windowed, int window_size) {
  if (type == 1) return 1.0;
  if (is_windowed && type >= 10 && (int)type == (window_size + 10)) return -1.0;
  return 0.0;
}

// --- Base Statistics ---

arma::uvec stat_dyadic_cov(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  return arma::uvec();
}
TERM_REGISTER("dyadic_cov", stat_dyadic_cov);

arma::uvec stat_intercept(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  return arma::uvec();
}
TERM_REGISTER("Intercept", stat_intercept);
TERM_REGISTER("intercept", stat_intercept);


// --- Common Partner Statistics ---

arma::uvec stat_general_degree_out_sender(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);
  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();

  arma::uvec res = object.current_stats.find_from(from);
  apply_update(object, res, col_number, val, transformation, data);
  return res;
}
TERM_REGISTER("general_degree_out_sender", stat_general_degree_out_sender);
TERM_REGISTER("degree_out_sender", stat_general_degree_out_sender);

arma::uvec stat_general_degree_out_receiver(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);
  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();

  arma::uvec res = object.current_stats.find_to(from);
  apply_update(object, res, col_number, val, transformation, data);
  return res;
}
TERM_REGISTER("general_degree_out_receiver", stat_general_degree_out_receiver);
TERM_REGISTER("degree_out_receiver", stat_general_degree_out_receiver);

arma::uvec stat_general_degree_in_sender(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);
  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();

  arma::uvec res = object.current_stats.find_from(to);
  apply_update(object, res, col_number, val, transformation, data);
  return res;
}
TERM_REGISTER("general_degree_in_sender", stat_general_degree_in_sender);
TERM_REGISTER("degree_in_sender", stat_general_degree_in_sender);

arma::uvec stat_general_degree_in_receiver(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);
  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();

  arma::uvec res = object.current_stats.find_to(to);
  apply_update(object, res, col_number, val, transformation, data);
  return res;
}
TERM_REGISTER("general_degree_in_receiver", stat_general_degree_in_receiver);
TERM_REGISTER("degree_in_receiver", stat_general_degree_in_receiver);

arma::uvec stat_general_degree_sum(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);
  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();

  arma::uvec res_from = object.current_stats.find_involved(from);
  arma::uvec res_to = object.current_stats.find_involved(to);
  apply_update(object, res_from, col_number, val, transformation, data);
  apply_update(object, res_to, col_number, val, transformation, data);
  return arma::join_cols(res_from, res_to);
}
TERM_REGISTER("general_degree_sum", stat_general_degree_sum);
TERM_REGISTER("degree_sum", stat_general_degree_sum);

arma::uvec stat_general_degree_absdiff(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);
  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();

  unsigned int u = from;
  unsigned int v = to;
  unsigned int focal_idx = object.current_stats.find_from_to(u, v);
  unsigned int deg_u = history->get_degree(u);
  unsigned int deg_v = history->get_degree(v);
  std::vector<unsigned int> pos, neg;
  for (unsigned int idx : object.current_stats.find_involved(u)) {
    if (idx == focal_idx) continue;
    unsigned int node_i = (unsigned int)object.current_stats.data.at(idx, 3);
    unsigned int node_j = (unsigned int)object.current_stats.data.at(idx, 4);
    unsigned int other = (node_i == u) ? node_j : node_i;
    unsigned int deg_other = history->get_degree(other);
    if (deg_u >= deg_other) pos.push_back(idx); else neg.push_back(idx);
  }
  for (unsigned int idx : object.current_stats.find_involved(v)) {
    if (idx == focal_idx) continue;
    unsigned int node_i = (unsigned int)object.current_stats.data.at(idx, 3);
    unsigned int node_j = (unsigned int)object.current_stats.data.at(idx, 4);
    unsigned int other = (node_i == v) ? node_j : node_i;
    unsigned int deg_other = history->get_degree(other);
    if (deg_v >= deg_other) pos.push_back(idx); else neg.push_back(idx);
  }
  if (!pos.empty()) apply_update(object, arma::uvec(pos), col_number, val, transformation, data);
  if (!neg.empty()) apply_update(object, arma::uvec(neg), col_number, -val, transformation, data);
  return arma::unique(arma::join_cols(object.current_stats.find_involved(u), object.current_stats.find_involved(v)));
}
TERM_REGISTER("general_degree_absdiff", stat_general_degree_absdiff);
TERM_REGISTER("degree_absdiff", stat_general_degree_absdiff);


// --- Count-based Degree Statistics (General History) ---

arma::uvec stat_general_count_out_sender(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  select_history(object, col_number, is_windowed, window_size);
  double move = get_stat_move(type, is_windowed, window_size);
  if (move == 0.0) return arma::uvec();

  arma::uvec res = object.current_stats.find_from(from);
  apply_update(object, res, col_number, move, transformation, data);
  return res;
}
TERM_REGISTER("general_count_out_sender", stat_general_count_out_sender);

arma::uvec stat_general_count_out_receiver(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  select_history(object, col_number, is_windowed, window_size);
  double move = get_stat_move(type, is_windowed, window_size);
  if (move == 0.0) return arma::uvec();

  arma::uvec res = object.current_stats.find_to(from);
  apply_update(object, res, col_number, move, transformation, data);
  return res;
}
TERM_REGISTER("general_count_out_receiver", stat_general_count_out_receiver);

arma::uvec stat_general_count_in_sender(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  select_history(object, col_number, is_windowed, window_size);
  double move = get_stat_move(type, is_windowed, window_size);
  if (move == 0.0) return arma::uvec();

  arma::uvec res = object.current_stats.find_from(to);
  apply_update(object, res, col_number, move, transformation, data);
  return res;
}
TERM_REGISTER("general_count_in_sender", stat_general_count_in_sender);

arma::uvec stat_general_count_in_receiver(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  select_history(object, col_number, is_windowed, window_size);
  double move = get_stat_move(type, is_windowed, window_size);
  if (move == 0.0) return arma::uvec();

  arma::uvec res = object.current_stats.find_to(to);
  apply_update(object, res, col_number, move, transformation, data);
  return res;
}
TERM_REGISTER("general_count_in_receiver", stat_general_count_in_receiver);

arma::uvec stat_general_count_sum(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  select_history(object, col_number, is_windowed, window_size);
  double move = get_stat_move(type, is_windowed, window_size);
  if (move == 0.0) return arma::uvec();

  arma::uvec res_from = object.current_stats.find_involved(from);
  arma::uvec res_to = object.current_stats.find_involved(to);
  apply_update(object, res_from, col_number, move, transformation, data);
  apply_update(object, res_to, col_number, move, transformation, data);
  return arma::join_cols(res_from, res_to);
}
TERM_REGISTER("general_count_sum", stat_general_count_sum);

arma::uvec stat_general_count_absdiff(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);
  double move = get_stat_move(type, is_windowed, window_size);
  if (move == 0.0) return arma::uvec();
  unsigned int u = from;
  unsigned int v = to;
  unsigned int focal_idx = object.current_stats.find_from_to(u, v);
  unsigned int deg_u = history->get_degree(u, "out", true);
  unsigned int deg_v = history->get_degree(v, "out", true);
  std::vector<unsigned int> pos, neg;
  for (unsigned int idx : object.current_stats.find_involved(u)) {
    if (idx == focal_idx) continue;
    unsigned int node_i = (unsigned int)object.current_stats.data.at(idx, 3);
    unsigned int node_j = (unsigned int)object.current_stats.data.at(idx, 4);
    unsigned int other = (node_i == u) ? node_j : node_i;
    unsigned int deg_other = history->get_degree(other, "out", true);
    if (deg_u >= deg_other) pos.push_back(idx); else neg.push_back(idx);
  }
  for (unsigned int idx : object.current_stats.find_involved(v)) {
    if (idx == focal_idx) continue;
    unsigned int node_i = (unsigned int)object.current_stats.data.at(idx, 3);
    unsigned int node_j = (unsigned int)object.current_stats.data.at(idx, 4);
    unsigned int other = (node_i == v) ? node_j : node_i;
    unsigned int deg_other = history->get_degree(other, "out", true);
    if (deg_v >= deg_other) pos.push_back(idx); else neg.push_back(idx);
  }
  if (!pos.empty()) apply_update(object, arma::uvec(pos), col_number, move, transformation, data);
  if (!neg.empty()) apply_update(object, arma::uvec(neg), col_number, -move, transformation, data);
  return arma::unique(arma::join_cols(object.current_stats.find_involved(u), object.current_stats.find_involved(v)));
}
TERM_REGISTER("general_count_absdiff", stat_general_count_absdiff);


// --- Current Degree/Count Statistics (Active interactions) ---

arma::uvec stat_current_degree_out_sender(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  double val = 0;
  if (type == 1) {
    if (!object.current_interactions.have_interacted(from, to)) val = 1.0;
  } else {
    if (object.current_interactions.get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();
  arma::uvec res = object.current_stats.find_from(from);
  apply_update(object, res, col_number, val, transformation, data);
  return res;
}
TERM_REGISTER("current_degree_out_sender", stat_current_degree_out_sender);

arma::uvec stat_current_degree_out_receiver(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  double val = 0;
  if (type == 1) {
    if (!object.current_interactions.have_interacted(from, to)) val = 1.0;
  } else {
    if (object.current_interactions.get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();
  arma::uvec res = object.current_stats.find_to(from);
  apply_update(object, res, col_number, val, transformation, data);
  return res;
}
TERM_REGISTER("current_degree_out_receiver", stat_current_degree_out_receiver);

arma::uvec stat_current_degree_in_sender(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  double val = 0;
  if (type == 1) {
    if (!object.current_interactions.have_interacted(from, to)) val = 1.0;
  } else {
    if (object.current_interactions.get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();
  arma::uvec res = object.current_stats.find_from(to);
  apply_update(object, res, col_number, val, transformation, data);
  return res;
}
TERM_REGISTER("current_degree_in_sender", stat_current_degree_in_sender);

arma::uvec stat_current_degree_in_receiver(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  double val = 0;
  if (type == 1) {
    if (!object.current_interactions.have_interacted(from, to)) val = 1.0;
  } else {
    if (object.current_interactions.get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();
  arma::uvec res = object.current_stats.find_to(to);
  apply_update(object, res, col_number, val, transformation, data);
  return res;
}
TERM_REGISTER("current_degree_in_receiver", stat_current_degree_in_receiver);

arma::uvec stat_current_count_out_sender(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  arma::uvec res = object.current_stats.find_from(from);
  apply_update(object, res, col_number, type ? 1.0 : -1.0, transformation, data);
  return res;
}
TERM_REGISTER("current_count_out_sender", stat_current_count_out_sender);

arma::uvec stat_current_count_out_receiver(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  arma::uvec res = object.current_stats.find_to(from);
  apply_update(object, res, col_number, type ? 1.0 : -1.0, transformation, data);
  return res;
}
TERM_REGISTER("current_count_out_receiver", stat_current_count_out_receiver);

arma::uvec stat_current_count_in_sender(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  arma::uvec res = object.current_stats.find_from(to);
  apply_update(object, res, col_number, type ? 1.0 : -1.0, transformation, data);
  return res;
}
TERM_REGISTER("current_count_in_sender", stat_current_count_in_sender);

arma::uvec stat_current_count_in_receiver(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  arma::uvec res = object.current_stats.find_to(to);
  apply_update(object, res, col_number, type ? 1.0 : -1.0, transformation, data);
  return res;
}
TERM_REGISTER("current_count_in_receiver", stat_current_count_in_receiver);

arma::uvec stat_current_degree_sum(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  double val = 0;
  if (type == 1) {
    if (!object.current_interactions.have_interacted(from, to)) val = 1.0;
  } else {
    if (object.current_interactions.get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();
  arma::uvec res_from = object.current_stats.find_involved(from);
  arma::uvec res_to = object.current_stats.find_involved(to);
  apply_update(object, res_from, col_number, val, transformation, data);
  apply_update(object, res_to, col_number, val, transformation, data);
  return arma::join_cols(res_from, res_to);
}
TERM_REGISTER("current_degree_sum", stat_current_degree_sum);

arma::uvec stat_current_count_sum(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  double val = type ? 1.0 : -1.0;
  arma::uvec res_from = object.current_stats.find_involved(from);
  arma::uvec res_to = object.current_stats.find_involved(to);
  apply_update(object, res_from, col_number, val, transformation, data);
  apply_update(object, res_to, col_number, val, transformation, data);
  return arma::join_cols(res_from, res_to);
}
TERM_REGISTER("current_count_sum", stat_current_count_sum);

arma::uvec stat_current_degree_absdiff(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  if (type == 1) {
    if (object.current_interactions.have_interacted(from, to)) return arma::uvec();
  } else {
    if (object.current_interactions.get_count(from, to) > 1) return arma::uvec();
  }
  unsigned int u = from; unsigned int v = to;
  unsigned int focal_idx = object.current_stats.find_from_to(u, v);
  unsigned int deg_u = object.current_interactions.get_degree(u, "out", false);
  unsigned int deg_v = object.current_interactions.get_degree(v, "out", false);
  std::vector<unsigned int> pos, neg;
  for (unsigned int idx : object.current_stats.find_involved(u)) {
    if (idx == focal_idx) continue;
    unsigned int node_i = (unsigned int)object.current_stats.data.at(idx, 3);
    unsigned int node_j = (unsigned int)object.current_stats.data.at(idx, 4);
    unsigned int other = (node_i == u) ? node_j : node_i;
    unsigned int deg_other = object.current_interactions.get_degree(other, "out", false);
    if (type ? (deg_u >= deg_other) : (deg_u > deg_other)) pos.push_back(idx); else neg.push_back(idx);
  }
  for (unsigned int idx : object.current_stats.find_involved(v)) {
    if (idx == focal_idx) continue;
    unsigned int node_i = (unsigned int)object.current_stats.data.at(idx, 3);
    unsigned int node_j = (unsigned int)object.current_stats.data.at(idx, 4);
    unsigned int other = (node_i == v) ? node_j : node_i;
    unsigned int deg_other = object.current_interactions.get_degree(other, "out", false);
    if (type ? (deg_v >= deg_other) : (deg_v > deg_other)) pos.push_back(idx); else neg.push_back(idx);
  }
  double delta = type ? 1.0 : -1.0;
  if (!pos.empty()) apply_update(object, arma::uvec(pos), col_number, delta, transformation, data);
  if (!neg.empty()) apply_update(object, arma::uvec(neg), col_number, -delta, transformation, data);
  return arma::unique(arma::join_cols(object.current_stats.find_involved(u), object.current_stats.find_involved(v)));
}
TERM_REGISTER("current_degree_absdiff", stat_current_degree_absdiff);

arma::uvec stat_current_count_absdiff(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  unsigned int u = from; unsigned int v = to;
  unsigned int focal_idx = object.current_stats.find_from_to(u, v);
  unsigned int deg_u = object.current_interactions.get_degree(u, "out", true);
  unsigned int deg_v = object.current_interactions.get_degree(v, "out", true);
  std::vector<unsigned int> pos, neg;
  for (unsigned int idx : object.current_stats.find_involved(u)) {
    if (idx == focal_idx) continue;
    unsigned int other = ((unsigned int)object.current_stats.data.at(idx, 3) == u) ? (unsigned int)object.current_stats.data.at(idx, 4) : (unsigned int)object.current_stats.data.at(idx, 3);
    unsigned int deg_other = object.current_interactions.get_degree(other, "out", true);
    if (type ? (deg_u >= deg_other) : (deg_u > deg_other)) pos.push_back(idx); else neg.push_back(idx);
  }
  for (unsigned int idx : object.current_stats.find_involved(v)) {
    if (idx == focal_idx) continue;
    unsigned int other = ((unsigned int)object.current_stats.data.at(idx, 3) == v) ? (unsigned int)object.current_stats.data.at(idx, 4) : (unsigned int)object.current_stats.data.at(idx, 3);
    unsigned int deg_other = object.current_interactions.get_degree(other, "out", true);
    if (type ? (deg_v >= deg_other) : (deg_v > deg_other)) pos.push_back(idx); else neg.push_back(idx);
  }
  if (!pos.empty()) apply_update(object, arma::uvec(pos), col_number, type ? 1.0 : -1.0, transformation, data);
  if (!neg.empty()) apply_update(object, arma::uvec(neg), col_number, type ? -1.0 : 1.0, transformation, data);
  return arma::unique(arma::join_cols(object.current_stats.find_involved(u), object.current_stats.find_involved(v)));
}
TERM_REGISTER("current_count_absdiff", stat_current_count_absdiff);


// --- Triadic / Triangle Statistics ---

// --- Undirected Common Partner ---

arma::uvec stat_current_common_partner(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  double move = (type) ? 1.0 : -1.0;

  arma::uvec ps_from = object.current_interactions.get_partners(from, "out");
  arma::uvec ps_to = object.current_interactions.get_partners(to, "out");

  arma::uvec js_to_update(ps_from.n_elem + ps_to.n_elem);
  unsigned int count = 0;
  for (unsigned int j : ps_from) {
    if (j != from && j != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(to, j);
    }
  }
  for (unsigned int j : ps_to) {
    if (j != from && j != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(from, j);
    }
  }

  if (count > 0) {
    arma::uvec res = js_to_update.head(count);
    apply_update(object, res, col_number, move, transformation, data);
    return res;
  }

  return arma::uvec();
}
TERM_REGISTER("current_common_partner", stat_current_common_partner);

arma::uvec stat_general_common_partner(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);
  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();

  arma::uvec ps_from = history->get_partners(from, "out");
  arma::uvec ps_to = history->get_partners(to, "out");

  arma::uvec js_to_update(ps_from.n_elem + ps_to.n_elem);
  unsigned int count = 0;
  for (unsigned int j : ps_from) {
    if (j != from && j != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(to, j);
    }
  }
  for (unsigned int j : ps_to) {
    if (j != from && j != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(from, j);
    }
  }

  if (count > 0) {
    arma::uvec res = arma::unique(js_to_update.head(count));
    apply_update(object, res, col_number, val, transformation, data);
    return res;
  }

  return arma::uvec();
}
TERM_REGISTER("general_common_partner", stat_general_common_partner);


// --- Directed OSP (Outgoing Shared Partner) Common Partners ---

arma::uvec stat_current_common_partner_OSP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  double move = (type) ? 1.0 : -1.0;

  arma::uvec ps_to = object.current_interactions.get_partners(to, "in");

  arma::uvec js_to_update(ps_to.n_elem * 2);
  unsigned int count = 0;
  for (unsigned int j : ps_to) {
    if (j != from && j != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(from, j);
      js_to_update.at(count++) = object.current_stats.find_from_to(j, from);
    }
  }

  if (count > 0) {
    arma::uvec unique_indices = arma::unique(js_to_update.head(count));
    apply_update(object, unique_indices, col_number, move, transformation, data);
    return unique_indices;
  }

  return arma::uvec();
}
TERM_REGISTER("current_common_partner_OSP", stat_current_common_partner_OSP);

arma::uvec stat_general_common_partner_OSP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);
  
  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();

  arma::uvec ps_to = history->get_partners(to, "in");

  arma::uvec js_to_update(ps_to.n_elem * 2);
  unsigned int count = 0;

  for (unsigned int j : ps_to) {
    if (j != from && j != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(from, j);
      js_to_update.at(count++) = object.current_stats.find_from_to(j, from);
    }
  }

  if (count > 0) {
    arma::uvec res = arma::unique(js_to_update.head(count));
    apply_update(object, res, col_number, val, transformation, data);
    return res;
  }

  return arma::uvec();
}
TERM_REGISTER("general_common_partner_OSP", stat_general_common_partner_OSP);


// --- Directed OSP (Outgoing Shared Partner) Triangles ---

arma::uvec stat_current_triangle_OSP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  std::vector<unsigned int> changed_indices;

  unsigned int focal_idx = object.current_stats.find_from_to(from, to);
  unsigned int cp_focal = object.current_interactions.get_common_partners(from, to, "OSP").n_elem;

  if (type == 1) {
    if (cp_focal > 0) {
      apply_update(object, arma::uvec{focal_idx}, col_number, (double)cp_focal, transformation, data);
      changed_indices.push_back(focal_idx);
    }
  } else {
    if (cp_focal > 0) {
      apply_update(object, arma::uvec{focal_idx}, col_number, -(double)cp_focal, transformation, data);
      changed_indices.push_back(focal_idx);
    }
  }

  double move = (type == 1) ? 1.0 : -1.0;
  arma::uvec ps_to = object.current_interactions.get_partners(to, "in");

  for (unsigned int j : ps_to) {
    if (j != from && j != to) {
      if (object.current_interactions.have_interacted(from, j)) {
        unsigned int idx = object.current_stats.find_from_to(from, j);
        apply_update(object, arma::uvec{idx}, col_number, move, transformation, data);
        changed_indices.push_back(idx);
      }
      if (object.current_interactions.have_interacted(j, from)) {
        unsigned int idx = object.current_stats.find_from_to(j, from);
        apply_update(object, arma::uvec{idx}, col_number, move, transformation, data);
        changed_indices.push_back(idx);
      }
    }
  }

  if (changed_indices.empty()) return arma::uvec();
  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("current_triangle_OSP", stat_current_triangle_OSP);

arma::uvec stat_general_triangle_OSP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);

  std::vector<unsigned int> changed_indices;
  unsigned int focal_idx = object.current_stats.find_from_to(from, to);
  unsigned int cp_focal = history->get_common_partners(from, to, "OSP").n_elem;

  if (type == 1) {
    if (!history->have_interacted(from, to)) {
      if (cp_focal > 0) {
        apply_update(object, arma::uvec{focal_idx}, col_number, (double)cp_focal, transformation, data);
        changed_indices.push_back(focal_idx);
      }
    }
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) {
      if (cp_focal > 0) {
        apply_update(object, arma::uvec{focal_idx}, col_number, -(double)cp_focal, transformation, data);
        changed_indices.push_back(focal_idx);
      }
    }
  }

  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }

  if (val != 0) {
    arma::uvec ps_to = history->get_partners(to, "in");
    for (unsigned int j : ps_to) {
      if (j != from && j != to) {
        if (history->have_interacted(from, j)) {
          unsigned int idx = object.current_stats.find_from_to(from, j);
          apply_update(object, arma::uvec{idx}, col_number, val, transformation, data);
          changed_indices.push_back(idx);
        }
        if (history->have_interacted(j, from)) {
          unsigned int idx = object.current_stats.find_from_to(j, from);
          apply_update(object, arma::uvec{idx}, col_number, val, transformation, data);
          changed_indices.push_back(idx);
        }
      }
    }
  }

  if (changed_indices.empty()) return arma::uvec();
  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("general_triangle_OSP", stat_general_triangle_OSP);


// --- Directed ISP (Incoming Shared Partner) Common Partners ---

arma::uvec stat_current_common_partner_ISP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  double move = (type) ? 1.0 : -1.0;

  arma::uvec js_out_from = object.current_interactions.get_partners(from, "out");
  
  arma::uvec js_to_update(js_out_from.n_elem * 2);
  unsigned int count = 0;
  for (unsigned int j : js_out_from) {
    if (j != from && j != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(to, j);
      js_to_update.at(count++) = object.current_stats.find_from_to(j, to);
    }
  }

  if (count > 0) {
    arma::uvec unique_indices = arma::unique(js_to_update.head(count));
    apply_update(object, unique_indices, col_number, move, transformation, data);
    return unique_indices;
  }

  return arma::uvec();
}
TERM_REGISTER("current_common_partner_ISP", stat_current_common_partner_ISP);

arma::uvec stat_general_common_partner_ISP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);
  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();

  arma::uvec js_out_from = history->get_partners(from, "out");

  arma::uvec js_to_update(js_out_from.n_elem * 2);
  unsigned int count = 0;
  for (unsigned int j : js_out_from) {
    if (j != from && j != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(to, j);
      js_to_update.at(count++) = object.current_stats.find_from_to(j, to);
    }
  }

  if (count > 0) {
    arma::uvec res = arma::unique(js_to_update.head(count));
    apply_update(object, res, col_number, val, transformation, data);
    return res;
  }

  return arma::uvec();
}
TERM_REGISTER("general_common_partner_ISP", stat_general_common_partner_ISP);


// --- Directed ISP (Incoming Shared Partner) Triangles ---

arma::uvec stat_current_triangle_ISP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  std::vector<unsigned int> changed_indices;

  unsigned int focal_idx = object.current_stats.find_from_to(from, to);
  unsigned int cp_focal = object.current_interactions.get_common_partners(from, to, "ISP").n_elem;

  if (type == 1) {
    if (cp_focal > 0) {
      apply_update(object, arma::uvec{focal_idx}, col_number, (double)cp_focal, transformation, data);
      changed_indices.push_back(focal_idx);
    }
  } else {
    if (cp_focal > 0) {
      apply_update(object, arma::uvec{focal_idx}, col_number, -(double)cp_focal, transformation, data);
      changed_indices.push_back(focal_idx);
    }
  }

  double move = (type == 1) ? 1.0 : -1.0;
  arma::uvec js_out_from = object.current_interactions.get_partners(from, "out");

  for (unsigned int j : js_out_from) {
    if (j != from && j != to) {
      if (object.current_interactions.have_interacted(to, j)) {
        unsigned int idx = object.current_stats.find_from_to(to, j);
        apply_update(object, arma::uvec{idx}, col_number, move, transformation, data);
        changed_indices.push_back(idx);
      }
      if (object.current_interactions.have_interacted(j, to)) {
        unsigned int idx = object.current_stats.find_from_to(j, to);
        apply_update(object, arma::uvec{idx}, col_number, move, transformation, data);
        changed_indices.push_back(idx);
      }
    }
  }

  if (changed_indices.empty()) return arma::uvec();
  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("current_triangle_ISP", stat_current_triangle_ISP);

arma::uvec stat_general_triangle_ISP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);

  std::vector<unsigned int> changed_indices;
  unsigned int focal_idx = object.current_stats.find_from_to(from, to);
  unsigned int cp_focal = history->get_common_partners(from, to, "ISP").n_elem;

  if (type == 1) {
    if (!history->have_interacted(from, to)) {
      if (cp_focal > 0) {
        apply_update(object, arma::uvec{focal_idx}, col_number, (double)cp_focal, transformation, data);
        changed_indices.push_back(focal_idx);
      }
    }
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) {
      if (cp_focal > 0) {
        apply_update(object, arma::uvec{focal_idx}, col_number, -(double)cp_focal, transformation, data);
        changed_indices.push_back(focal_idx);
      }
    }
  }

  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }

  if (val != 0) {
    arma::uvec js_out_from = history->get_partners(from, "out");
    for (unsigned int j : js_out_from) {
      if (j != from && j != to) {
        if (history->have_interacted(to, j)) {
          unsigned int idx = object.current_stats.find_from_to(to, j);
          apply_update(object, arma::uvec{idx}, col_number, val, transformation, data);
          changed_indices.push_back(idx);
        }
        if (history->have_interacted(j, to)) {
          unsigned int idx = object.current_stats.find_from_to(j, to);
          apply_update(object, arma::uvec{idx}, col_number, val, transformation, data);
          changed_indices.push_back(idx);
        }
      }
    }
  }

  if (changed_indices.empty()) return arma::uvec();
  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("general_triangle_ISP", stat_general_triangle_ISP);


// --- Directed OTP (Outgoing Two-Path) Common Partners ---

arma::uvec stat_current_common_partner_OTP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  double move = (type) ? 1.0 : -1.0;

  arma::uvec js_out_to = object.current_interactions.get_partners(to, "out");
  arma::uvec is_in_from = object.current_interactions.get_partners(from, "in");

  arma::uvec js_to_update(js_out_to.n_elem + is_in_from.n_elem);
  unsigned int count = 0;
  for (unsigned int j : js_out_to) {
    if (j != from && j != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(from, j);
    }
  }
  for (unsigned int i : is_in_from) {
    if (i != from && i != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(i, to);
    }
  }
  if (count > 0) {
    arma::uvec unique_indices = arma::unique(js_to_update.head(count));
    apply_update(object, unique_indices, col_number, move, transformation, data);
    return unique_indices;
  }

  return arma::uvec();
}
TERM_REGISTER("current_common_partner_OTP", stat_current_common_partner_OTP);

arma::uvec stat_general_common_partner_OTP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);
  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();

  arma::uvec js_out_to = history->get_partners(to, "out");
  arma::uvec is_in_from = history->get_partners(from, "in");

  arma::uvec js_to_update(js_out_to.n_elem + is_in_from.n_elem);
  unsigned int count = 0;
  for (unsigned int j : js_out_to) {
    if (j != from && j != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(from, j);
    }
  }
  for (unsigned int i : is_in_from) {
    if (i != from && i != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(i, to);
    }
  }
  if (count > 0) {
    arma::uvec unique_indices = arma::unique(js_to_update.head(count));
    apply_update(object, unique_indices, col_number, val, transformation, data);
    return unique_indices;
  }

  return arma::uvec();
}
TERM_REGISTER("general_common_partner_OTP", stat_general_common_partner_OTP);


// --- Directed OTP (Outgoing Two-Path) Triangles ---

arma::uvec stat_current_triangle_OTP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  std::vector<unsigned int> changed_indices;

  unsigned int focal_idx = object.current_stats.find_from_to(from, to);
  unsigned int cp_focal = object.current_interactions.get_common_partners(from, to, "OTP").n_elem;

  if (type == 1) {
    if (cp_focal > 0) {
      apply_update(object, arma::uvec{focal_idx}, col_number, (double)cp_focal, transformation, data);
      changed_indices.push_back(focal_idx);
    }
  } else {
    if (cp_focal > 0) {
      apply_update(object, arma::uvec{focal_idx}, col_number, -(double)cp_focal, transformation, data);
      changed_indices.push_back(focal_idx);
    }
  }

  double move = (type == 1) ? 1.0 : -1.0;
  arma::uvec js_out_to = object.current_interactions.get_partners(to, "out");
  arma::uvec is_in_from = object.current_interactions.get_partners(from, "in");

  for (unsigned int j : js_out_to) {
    if (j != from && j != to) {
      if (object.current_interactions.have_interacted(from, j)) {
        unsigned int idx = object.current_stats.find_from_to(from, j);
        apply_update(object, arma::uvec{idx}, col_number, move, transformation, data);
        changed_indices.push_back(idx);
      }
    }
  }
  for (unsigned int i : is_in_from) {
    if (i != from && i != to) {
      if (object.current_interactions.have_interacted(i, to)) {
        unsigned int idx = object.current_stats.find_from_to(i, to);
        apply_update(object, arma::uvec{idx}, col_number, move, transformation, data);
        changed_indices.push_back(idx);
      }
    }
  }

  if (changed_indices.empty()) return arma::uvec();
  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("current_triangle_OTP", stat_current_triangle_OTP);

arma::uvec stat_general_triangle_OTP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);

  std::vector<unsigned int> changed_indices;
  unsigned int focal_idx = object.current_stats.find_from_to(from, to);
  unsigned int cp_focal = history->get_common_partners(from, to, "OTP").n_elem;

  if (type == 1) {
    if (!history->have_interacted(from, to)) {
      if (cp_focal > 0) {
        apply_update(object, arma::uvec{focal_idx}, col_number, (double)cp_focal, transformation, data);
        changed_indices.push_back(focal_idx);
      }
    }
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) {
      if (cp_focal > 0) {
        apply_update(object, arma::uvec{focal_idx}, col_number, -(double)cp_focal, transformation, data);
        changed_indices.push_back(focal_idx);
      }
    }
  }

  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }

  if (val != 0) {
    arma::uvec js_out_to = history->get_partners(to, "out");
    arma::uvec is_in_from = history->get_partners(from, "in");

    for (unsigned int j : js_out_to) {
      if (j != from && j != to) {
        if (history->have_interacted(from, j)) {
          unsigned int idx = object.current_stats.find_from_to(from, j);
          apply_update(object, arma::uvec{idx}, col_number, val, transformation, data);
          changed_indices.push_back(idx);
        }
      }
    }
    for (unsigned int i : is_in_from) {
      if (i != from && i != to) {
        if (history->have_interacted(i, to)) {
          unsigned int idx = object.current_stats.find_from_to(i, to);
          apply_update(object, arma::uvec{idx}, col_number, val, transformation, data);
          changed_indices.push_back(idx);
        }
      }
    }
  }

  if (changed_indices.empty()) return arma::uvec();
  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("general_triangle_OTP", stat_general_triangle_OTP);


// --- Directed ITP (Incoming Two-Path) Common Partners ---

arma::uvec stat_current_common_partner_ITP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  double move = (type) ? 1.0 : -1.0;

  arma::uvec out_to_edge = object.current_interactions.get_partners(to, "out");
  arma::uvec in_from_edge = object.current_interactions.get_partners(from, "in");

  arma::uvec js_to_update(out_to_edge.n_elem + in_from_edge.n_elem);
  unsigned int count = 0;
  for (unsigned int i : out_to_edge) {
    if (i != from && i != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(i, from);
    }
  }
  for (unsigned int j : in_from_edge) {
    if (j != from && j != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(to, j);
    }
  }
  if (count > 0) {
    arma::uvec unique_indices = arma::unique(js_to_update.head(count));
    apply_update(object, unique_indices, col_number, move, transformation, data);
    return unique_indices;
  }

  return arma::uvec();
}
TERM_REGISTER("current_common_partner_ITP", stat_current_common_partner_ITP);

arma::uvec stat_general_common_partner_ITP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);
  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }
  if (val == 0) return arma::uvec();

  arma::uvec is_out_to = history->get_partners(to, "out");
  arma::uvec js_in_from = history->get_partners(from, "in");

  arma::uvec js_to_update(is_out_to.n_elem + js_in_from.n_elem);
  unsigned int count = 0;
  for (unsigned int i : is_out_to) {
    if (i != from && i != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(i, from);
    }
  }
  for (unsigned int j : js_in_from) {
    if (j != from && j != to) {
      js_to_update.at(count++) = object.current_stats.find_from_to(to, j);
    }
  }
  if (count > 0) {
    arma::uvec unique_indices = arma::unique(js_to_update.head(count));
    apply_update(object, unique_indices, col_number, val, transformation, data);
    return unique_indices;
  }

  return arma::uvec();
}
TERM_REGISTER("general_common_partner_ITP", stat_general_common_partner_ITP);


// --- Directed ITP (Incoming Two-Path) Triangles ---

arma::uvec stat_current_triangle_ITP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  std::vector<unsigned int> changed_indices;

  unsigned int focal_idx = object.current_stats.find_from_to(from, to);
  unsigned int cp_focal = object.current_interactions.get_common_partners(from, to, "ITP").n_elem;

  if (type == 1) {
    if (cp_focal > 0) {
      apply_update(object, arma::uvec{focal_idx}, col_number, (double)cp_focal, transformation, data);
      changed_indices.push_back(focal_idx);
    }
  } else {
    if (cp_focal > 0) {
      apply_update(object, arma::uvec{focal_idx}, col_number, -(double)cp_focal, transformation, data);
      changed_indices.push_back(focal_idx);
    }
  }

  double move = (type == 1) ? 1.0 : -1.0;
  arma::uvec out_to_edge = object.current_interactions.get_partners(to, "out");
  arma::uvec in_from_edge = object.current_interactions.get_partners(from, "in");

  for (unsigned int i : out_to_edge) {
    if (i != from && i != to) {
      if (object.current_interactions.have_interacted(i, from)) {
        unsigned int idx = object.current_stats.find_from_to(i, from);
        apply_update(object, arma::uvec{idx}, col_number, move, transformation, data);
        changed_indices.push_back(idx);
      }
    }
  }
  for (unsigned int j : in_from_edge) {
    if (j != from && j != to) {
      if (object.current_interactions.have_interacted(to, j)) {
        unsigned int idx = object.current_stats.find_from_to(to, j);
        apply_update(object, arma::uvec{idx}, col_number, move, transformation, data);
        changed_indices.push_back(idx);
      }
    }
  }

  if (changed_indices.empty()) return arma::uvec();
  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("current_triangle_ITP", stat_current_triangle_ITP);

arma::uvec stat_general_triangle_ITP(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);

  std::vector<unsigned int> changed_indices;
  unsigned int focal_idx = object.current_stats.find_from_to(from, to);
  unsigned int cp_focal = history->get_common_partners(from, to, "ITP").n_elem;

  if (type == 1) {
    if (!history->have_interacted(from, to)) {
      if (cp_focal > 0) {
        apply_update(object, arma::uvec{focal_idx}, col_number, (double)cp_focal, transformation, data);
        changed_indices.push_back(focal_idx);
      }
    }
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) {
      if (cp_focal > 0) {
        apply_update(object, arma::uvec{focal_idx}, col_number, -(double)cp_focal, transformation, data);
        changed_indices.push_back(focal_idx);
      }
    }
  }

  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }

  if (val != 0) {
    arma::uvec out_to_edge = history->get_partners(to, "out");
    arma::uvec in_from_edge = history->get_partners(from, "in");

    for (unsigned int i : out_to_edge) {
      if (i != from && i != to) {
        if (history->have_interacted(i, from)) {
          unsigned int idx = object.current_stats.find_from_to(i, from);
          apply_update(object, arma::uvec{idx}, col_number, val, transformation, data);
          changed_indices.push_back(idx);
        }
      }
    }
    for (unsigned int j : in_from_edge) {
      if (j != from && j != to) {
        if (history->have_interacted(to, j)) {
          unsigned int idx = object.current_stats.find_from_to(to, j);
          apply_update(object, arma::uvec{idx}, col_number, val, transformation, data);
          changed_indices.push_back(idx);
        }
      }
    }
  }

  if (changed_indices.empty()) return arma::uvec();
  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("general_triangle_ITP", stat_general_triangle_ITP);


// --- Undirected Triangle ---

arma::uvec stat_current_triangle(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  std::vector<unsigned int> changed_indices;

  unsigned int focal_idx = object.current_stats.find_from_to(from, to);
  unsigned int cp_focal = object.current_interactions.get_common_partners(from, to, "OSP").n_elem;

  if (type == 1) {
    if (cp_focal > 0) {
      apply_update(object, arma::uvec{focal_idx}, col_number, (double)cp_focal, transformation, data);
      changed_indices.push_back(focal_idx);
    }
  } else {
    if (cp_focal > 0) {
      apply_update(object, arma::uvec{focal_idx}, col_number, -(double)cp_focal, transformation, data);
      changed_indices.push_back(focal_idx);
    }
  }

  double move = (type == 1) ? 1.0 : -1.0;
  arma::uvec ps_from = object.current_interactions.get_partners(from, "out");
  arma::uvec ps_to = object.current_interactions.get_partners(to, "out");

  for (unsigned int j : ps_from) {
    if (j != from && j != to) {
      if (object.current_interactions.have_interacted(to, j)) {
        unsigned int idx = object.current_stats.find_from_to(to, j);
        apply_update(object, arma::uvec{idx}, col_number, move, transformation, data);
        changed_indices.push_back(idx);
      }
    }
  }
  for (unsigned int j : ps_to) {
    if (j != from && j != to) {
      if (object.current_interactions.have_interacted(from, j)) {
        unsigned int idx = object.current_stats.find_from_to(from, j);
        apply_update(object, arma::uvec{idx}, col_number, move, transformation, data);
        changed_indices.push_back(idx);
      }
    }
  }

  if (changed_indices.empty()) return arma::uvec();
  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("current_triangle", stat_current_triangle);

arma::uvec stat_general_triangle(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  Hist_Events* history = select_history(object, col_number, is_windowed, window_size);

  std::vector<unsigned int> changed_indices;
  unsigned int focal_idx = object.current_stats.find_from_to(from, to);
  unsigned int cp_focal = history->get_common_partners(from, to, "OSP").n_elem;

  if (type == 1) {
    if (!history->have_interacted(from, to)) {
      if (cp_focal > 0) {
        apply_update(object, arma::uvec{focal_idx}, col_number, (double)cp_focal, transformation, data);
        changed_indices.push_back(focal_idx);
      }
    }
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) {
      if (cp_focal > 0) {
        apply_update(object, arma::uvec{focal_idx}, col_number, -(double)cp_focal, transformation, data);
        changed_indices.push_back(focal_idx);
      }
    }
  }

  double val = 0;
  if (type == 1) {
    if (!history->have_interacted(from, to)) val = 1.0;
  } else if (is_windowed && type >= 10 && (int)type == (window_size + 10)) {
    if (history->get_count(from, to) == 1) val = -1.0;
  }

  if (val != 0) {
    arma::uvec ps_from = history->get_partners(from, "out");
    arma::uvec ps_to = history->get_partners(to, "out");

    for (unsigned int j : ps_from) {
      if (j != from && j != to) {
        if (history->have_interacted(to, j)) {
          unsigned int idx = object.current_stats.find_from_to(to, j);
          apply_update(object, arma::uvec{idx}, col_number, val, transformation, data);
          changed_indices.push_back(idx);
        }
      }
    }
    for (unsigned int j : ps_to) {
      if (j != from && j != to) {
        if (history->have_interacted(from, j)) {
          unsigned int idx = object.current_stats.find_from_to(from, j);
          apply_update(object, arma::uvec{idx}, col_number, val, transformation, data);
          changed_indices.push_back(idx);
        }
      }
    }
  }

  if (changed_indices.empty()) return arma::uvec();
  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("general_triangle", stat_general_triangle);



// --- Interaction / Structural Statistics ---

arma::uvec stat_current_interaction(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0) {
    arma::uvec interacting = object.current_stats.get_currently_interacting();
    double change = 0;
    if (number_event < object.changepoints.size()) {
        if (number_event == 0) {
            change = object.changepoints.at(0);
        } else {
            change = object.changepoints.at(number_event) - object.changepoints.at(number_event - 1);
        }
    }
    if (change > 0 && !interacting.is_empty()) apply_update(object, interacting, col_number, change, transformation, data);
    return interacting;
  }
  if (to == 0) return arma::uvec();
  unsigned int idx = object.current_stats.find_from_to(from, to);
  double neutral = (transformation == "recip" || transformation == "sig") ? 1.0 : 0.0;
  object.current_stats.set_stat(idx, col_number, neutral);
  return arma::uvec{idx};
}
TERM_REGISTER("current_interaction", stat_current_interaction);

arma::uvec stat_number_interaction(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  select_history(object, col_number, is_windowed, window_size);
  double move = get_stat_move(type, is_windowed, window_size);
  if (object.is_dem && !is_windowed) {
    if (type == 1) move = 0.0;
    if (type == 0) move = 1.0;
  }
  if (move == 0.0) return arma::uvec();

  unsigned int idx = object.current_stats.find_from_to(from, to);
  apply_update(object, arma::uvec{idx}, col_number, move, transformation, data);
  return arma::uvec{idx};
}
TERM_REGISTER("number_interaction", stat_number_interaction);
TERM_REGISTER("inertia", stat_number_interaction);

arma::uvec stat_reciprocity(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (from == 0 || to == 0) return arma::uvec();
  bool is_windowed = false;
  int window_size = -1;
  select_history(object, col_number, is_windowed, window_size);
  double move = get_stat_move(type, is_windowed, window_size);
  if (object.is_dem && !is_windowed) {
    if (type == 1) move = 0.0;
    if (type == 0) move = 1.0;
  }
  if (move == 0.0) return arma::uvec();

  unsigned int idx = object.current_stats.find_from_to(to, from);
  apply_update(object, arma::uvec{idx}, col_number, move, transformation, data);
  return arma::uvec{idx};
}
TERM_REGISTER("reciprocity", stat_reciprocity);

// --- Participation Shifts (P-shifts) ---
// Participation shifts (Gibson, 2003) capture the sequential dependencies between 
// consecutive events. Assuming the previous event was A -> B:

// 1. psABBA: Reciprocation (AB -> BA)
// Measures B's tendency to immediately respond to A. The active dyad is B -> A.
arma::uvec stat_psABBA(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (type >= 10 || type == 3) return arma::uvec();

  unsigned int A = object.last_sender;
  unsigned int B = object.last_receiver;

  std::vector<unsigned int> changed_indices;

  // Reset old active dyad (B -> A)
  if (A != 0 && B != 0) {
    unsigned int old_idx = object.current_stats.find_from_to(B, A);
    object.current_stats.set_stat(old_idx, col_number, 0.0);
    changed_indices.push_back(old_idx);
  }

  // Set new active dyad (to -> from)
  if (from != 0 && to != 0) {
    unsigned int new_idx = object.current_stats.find_from_to(to, from);
    object.current_stats.set_stat(new_idx, col_number, 1.0);
    changed_indices.push_back(new_idx);
  }

  return arma::uvec(changed_indices);
}
TERM_REGISTER("psABBA", stat_psABBA);

// 2. psABBY: Turn-continuing / sender-to-other (AB -> BY)
// B takes the turn and addresses a new actor Y (Y != A, B).
arma::uvec stat_psABBY(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (type >= 10 || type == 3) return arma::uvec();

  unsigned int A = object.last_sender;
  unsigned int B = object.last_receiver;

  std::vector<unsigned int> changed_indices;

  // Reset old active dyads (B -> Y, Y != A)
  if (A != 0 && B != 0) {
    arma::uvec old_dyads = object.current_stats.find_from(B);
    for (unsigned int idx : old_dyads) {
      unsigned int receiver = (unsigned int)object.current_stats.data.at(idx, 4);
      if (receiver != A) {
        object.current_stats.set_stat(idx, col_number, 0.0);
        changed_indices.push_back(idx);
      }
    }
  }

  // Set new active dyads (to -> Y, Y != from)
  if (from != 0 && to != 0) {
    arma::uvec new_dyads = object.current_stats.find_from(to);
    for (unsigned int idx : new_dyads) {
      unsigned int receiver = (unsigned int)object.current_stats.data.at(idx, 4);
      if (receiver != from) {
        object.current_stats.set_stat(idx, col_number, 1.0);
        changed_indices.push_back(idx);
      }
    }
  }

  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("psABBY", stat_psABBY);

// 3. psABAY: Turn-continuing / receiver-to-other (AB -> AY)
// A keeps the turn and addresses a new actor Y (Y != A, B).
arma::uvec stat_psABAY(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (type >= 10 || type == 3) return arma::uvec();

  unsigned int A = object.last_sender;
  unsigned int B = object.last_receiver;

  std::vector<unsigned int> changed_indices;

  // Reset old active dyads (A -> Y, Y != B)
  if (A != 0 && B != 0) {
    arma::uvec old_dyads = object.current_stats.find_from(A);
    for (unsigned int idx : old_dyads) {
      unsigned int receiver = (unsigned int)object.current_stats.data.at(idx, 4);
      if (receiver != B) {
        object.current_stats.set_stat(idx, col_number, 0.0);
        changed_indices.push_back(idx);
      }
    }
  }

  // Set new active dyads (from -> Y, Y != to)
  if (from != 0 && to != 0) {
    arma::uvec new_dyads = object.current_stats.find_from(from);
    for (unsigned int idx : new_dyads) {
      unsigned int receiver = (unsigned int)object.current_stats.data.at(idx, 4);
      if (receiver != to) {
        object.current_stats.set_stat(idx, col_number, 1.0);
        changed_indices.push_back(idx);
      }
    }
  }

  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("psABAY", stat_psABAY);

// 4. psABXA: Turn-usurping / other-to-sender (AB -> XA)
// A new actor X (X != A, B) addresses the original sender A.
arma::uvec stat_psABXA(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (type >= 10 || type == 3) return arma::uvec();

  unsigned int A = object.last_sender;
  unsigned int B = object.last_receiver;

  std::vector<unsigned int> changed_indices;

  // Reset old active dyads (X -> A, X != B)
  if (A != 0 && B != 0) {
    arma::uvec old_dyads = object.current_stats.find_to(A);
    for (unsigned int idx : old_dyads) {
      unsigned int sender = (unsigned int)object.current_stats.data.at(idx, 3);
      if (sender != B) {
        object.current_stats.set_stat(idx, col_number, 0.0);
        changed_indices.push_back(idx);
      }
    }
  }

  // Set new active dyads (X -> from, X != to)
  if (from != 0 && to != 0) {
    arma::uvec new_dyads = object.current_stats.find_to(from);
    for (unsigned int idx : new_dyads) {
      unsigned int sender = (unsigned int)object.current_stats.data.at(idx, 3);
      if (sender != to) {
        object.current_stats.set_stat(idx, col_number, 1.0);
        changed_indices.push_back(idx);
      }
    }
  }

  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("psABXA", stat_psABXA);

// 5. psABXB: Turn-usurping / other-to-receiver (AB -> XB)
// A new actor X (X != A, B) addresses the original receiver B.
arma::uvec stat_psABXB(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (type >= 10 || type == 3) return arma::uvec();

  unsigned int A = object.last_sender;
  unsigned int B = object.last_receiver;

  std::vector<unsigned int> changed_indices;

  // Reset old active dyads (X -> B, X != A)
  if (A != 0 && B != 0) {
    arma::uvec old_dyads = object.current_stats.find_to(B);
    for (unsigned int idx : old_dyads) {
      unsigned int sender = (unsigned int)object.current_stats.data.at(idx, 3);
      if (sender != A) {
        object.current_stats.set_stat(idx, col_number, 0.0);
        changed_indices.push_back(idx);
      }
    }
  }

  // Set new active dyads (X -> to, X != from)
  if (from != 0 && to != 0) {
    arma::uvec new_dyads = object.current_stats.find_to(to);
    for (unsigned int idx : new_dyads) {
      unsigned int sender = (unsigned int)object.current_stats.data.at(idx, 3);
      if (sender != from) {
        object.current_stats.set_stat(idx, col_number, 1.0);
        changed_indices.push_back(idx);
      }
    }
  }

  return arma::unique(arma::uvec(changed_indices));
}
TERM_REGISTER("psABXB", stat_psABXB);

// 6. psABXY: Full shift / usurpation (AB -> XY)
// A completely new interaction between two new actors X and Y (X, Y != A, B and X != Y).
arma::uvec stat_psABXY(Data_DEM &object, arma::mat &data, unsigned int &from, unsigned int &to, unsigned int &number_event, unsigned int col_number, std::string transformation, unsigned int type) {
  if (type >= 10 || type == 3) return arma::uvec();

  std::vector<unsigned int> changed_indices;

  const unsigned int n_rows = object.current_stats.data.n_rows;
  for (unsigned int idx = 0; idx < n_rows; ++idx) {
    unsigned int sender = (unsigned int)object.current_stats.data.at(idx, 3);
    unsigned int receiver = (unsigned int)object.current_stats.data.at(idx, 4);

    double old_val = object.current_stats.data.at(idx, col_number);
    double new_val = 0.0;
    if (from != 0 && to != 0) {
      if (sender != from && sender != to && receiver != from && receiver != to) {
        new_val = 1.0;
      }
    }

    if (old_val != new_val) {
      object.current_stats.set_stat(idx, col_number, new_val);
      changed_indices.push_back(idx);
    }
  }

  return arma::uvec(changed_indices);
}
TERM_REGISTER("psABXY", stat_psABXY);
