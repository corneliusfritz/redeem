#include "data_dem.h"
// #include "sufficient_statistics.h"
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>


Data_DEM::Data_DEM(unsigned int n_nodes_, unsigned int n_statistics_, bool directed_,
                   std::vector<std::string> terms_, std::vector<arma::mat> data_list_,
                   std::vector<std::string> transformations_, Rcpp::List window_info):
  current_stats(n_nodes_, n_statistics_, directed_, data_list_),
  current_interactions(n_nodes_, directed_),
  general_interactions(n_nodes_, directed_) {
  
  n_nodes = n_nodes_;
  directed = directed_;
  n_statistics = n_statistics_;
  term_names = terms_;
  is_dem = false;
  last_sender = 0;
  last_receiver = 0;

  if (window_info.size() > 0) {
    Rcpp::CharacterVector names = window_info.names();
    for (int i = 0; i < window_info.size(); ++i) {
      int w_id = std::stoi(std::string(names[i]));
      double w_len = Rcpp::as<double>(window_info[i]);
      window_lengths[w_id] = w_len;
      windowed_history.emplace(w_id, Hist_Events(n_nodes_, directed_));
    }
  }

  has_time_dependent_stats = false;
  for (const auto& t : terms_) {
    if (t.find("current_interaction") != std::string::npos) {
      has_time_dependent_stats = true;
      break;
    }
  }
}

void Data_DEM::reinitialize(){
  current_interactions.clear();
  general_interactions.clear();
  for (auto& pair : windowed_history) {
    pair.second.clear();
  }
  current_stats.clear();
  last_sender = 0;
  last_receiver = 0;
}
void Data_DEM::set_changepoints(std::vector<double> changepoints_) {
  changepoints = changepoints_;
}

void Data_DEM::initialize_changepoints() {
  changepoints = std::vector<double>();
}

void Data_DEM::add_changepoints(double what) {
  changepoints.push_back(what);
}
// // [[Rcpp::export]]
// void trying(arma::mat edgelist_, unsigned int n_nodes_){
//   Data_DEM tmp(edgelist_, n_nodes_);
//
//   // tmp.nodes;
// }
