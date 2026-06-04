#pragma once
#include <RcppArmadillo.h>
#include <iostream>
#include <map>
#include "current_stat.h"
#include "hist_events.h"

class Data_DEM
{
public:
  // Needed for construction is
  // 1. Unsigned number of actors
  // 2. Number of statistics
  // 3. Directed or not?
  // 4. String of names of statistics
  // 5. Matrix of exogenous statistics for each term
  // 6. String for each term indicating the employed transformation
  // 7. Vector of integers indicating the window size for each term
  Data_DEM(unsigned int,unsigned int,bool, std::vector<std::string>, std::vector<arma::mat>, std::vector<std::string>, Rcpp::List);
  std::map<int, double> window_lengths;
  int n_nodes;
  bool directed;
  Current_Stat current_stats;
  // Current_Stat previous_stats;
  // Info given to the sufficient statistics
  // Which actors are currently connected?
  Hist_Events current_interactions;
  // Which actors have been connected sometimes in the past?
  Hist_Events general_interactions;
  std::map<int, Hist_Events> windowed_history;
  std::vector<double> changepoints;
  std::vector<std::string> term_names;
  int n_statistics;
  bool has_time_dependent_stats;
  bool is_dem;
  unsigned int last_sender = 0;
  unsigned int last_receiver = 0;
  void set_changepoints(std::vector<double>);
  void add_changepoints(double);
  void initialize_changepoints();
  // To set everything back to the original value (all current and general statistics to 0)
  void reinitialize();
protected:
private:
};



