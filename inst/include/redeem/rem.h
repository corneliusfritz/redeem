#pragma once
#include <RcppArmadillo.h>
#include <iostream>
#include <random>
#include <math.h>
#include "data_dem.h"
#include "sufficient_statistics.h"
#include "helper_functions.h"
#include <queue>

class REM
{
public:
  // Needed for construction is
  // 1. A matrix with the structure (time = When?, from= Sender of Event, to = Receiver of Event, type = Beginning (1) or Ending (0) of Interaction? )
  // 2. String of names of statistics
  // 3. Unsigned number of actors
  // 4. Directed or not?
  // 5. Verbose or not?
  // 6. Matrix of exogenous statistics for each term
  // 7. String for each term indicating the employed transformation
  REM(arma::mat, std::vector<std::string>, unsigned int,bool,bool,Rcpp::List,std::vector<arma::mat>, std::vector<std::string>, Rcpp::List);

  //   Another constructor where some block informaiton is provided and only within edges are simulated
  REM(arma::mat, std::vector<std::string>, unsigned int,bool,bool,Rcpp::List,std::vector<arma::mat>, std::vector<std::string>, arma::vec, Rcpp::List);

  // Observed quantities
  bool verbose;
  bool directed;
  bool clustered_times;
  bool between_blocks_only;

  // This is a vector of functions
  std::vector<ValidateFunction> terms;
  arma::mat edgelist;
  unsigned int n_nodes;
  unsigned int n_entries;
  unsigned int n_changepoints;
  // When do changes occur (this might become useful when we have clustered events)
  std::vector<double> changepoints;
  // All important data needed for the preprocessing (previous and general events, current statistics)
  Data_DEM data_dem;
  std::vector<arma::mat> data_list;
  Rcpp::List original_data_list;
  std::vector<std::string> transformations;
  std::set<double> covariate_changepoints;
  std::vector<std::vector<double>> covariate_times_list;
  std::vector<int> last_covariate_indices;
  arma::uvec unavailable_nodes_indices;
  // The processed long format, which we need for estimating
  arma::mat preprocessed;
  
  // Update exogenous covariates to the values at the specific time point
  bool update_covariates(double time);

  // Sampling quantities
  // Where the events are going to be saved
  arma::mat edgelist_sample;
  double time_old;
  std::priority_queue<ScheduledEvent, std::vector<ScheduledEvent>, std::greater<ScheduledEvent>> scheduled_events;


  // Functions
  arma::mat preprocess(bool, unsigned int, unsigned int,  unsigned int, double);
  void preprocess_build(unsigned int, unsigned int);
  arma::mat preprocess_clustered(bool, unsigned int, unsigned int, double);
  arma::mat sample(unsigned int n_events, double time, arma::vec coef,
                   arma::vec coef_degree, unsigned int max_events = 400000);
  Rcpp::List sample_time_varying_baseline(arma::vec time_changepoints, 
                                           arma::vec baseline, 
                                           arma::vec coef, 
                                           arma::vec coef_degree);
  // Function that takes non-degree, degree, baseline coefficients.
  // baseline_0_1 is a per-unique-test-time vector (length 1 = scalar, backward-compatible).
  Rcpp::List get_probabilities_per_test_event(arma::vec, arma::vec,
                                               arma::mat, bool, int,
                                               arma::vec baseline_0_1);
  arma::vec get_oos_likelihood(arma::vec, arma::vec,
                               arma::mat, arma::vec baseline_0_1);
protected:
private:
};
