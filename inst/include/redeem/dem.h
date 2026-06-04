#pragma once
#include <RcppArmadillo.h>
#include <iostream>
#include <random>
#include <math.h>
#include "data_dem.h"
#include "sufficient_statistics.h"
#include "helper_functions.h"
#include <queue>

class DEM
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
  // 8. Simultaneous interactions or not?
  DEM(arma::mat, std::vector<std::string>, unsigned int,bool,bool,Rcpp::List,std::vector<arma::mat>, std::vector<std::string>, bool, Rcpp::List);

  // Observed quantities
  bool verbose;
  bool directed;
  bool clustered_times;
  bool simultaneous_interactions;

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
  // arma::mat preprocess_nonclustered(bool, unsigned int, unsigned int);
  // arma::mat preprocess_alt(bool, unsigned int, unsigned int);
  // arma::mat preprocess_batch(bool, unsigned int);
  // u int of the number of events, double until specific time,
  // unsigned int max number of events, arma::vec coef_0_1,
  // arma::vec coef_1_0 and then same for degree
  arma::mat sample(unsigned int n_events, double time, arma::vec coef_0_1, arma::vec coef_1_0, arma::vec coef_0_1_degree, arma::vec coef_1_0_degree, unsigned int max_events = 400000);
  // Similar to above just that you provide the time varying baseline as two vectors,
  // the first vector says something about the change points (assuming you start at 0)
  // and the second vector says the value of the baseline
  // (actually two since you have one for the 0-1 and one for 1-0)
  Rcpp::List sample_time_varying_baseline(arma::vec time_changepoints, 
                                           arma::vec baseline_0_1, 
                                           arma::vec baseline_1_0, 
                                           unsigned int n_events, 
                                           arma::vec coef_0_1, 
                                           arma::vec coef_1_0, 
                                           arma::vec coef_0_1_degree, 
                                           arma::vec coef_1_0_degree);
  // Function that takes non-degree, degree, baseline coefficients for 0_1 and 1_0 model.
  // baseline_0_1 and baseline_1_0 are per-unique-test-time vectors (one value per unique timestamp
  // in the test edgelist, in the order they appear in arma::unique(test_events.col(0))).
  // If length == 1, the single value is used for all test times (backward-compatible).
  Rcpp::List get_probabilities_per_test_event(arma::vec,arma::vec,arma::vec, arma::vec, arma::mat, bool, int,
                                               arma::vec baseline_0_1, arma::vec baseline_1_0);
  arma::vec get_oos_likelihood(arma::vec, arma::vec, arma::vec, arma::vec,
                               arma::mat, bool, arma::vec baseline_0_1, arma::vec baseline_1_0);
protected:
private:
};




