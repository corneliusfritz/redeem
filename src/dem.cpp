#include "helper_functions.h"
// [[Rcpp::depends(RcppProgress)]]
#include "dem.h"
#include "helper_functions.h"
#include <RcppArmadillo.h>
#include <random>
#include <set>
#include <unordered_map>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <queue>
#include <utility>
#include <limits>

DEM::DEM(arma::mat edgelist_, std::vector<std::string> terms_, unsigned int n_nodes_,bool directed_,bool verbose_,
         Rcpp::List data_list_raw, std::vector<arma::mat> data_list_, std::vector<std::string> transformations_,
         bool simultaneous_interactions_, Rcpp::List window_info):
  data_dem(n_nodes_,(unsigned int)terms_.size(), directed_, terms_, data_list_, transformations_, window_info){
  data_dem.is_dem = true;
  data_list = data_list_;
  original_data_list = data_list_raw;
  time_old = 0.0;

  transformations = transformations_;
  terms = change_statistics_generate(terms_);
  n_nodes = n_nodes_;
  edgelist = edgelist_;
  n_entries = edgelist_.n_rows;

  // 1. Collect all potential change points
  std::set<double> all_cp;

  // From edgelist
  for (unsigned int i = 0; i < edgelist_.n_rows; ++i) all_cp.insert(edgelist_(i, 0));

  // From covariate measurements
  for (R_xlen_t i = 0; i < original_data_list.size(); ++i) {
    std::vector<double> times_vec;
    if (Rcpp::is<Rcpp::List>(original_data_list[i])) {
      Rcpp::List tv_list = Rcpp::as<Rcpp::List>(original_data_list[i]);
      Rcpp::CharacterVector times = tv_list.names();
      for (R_xlen_t k = 0; k < times.size(); ++k) {
        double cp = safe_stod(std::string(times[k]), "parsing covariate names in DEM constructor");
        all_cp.insert(cp);
        covariate_changepoints.insert(cp);
        times_vec.push_back(cp);
      }
    }
    covariate_times_list.push_back(times_vec);
    last_covariate_indices.push_back(0);
  }

  // 2. Update exogenous stats with initial data_list
  for (unsigned int i = 0; i < (unsigned int)data_list.size(); ++i) {
    data_dem.current_stats.update_exogenous(i, data_list[i]);
  }

  changepoints = std::vector<double>(all_cp.begin(), all_cp.end());
  n_changepoints = changepoints.size();

  // Set clustered_times to true if we have more changepoints than edgelist rows
  // Or if multiple edgelist rows have the same time
  std::set<double> unique_edge_times;
  for (unsigned int i = 0; i < edgelist_.n_rows; ++i) unique_edge_times.insert(edgelist_(i, 0));

  data_dem.set_changepoints(changepoints);
  clustered_times = (n_changepoints != edgelist_.n_rows) || (unique_edge_times.size() != edgelist_.n_rows);
  directed = directed_;
  verbose = verbose_;
  simultaneous_interactions = simultaneous_interactions_;
}

// Update exogenous covariates to the values at the specific time point
bool DEM::update_covariates(double time) {
  bool changed = false;
  // Search for the latest time point for each covariate that is less than or equal to the current time
  // (to find the correct time step of the exogenous covariates given the current time).
  for (unsigned int i = 0; i < (unsigned int)covariate_times_list.size(); ++i) {
    // Skip if the covariate has no time points.
    if (covariate_times_list[i].empty()) continue;

    int current_idx = last_covariate_indices[i];
    const std::vector<double>& times = covariate_times_list[i];

    // If the current time is larger than the last recorded time for this covariate,
    // we might need to update the covariate values to the values at the current time point.
    // If the new time point is larger than the current time point, we need to update the covariate values.
    // We continue to update the covariate values until the current time point is larger than the last recorded time point for this covariate.
    // This is necessary because the covariates are updated at each changepoint, and the changepoints are not necessarily equally spaced.
    bool advanced = false;
    while (current_idx + 1 < (int)times.size() && times[current_idx + 1] <= time + 1e-10) {
      current_idx++;
      advanced = true;
    }
    // If the covariate values have advanced, update the covariate values.
    if (advanced) {
      last_covariate_indices[i] = current_idx;
      Rcpp::List tv_list = Rcpp::as<Rcpp::List>(original_data_list[i]);
      arma::mat raw_mat = Rcpp::as<arma::mat>(tv_list[current_idx]);
      double K = (raw_mat.n_elem > 0) ? raw_mat.at(0, 0) : 1.0;
      arma::mat new_mat = apply_transformation(raw_mat, transformations[i], K);
      if (data_dem.current_stats.update_exogenous(i, new_mat)) {
        data_list[i] = new_mat;
        changed = true;
      }
    }
  }
  return changed;
}


Rcpp::List DEM::sample_time_varying_baseline(arma::vec time_changepoints,
                                             arma::vec baseline_0_1,
                                             arma::vec baseline_1_0,
                                             unsigned int n_events,
                                             arma::vec coef_0_1,
                                             arma::vec coef_1_0,
                                             arma::vec coef_0_1_degree ,
                                             arma::vec coef_1_0_degree){
  Rcpp::List res (time_changepoints.size());
  // It is assumed that the simulation starts at 0 and that all parameters are scaled to be a difference to the first time point
  arma::vec new_vec = arma::vec(1).fill(0.0);
  time_changepoints = arma::join_vert(new_vec, time_changepoints);
  baseline_0_1 = arma::join_vert(new_vec, baseline_0_1);
  baseline_1_0 = arma::join_vert(new_vec, baseline_1_0);
  // This is an indicator whether or not the sample includes degree parameters
  // bool non_degree = coef_0_1_degree.size() == 1;
  arma::vec eff_coef_0_1 = coef_0_1;
  arma::vec eff_coef_1_0= coef_1_0;
  for(unsigned int i = 1; i < time_changepoints.size(); i++){
    // If no degree terms are included, we assume that the first coefficient relates to the intercept
    arma::vec current_coef_0_1 = coef_0_1;
    arma::vec current_coef_1_0 = coef_1_0;
    current_coef_0_1.at(0) += baseline_0_1.at(i-1);
    current_coef_1_0.at(0) += baseline_1_0.at(i-1);
    res.at(i-1) = sample(n_events, time_changepoints.at(i) - time_changepoints.at(i-1),
           current_coef_0_1, current_coef_1_0, coef_0_1_degree, coef_1_0_degree);
  }
  return(res);
}


arma::mat DEM::sample(unsigned int n_events, double time, arma::vec coef_0_1, arma::vec coef_1_0, arma::vec coef_0_1_degree, arma::vec coef_1_0_degree, unsigned int max_events) {
  // Determine buffer size
  unsigned int buffer_size = (n_events > 0) ? n_events : 40000;
  arma::mat event_mat(buffer_size, 4);

  arma::uvec currently_interacting(data_dem.current_stats.data.n_rows), currently_noninteracting(data_dem.current_stats.data.n_rows);

  arma::uvec tmp_changes;
  unsigned int from_, to_, type_, tmp_pair_; double time_;
  double time_limit = time_old + time;
  Progress p(n_events, verbose);
  data_dem.initialize_changepoints();

  arma::vec coef_degree_1_0_from, coef_degree_1_0_to;
  arma::vec coef_degree_0_1_from, coef_degree_0_1_to;

  if(coef_0_1_degree.size() > 1){
    coef_degree_0_1_from = coef_0_1_degree.rows(arma::conv_to<arma::uvec>::from(data_dem.current_stats.data.col(3)-1));
    coef_degree_0_1_to = coef_0_1_degree.rows(arma::conv_to<arma::uvec>::from(data_dem.current_stats.data.col(4)-1 + (directed ? n_nodes : 0)));
  } else {
    coef_degree_0_1_from = arma::zeros(data_dem.current_stats.data.n_rows);
    coef_degree_0_1_to = arma::zeros(data_dem.current_stats.data.n_rows);
  }
  if(coef_1_0_degree.size() > 1){
    coef_degree_1_0_from = coef_1_0_degree.rows(arma::conv_to<arma::uvec>::from(data_dem.current_stats.data.col(3)-1));
    coef_degree_1_0_to = coef_1_0_degree.rows(arma::conv_to<arma::uvec>::from(data_dem.current_stats.data.col(4)-1 + directed*n_nodes));
  } else {
    coef_degree_1_0_from = arma::zeros(data_dem.current_stats.data.n_rows);
    coef_degree_1_0_to = arma::zeros(data_dem.current_stats.data.n_rows);
  }
  unsigned int m = 0;
  // Initialize fast RNG seeded from R
  GetRNGstate();
  unsigned int master_seed = (unsigned int)R::runif(0.0, 4294967295.0);
  PutRNGstate();
  std::mt19937 local_gen(master_seed);
  std::uniform_real_distribution<double> local_dist(0.0, 1.0);

  arma::vec intensity(data_dem.current_stats.data.n_rows);

  // Initialize scalable intensities once
  data_dem.current_stats.initialize_intensities(coef_0_1, coef_0_1_degree, coef_1_0, coef_1_0_degree, simultaneous_interactions);

  unsigned int target_m = (n_events > 0) ? n_events : max_events;
  while (m < target_m && (n_events > 0 || (n_events == 0 && time_old < time_limit))) {
    Rcpp::checkUserInterrupt();
    if (verbose && n_events > 0) p.increment();

    // 1. Get intensities from stateful cache (O(1) access)
    intensity = data_dem.current_stats.combined_intensity;
    // If we observe any values that are not finite, we have a problem and should inspect it
    if(intensity.is_finite() == 0){
      if (verbose) Rcpp::Rcout << "Caution! Some sufficient statistics were infinite ... " << std::endl;
    }

    double t_next_scheduled = scheduled_events.empty() ? R_PosInf : (scheduled_events.top().time - time_old);

    // 3. Sample a Exponential-distributed RV where lambda is sum(intensities) -> This will be the interevent time
    double total_intensity = arma::sum(intensity);
    if (total_intensity < 1e-15 && scheduled_events.empty()) {
       if (time > 0) time_old = time_limit;
       break;
    }


    if (total_intensity < 1e-15) {
        time_ = R_PosInf;
    } else {
        double u_time = local_dist(local_gen);
        if (u_time < 1e-15) u_time = 1e-15;
        time_ = -log(u_time) / total_intensity;
    }

    if (time_ < t_next_scheduled) {
        if(time > 0 && time_old + time_ > time_limit){
            time_old = time_limit;
            break;
        }

        data_dem.add_changepoints(time_ + time_old);
        time_old += time_;

        // Advance time-dependent statistics (durations) for all interacting pairs
        if (data_dem.has_time_dependent_stats) {
          advance_delta_time(data_dem, m, data_list, transformations, terms, tmp_changes);
          if (!tmp_changes.is_empty()) {
            data_dem.current_stats.update_intensities_at_indices(tmp_changes, coef_0_1, coef_1_0);
          }
        }

        // 4. Sample from multinomial distribution to determine which event actually occurs
        arma::vec tmp = intensity/arma::sum(intensity);
        double u = local_dist(local_gen);
        arma::vec cumulative = arma::cumsum(tmp);
        tmp_pair_ = std::lower_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin();
        if (tmp_pair_ >= data_dem.current_stats.data.n_rows) {
            tmp_pair_ = data_dem.current_stats.data.n_rows - 1;
        }

        from_ = data_dem.current_stats.data(tmp_pair_, 3);
        to_ = data_dem.current_stats.data(tmp_pair_, 4);
        type_ = 1 - data_dem.current_stats.get_status(from_, to_);
        // Add sampled event to the event matrix
        if (m >= buffer_size) {
            event_mat.insert_rows(m, 1000); // Dynamic expansion
            buffer_size += 1000;
        }
        event_mat.row(m) = {time_old, static_cast<double>(from_), static_cast<double>(to_), static_cast<double>(type_)};
        m++;


        // 5. Update the statistics accordingly (just repeat what was implemented in the preprocess step)
        if (type_ != 3) {
           compute_changes(from_, to_, m, type_, data_dem, data_list, transformations, terms, simultaneous_interactions, tmp_changes);
        }

        if(type_ == 0){
          // If the observed event was an end to an interaction it will be excluded from the current_interactions
          data_dem.current_interactions.exclude_event(from_, to_);
          // Now an interaction between from_ and to_ ends -> set corresponding status to 0
          data_dem.current_stats.set_status(from_, to_, 0);
          // Here implement that the avail values are updated according to the end event -> set avail or from_ and to_ to 1
          if(!simultaneous_interactions){
            data_dem.current_stats.now_avail(from_);
            data_dem.current_stats.now_avail(to_);
          }
        } else if (type_ == 1) {
          // If the observed event was a beginning to an interaction it will be added from the current_interactions and general_interactions
          data_dem.current_interactions.add_event(from_, to_);
          data_dem.general_interactions.add_event(from_, to_);

          // Add to windowed histories and schedule dissolution
          for (auto& pair : data_dem.windowed_history) {
              pair.second.add_event(from_, to_);
              unsigned int w_id = pair.first;
              // Use +10 offset for windowed expirations to distinguish them from formation (1) and dissolution (0)
              scheduled_events.push({time_old + data_dem.window_lengths[w_id], (unsigned int)from_, (unsigned int)to_, w_id + 10});
          }
          // Now an interaction between from_ and to_ begins -> set corresponding status to 1
          data_dem.current_stats.set_status(from_, to_, 1);
          // Here implement that the avail values are updated according to the start event -> set avail or from_ and to_ to 0
          if(!simultaneous_interactions){
            data_dem.current_stats.not_avail(from_);
            data_dem.current_stats.not_avail(to_);
          }
        }

        // Only update the dyads that actually changed, AFTER status and interaction updates
        if (type_ != 3) {
          data_dem.current_stats.update_intensities_at_indices(tmp_changes, coef_0_1, coef_1_0);
        }
    } else {
        // Scheduled event occurs
        if(time > 0 && time_old + t_next_scheduled > time_limit){
            time_old = time_limit;
            break;
        }

        data_dem.add_changepoints(t_next_scheduled + time_old);
        time_old += t_next_scheduled;

        // Advance time-dependent statistics (durations) for all interacting pairs
        if (data_dem.has_time_dependent_stats) {
          advance_delta_time(data_dem, m, data_list, transformations, terms, tmp_changes);
          if (!tmp_changes.is_empty()) {
            data_dem.current_stats.update_intensities_at_indices(tmp_changes, coef_0_1, coef_1_0);
          }
        }

        ScheduledEvent se = scheduled_events.top();
        scheduled_events.pop();

         compute_changes(se.from, se.to, m, se.type, data_dem, data_list, transformations, terms, simultaneous_interactions, tmp_changes);

        // Update global state before intensities
        if (se.type >= 10) {
            // Subtract 10 to get the original window ID for history lookup
            auto it = data_dem.windowed_history.find(se.type - 10);
            if (it != data_dem.windowed_history.end()) {
                it->second.exclude_event(se.from, se.to);
            }
        }

        data_dem.current_stats.update_intensities_at_indices(tmp_changes, coef_0_1, coef_1_0);

    }
  }


  if(m == 0){
    return(arma::mat(0, 4));
  } else {
    return(event_mat.rows(0, m-1));
  }
}
// arma::mat dem::preprocess_batch(bool reset, unsigned int batch_size){
//   int start = 0;
//   while(start < n_entries -batch_size) {
//     preprocess(false,start,start + batch_size-1 );
//     start += batch_size ;
//   }
//   if(start!= n_entries){
//     preprocess(reset,start,n_entries);
//   }
//   return(preprocessed);
// }

arma::mat DEM::preprocess(bool reset, unsigned int start, unsigned int end, unsigned int train_start, double max_time){
  arma::mat res;
  if(reset) data_dem.reinitialize();
  if(train_start == 0) {
    if(clustered_times){
      if(end== edgelist.n_rows) {
        res = preprocess_clustered(reset, start, n_changepoints, max_time);
      } else {
        res = preprocess_clustered(reset, start, end, max_time);
      }
    } else {
      res = preprocess_clustered(reset, start, end, max_time);
    }
  } else {
    if(clustered_times){
      if(end== edgelist.n_rows) {
        preprocess_build(start, train_start);
        res = preprocess_clustered(reset, train_start, n_changepoints, max_time);
      } else {

        preprocess_build(start, train_start);
        res = preprocess_clustered(reset, train_start, (end == edgelist.n_rows) ? n_changepoints : end, max_time);
      }

    } else {

      preprocess_build(start, train_start);
      res = preprocess_clustered(reset, train_start, (end == edgelist.n_rows) ? n_changepoints : end, max_time);
    }
  }
  return(res);
}

void DEM:: preprocess_build(unsigned int start, unsigned int end){
  arma::uvec tmp_changes;
  Progress p(end - start, verbose);
  double time_;
  // If the preprocessing starts at the very beginning, we set the previous time to 0 (in other words, we assume that the data collection starts at 0)
  if(start == 0){
    time_ = 0;
  } else {
    time_ = data_dem.changepoints.at(start-1);
  }
  // data_dem.changepoints>;
  // We go through all events, calculate the changes of it on current_stats, save all the pairs where covariates changed
  for(unsigned int i = start; i < end; i++) {
    Rcpp::checkUserInterrupt();
    p.increment(); // update progress
    time_ = data_dem.changepoints.at(i);

    // Pass 0: Advance time-dependent statistics (durations) for all interacting pairs
    if (data_dem.has_time_dependent_stats) {
      advance_delta_time(data_dem, i, data_list, transformations, terms);
    }

    // Pass 1: Handle covariate updates (Refresh state)
    if (covariate_changepoints.count(time_)) {
      update_covariates(time_);
    }

    arma::uvec tmp_ind = arma::find(edgelist.col(0) == time_);

    for(unsigned int j = 0; j < tmp_ind.size(); j++) {
      unsigned int from_ = (unsigned int) edgelist.at(tmp_ind.at(j),1);
      unsigned int to_ = (unsigned int) edgelist.at(tmp_ind.at(j),2);
      unsigned int type_= (unsigned int) edgelist.at(tmp_ind.at(j),3);
      // If event is of type 3, this indicates that there is an exogenous change to the baseline intensity
      // -> no covariates are changed but all are saved again (only the time changes)
      if(type_ != 3) {
        // Call the function that changes the statistics in current_stats and return the indicators of pairs where statistics change
         compute_changes(from_,to_,i,type_,data_dem, data_list, transformations, terms, simultaneous_interactions, tmp_changes);
      }
      // Update the event histories
      if(type_ >= 10){
        // Windowed dissolution
        auto it = data_dem.windowed_history.find(type_ - 10);
        if (it != data_dem.windowed_history.end()) {
          it->second.exclude_event(from_, to_);
        }
      } else if(type_ == 0){
        // If the observed event was an end to an interaction it will be excluded from the current_interactions
        data_dem.current_interactions.exclude_event(from_, to_);
        // Now an interaction between from_ and to_ ends -> set corresponding status to 0
        data_dem.current_stats.set_status(from_, to_, 0);
        // Here implement that the avail values are updated according to the end event -> set avail or from_ and to_ to 1
        if(!simultaneous_interactions){
          data_dem.current_stats.now_avail(from_);
          data_dem.current_stats.now_avail(to_);
        }
      } else if(type_ == 1){
        // If the observed event was a beginning to an interaction it will be added from the current_interactions and general_interactions
        data_dem.current_interactions.add_event(from_, to_);
        data_dem.general_interactions.add_event(from_, to_);
        // Add to all windowed histories
        for (auto& pair : data_dem.windowed_history) {
          pair.second.add_event(from_, to_);
        }
        // Now an interaction between from_ and to_ begins -> set corresponding status to 1
        data_dem.current_stats.set_status(from_, to_, 1);
        if(!simultaneous_interactions){
          data_dem.current_stats.not_avail(from_);
          data_dem.current_stats.not_avail(to_);
        }
      }
    }
    // Set the boolean variable indicating which event was actually happening back to the starting point (where nothing happened)
    for(unsigned int j = 0; j < tmp_ind.size(); j++) {
      unsigned int from_ =  (unsigned int) edgelist.at(tmp_ind.at(j),1);
      unsigned int to_ =  (unsigned int) edgelist.at(tmp_ind.at(j),2);
      data_dem.current_stats.reset_event(from_,to_);
    }
  }
}



arma::mat DEM::preprocess_clustered(bool reset, unsigned int start, unsigned int end, double max_time){
  arma::uvec tmp_changes;
  arma::uvec all_changes;
  // cols is just a uvec enumerating all columns of the current_stats matrix (needed to copy them in the iterative scheme)
  arma::uvec cols = arma::regspace<arma::uvec>(0,data_dem.current_stats.data.n_cols-1);
  unsigned int n_entries_tmp = end - start;
  // We will save the results in a list object with one additional entry (which is the state of everything at the very beginning)
  List res(n_entries_tmp+1);
  double time_;
  // If the preprocessing starts at the very beginning, we set the previous time to 0 (in other words, we assume that the data collection starts at 0)
  if(start == 0){
    time_ = 0;
  } else {
    time_ = data_dem.changepoints.at(start-1);
  }
  // At the beginning, we need the info on the covariates at the very beginning
  arma::mat initial_data = data_dem.current_stats.data;

  if(data_dem.directed){
    arma::vec times(n_nodes*(n_nodes-1));
    times.fill(time_);
    res.at(0) = arma::join_rows(times, initial_data);
  } else {
    arma::vec times(n_nodes*(n_nodes-1)/2);
    times.fill(time_);
    res.at(0) = arma::join_rows(times, initial_data);
  }
  Progress p(end - start, verbose);

  // Optimization: Check which terms are actually time-dependent (current_interaction)
  bool any_time_dependent = false;
  for (const auto& name : data_dem.term_names) {
      if (name.find("current_interaction") != std::string::npos) {
          any_time_dependent = true;
          break;
      }
  }

  for(unsigned int i = start; i < end; i++) {
    Rcpp::checkUserInterrupt();
    p.increment(); // update progress
    // Here we now do not iterate over the events but the event times
    time_ =  data_dem.changepoints.at(i);

    // Pass 0: Advance time-dependent statistics (durations) for all interacting pairs
    if (any_time_dependent) {
      advance_delta_time(data_dem, i, data_list, transformations, terms, tmp_changes);
      if (!tmp_changes.is_empty()) {
        all_changes.insert_rows(all_changes.n_rows, tmp_changes);
      }
    }

    arma::uvec tmp_ind = arma::find(edgelist.col(0) == time_);

    // Pass 1: Handle covariate updates (Refresh state)
    if (covariate_changepoints.count(time_)) {
      update_covariates(time_);
      tmp_changes = arma::regspace<arma::uvec>(0, 1, data_dem.current_stats.data.n_rows - 1);
      if (tmp_changes.n_elem > 0) {
        all_changes.insert_rows(all_changes.n_rows, tmp_changes);
      }
    }

    // Pass 2: Process interaction events using updated state
    for(unsigned int j = 0; j < (unsigned int)tmp_ind.size(); j++) {
      unsigned int from_ = (unsigned int) edgelist.at(tmp_ind.at(j), 1);
      unsigned int to_ = (unsigned int) edgelist.at(tmp_ind.at(j), 2);
        unsigned int type_= (unsigned int) edgelist.at(tmp_ind.at(j), 3);

        // Record the focal dyad's current state regardless of whether any statistics changed.
        // This ensures the preprocessed output always includes a row for every interaction event.
        unsigned int pair_id = data_dem.current_stats.find_from_to(from_, to_);
        all_changes.insert_rows(all_changes.n_rows, arma::uvec{pair_id});

        if (type_ == 3) {
          // Type 3 represents an exogenous change to the baseline intensity or a time-varying coefficient changepoint.
          // These are hard boundaries for all dyads and should not be compressed.
          tmp_changes = arma::regspace<arma::uvec>(0, 1, data_dem.current_stats.data.n_rows - 1);
          if (tmp_changes.n_elem > 0) {
            // Mark all dyads as changed for this timestamp
            all_changes.insert_rows(all_changes.n_rows, tmp_changes);
            // We also record it in the data matrix
            data_dem.current_stats.data.col(2).fill(0.0);
          }
        } else {
          // Call the function that changes the statistics in current_stats and return the indicators of pairs where statistics change
           compute_changes(from_, to_, i, type_, data_dem, data_list, transformations, terms, simultaneous_interactions, tmp_changes);
          if (tmp_changes.n_elem > 0) {
            all_changes.insert_rows(all_changes.n_rows, tmp_changes);
          }

          // Mark the focal dyad as having a physical event (Start or End)
          if (type_ == 1 || type_ == 0) {
              data_dem.current_stats.data.at(pair_id, 2) = 1.0;
          } else {
              data_dem.current_stats.data.at(pair_id, 2) = 0.0;
          }

        // Update the event histories
        if (type_ >= 10) {
          // Windowed dissolution only affects the specific windowed history term
          auto it = data_dem.windowed_history.find(type_ - 10);
          if (it != data_dem.windowed_history.end()) {
            it->second.exclude_event(from_, to_);
          }
        } else if (type_ == 0) {
          // Observed interaction ends -> exclude from current_interactions
          data_dem.current_interactions.exclude_event(from_, to_);
          // Set corresponding status to 0
          data_dem.current_stats.set_status(from_, to_, 0);

          // Update availability according to the end event
          if (!simultaneous_interactions) {
            all_changes.insert_rows(all_changes.n_rows, data_dem.current_stats.now_avail(from_));
            all_changes.insert_rows(all_changes.n_rows, data_dem.current_stats.now_avail(to_));
          }
        } else if (type_ == 1) {
          // If the observed event was a beginning to an interaction it will be added from the current_interactions and general_interactions
          data_dem.current_interactions.add_event(from_, to_);
          data_dem.general_interactions.add_event(from_, to_);
          for (auto& pair : data_dem.windowed_history) {
            pair.second.add_event(from_, to_);
          }
          // Now an interaction between from_ and to_ begins -> set corresponding status to 1
          data_dem.current_stats.set_status(from_, to_, 1);

          // Here implement that the avail values are updated according to the begin event -> set avail or from_ and to_ to 0
          if (!simultaneous_interactions) {
            all_changes.insert_rows(all_changes.n_rows, data_dem.current_stats.not_avail(from_));
            all_changes.insert_rows(all_changes.n_rows, data_dem.current_stats.not_avail(to_));
          }
        }
      }
    }

    all_changes = arma::unique(all_changes);
    arma::vec times(all_changes.size(), arma::fill::value(time_));
    res.at(i - start + 1) = arma::join_rows(times, data_dem.current_stats.data.rows(all_changes));

    for(unsigned int j = 0; j < (unsigned int)tmp_ind.size(); j++) {
      unsigned int from_ = (unsigned int) edgelist.at(tmp_ind.at(j), 1);
      unsigned int to_ = (unsigned int) edgelist.at(tmp_ind.at(j), 2);
      data_dem.current_stats.reset_event(from_, to_);
    }

    // Reset events for next cluster
    data_dem.current_stats.data.col(2).zeros();

    all_changes.clear();
  }
  // Only return the rows that were actually used
  double m_time = std::max(time_, max_time);
  arma::mat res_final = do_call_full_cpp(res, m_time);

  //   // Just a small trick to order according to pairs and times simultaneous
  //   arma::uvec indices = sort_index(res_final.col(1) +res_final.col(0)/res_final.n_rows);
  //   res_final = res_final.rows(indices);
  //
  //   // Cut the first entry and replace all zeroes in the resulting matrix with the maximal time (which is time_)
  //   arma::vec time_new =  res_final.submat(1,0,res_final.n_rows-1,0);
  //   time_new.resize(time_new.size()+1);
  //   time_new.replace(0.0,time_);
  //   // Do the same with the event indicator
  //   arma::vec event_new =  res_final.submat(1,3,res_final.n_rows-1,3);
  //   event_new.resize(event_new.size()+1);
  //   res_final.col(3) = event_new;
  //
  //
  //   res_final = arma::join_horiz(time_new,res_final);
  // Find where res_final.col(0) != res_final.col(1) (only those values we need to return)
  // arma::uvec tmp_indicator = arma::find(res_final.col(0) != res_final.col(1));
  // If there is already a preprocessed file available then join them (needed for preprocessing in batches)
  if(preprocessed.n_elem== 0){

    preprocessed = res_final;
  } else {
    preprocessed = arma::join_vert(preprocessed,res_final);
    arma::uvec indices = arma::regspace<arma::uvec>(0, preprocessed.n_rows - 1);
    std::sort(indices.begin(), indices.end(), [&](arma::uword a, arma::uword b) {
      if (preprocessed(a, 1) != preprocessed(b, 1)) return preprocessed(a, 1) < preprocessed(b, 1);
      return preprocessed(a, 0) < preprocessed(b, 0);
    });
    preprocessed = preprocessed.rows(indices);
  }
  // We just reinitialize the internal data if needed
  // (this is important, if we first estimate and then want to sample from the beginning
  //  but is not needed if we then want to simulate the next events)
  if(reset){
    data_dem.reinitialize();
  }


  return(preprocessed);
}

Rcpp::List DEM::get_probabilities_per_test_event(arma::vec coef_0_1, arma::vec coef_0_1_degree,
                                                 arma::vec coef_1_0 , arma::vec coef_1_0_degree,
                                                 arma::mat test_events, bool simultaneous_interactions, int k,
                                                 arma::vec baseline_0_1, arma::vec baseline_1_0){


  // This is an indicator whether the sample should go on until a specific time is met or until a specific number of events are observed
  arma::vec lp(data_dem.current_stats.data.n_rows);
  // 7 and +6 should always point towards the first and last-to-first index that relates to the statistics
  arma::uvec cols_statistics;
  if (terms.size() > 0) {
    cols_statistics = arma::regspace<arma::uvec>(7, terms.size() + 6);
  }
  arma::uvec tmp_changes;

  double time_;

  arma::mat predicted, observed;
  // Go over each unique time in the test_events
  arma::vec unique_times = arma::unique(test_events.col(0));

  unsigned int to_, from_, type_, pair_;
  Rcpp::List res(unique_times.size());
  arma::uvec ord;
  Progress p(unique_times.size(), verbose);

  // Use the first baseline value for initialization
  double init_b_0_1 = (baseline_0_1.n_elem > 0) ? baseline_0_1.at(0) : 0.0;
  double init_b_1_0 = (baseline_1_0.n_elem > 0) ? baseline_1_0.at(0) : 0.0;

  // Initialize stateful intensity cache once
  data_dem.current_stats.initialize_intensities(coef_0_1, coef_0_1_degree, coef_1_0, coef_1_0_degree, simultaneous_interactions, init_b_0_1, init_b_1_0);
  // Initialize changepoints for time-dependent stats advancement
  data_dem.set_changepoints(arma::conv_to<std::vector<double>>::from(unique_times));

  for(unsigned int m = 0; m<unique_times.size(); m++){
    Rcpp::checkUserInterrupt();
    p.increment();
    time_ =  unique_times.at(m);

    // Update baseline for this test time (per-event baseline vector support)
    {
      unsigned int b_idx = std::min((unsigned int)m, (unsigned int)baseline_0_1.n_elem - 1);
      unsigned int b_idx_1 = std::min((unsigned int)m, (unsigned int)baseline_1_0.n_elem - 1);
      double b_0_1 = baseline_0_1.at(b_idx);
      double b_1_0 = baseline_1_0.at(b_idx_1);
      data_dem.current_stats.update_baseline(b_0_1, b_1_0, coef_0_1, coef_1_0);
    }

    // 0. Update time-varying covariates for this prediction interval
    bool changed = update_covariates(time_);

    // Get the index of the events that are currently predicted
    arma::uvec tmp_ind = arma::find(test_events.col(0) == time_);

    // 1. Get intensities from stateful cache (O(1) access)
    // Refresh intensities ONLY after covariate update OR if there are time-dependent stats
    if (data_dem.has_time_dependent_stats) {
        arma::uvec time_changes = advance_delta_time(data_dem, m, data_list, transformations, terms);
        if (changed) {
            data_dem.current_stats.update_intensities_at_indices(arma::regspace<arma::uvec>(0, data_dem.current_stats.data.n_rows - 1), coef_0_1, coef_1_0);
        } else if (!time_changes.is_empty()) {
            data_dem.current_stats.update_intensities_at_indices(time_changes, coef_0_1, coef_1_0);
        }
    } else if (changed) {
        data_dem.current_stats.update_intensities_at_indices(arma::regspace<arma::uvec>(0, data_dem.current_stats.data.n_rows - 1), coef_0_1, coef_1_0);
    }
    lp = data_dem.current_stats.lp_0_1;
    if (!simultaneous_interactions) {
        arma::uvec blocked = arma::find(data_dem.current_stats.data.col(5) == 0.0 || data_dem.current_stats.data.col(6) == 0.0);
        if (!blocked.is_empty()) {
            lp.elem(blocked).fill(-arma::datum::inf);
        }
    }
    arma::uvec interacting = data_dem.current_stats.get_currently_interacting();
    if (!interacting.is_empty()) {
        lp.elem(interacting) = data_dem.current_stats.lp_1_0.elem(interacting);
    }

    arma::uword K = k - 1 + tmp_ind.n_elem;
    ord = topk_heap(lp, K);

    predicted = data_dem.current_stats.data.submat(ord, arma::uvec{3u, 4u});
    predicted = arma::join_rows(predicted, lp.elem(ord));

    // Get the pair of events that interacted
    observed = test_events.submat(tmp_ind, arma::uvec{1u, 2u});
    res.at(m) = Rcpp::List::create(
      Rcpp::Named("observed")  = observed,
      Rcpp::Named("predicted") = predicted
    );
    tmp_changes.reset();
    for(unsigned int i =0; i<tmp_ind.size(); i++){
      from_ = test_events.at(tmp_ind.at(i),1);
      to_ =  test_events.at(tmp_ind.at(i),2);
      pair_ = data_dem.current_stats.find_from_to(from_, to_);
      tmp_changes.insert_rows(tmp_changes.n_rows, arma::uvec{pair_});
      type_= test_events.at(tmp_ind.at(i),3);
      // Update the statistics accordingly (just repeat what was implemented in the preprocess step)
      tmp_changes.insert_rows(tmp_changes.n_rows,
                              compute_changes(from_,to_,m,type_,data_dem, data_list, transformations, terms, simultaneous_interactions));
      if(type_ == 0){
        // If the observed event was an end to an interaction it will be excluded from the current_interactions
        data_dem.current_interactions.exclude_event(from_, to_);
        // Now an interaction between from_ and to_ ends -> set corresponding status to 0
        data_dem.current_stats.set_status(from_, to_, 0);
        // Here implement that the avail values are updated according to the end event -> set avail or from_ and to_ to 1
        if(!simultaneous_interactions){
          data_dem.current_stats.now_avail(from_);
          data_dem.current_stats.now_avail(to_);
        }
      } else if (type_ == 1) {
        // If the observed event was a beginning to an interaction it will be added from the current_interactions and general_interactions
        data_dem.current_interactions.add_event(from_, to_);
        data_dem.general_interactions.add_event(from_, to_);
        for (auto& pair : data_dem.windowed_history) {
          pair.second.add_event(from_, to_);
        }
        // Now an interaction between from_ and to_ begins -> set corresponding status to 1
        data_dem.current_stats.set_status(from_, to_, 1);
        // Here implement that the avail values are updated according to the start event -> set avail or from_ and to_ to 0
        if(!simultaneous_interactions){
          data_dem.current_stats.not_avail(from_);
          data_dem.current_stats.not_avail(to_);
        }
      } else if (type_ >= 10) {
        auto it = data_dem.windowed_history.find((unsigned int)type_ - 10);
        if (it != data_dem.windowed_history.end()) {
          it->second.exclude_event(from_, to_);
        }
      }
    }
    tmp_changes = arma::unique(tmp_changes);
    data_dem.current_stats.update_intensities_at_indices(tmp_changes, coef_0_1, coef_1_0);
  }

  return(res);



}

arma::vec DEM::get_oos_likelihood(arma::vec coef_0_1, arma::vec coef_0_1_degree,
                                  arma::vec coef_1_0, arma::vec coef_1_0_degree,
                                  arma::mat test_events, bool simultaneous_interactions,
                                  arma::vec baseline_0_1, arma::vec baseline_1_0) {
  // test_events columns: time, from, to, type
  arma::vec lp(data_dem.current_stats.data.n_rows);
  arma::vec unique_times = arma::unique(test_events.col(0));
  
  arma::vec event_log_likelihoods(test_events.n_rows);
  
  // Use the first baseline value for initialization
  double init_b_0_1 = (baseline_0_1.n_elem > 0) ? baseline_0_1.at(0) : 0.0;
  double init_b_1_0 = (baseline_1_0.n_elem > 0) ? baseline_1_0.at(0) : 0.0;

  // Initialize stateful intensity cache once
  data_dem.current_stats.initialize_intensities(coef_0_1, coef_0_1_degree, coef_1_0, coef_1_0_degree, simultaneous_interactions, init_b_0_1, init_b_1_0);
  data_dem.set_changepoints(arma::conv_to<std::vector<double>>::from(unique_times));
  
  arma::uvec tmp_changes;
  
  for(unsigned int m = 0; m < (unsigned int)unique_times.size(); m++){
    Rcpp::checkUserInterrupt();
    double time_ = unique_times.at(m);

    // Update baseline for this test time
    {
      unsigned int b_idx   = std::min((unsigned int)m, (unsigned int)baseline_0_1.n_elem - 1);
      unsigned int b_idx_1 = std::min((unsigned int)m, (unsigned int)baseline_1_0.n_elem - 1);
      data_dem.current_stats.update_baseline(baseline_0_1.at(b_idx), baseline_1_0.at(b_idx_1),
                                             coef_0_1, coef_1_0);
    }

    // Update time-varying covariates
    bool changed = update_covariates(time_);
    
    // Update intensities
    if (data_dem.has_time_dependent_stats) {
        arma::uvec time_changes = advance_delta_time(data_dem, m, data_list, transformations, terms);
        if (changed) {
            data_dem.current_stats.update_intensities_at_indices(arma::regspace<arma::uvec>(0, data_dem.current_stats.data.n_rows - 1), coef_0_1, coef_1_0);
        } else if (!time_changes.is_empty()) {
            data_dem.current_stats.update_intensities_at_indices(time_changes, coef_0_1, coef_1_0);
        }
    } else if (changed) {
        data_dem.current_stats.update_intensities_at_indices(arma::regspace<arma::uvec>(0, data_dem.current_stats.data.n_rows - 1), coef_0_1, coef_1_0);
    }
    
    lp = data_dem.current_stats.lp_0_1;
    if (!simultaneous_interactions) {
        arma::uvec blocked = arma::find(data_dem.current_stats.data.col(5) == 0.0 || data_dem.current_stats.data.col(6) == 0.0);
        if (!blocked.is_empty()) {
            lp.elem(blocked).fill(-arma::datum::inf);
        }
    }
    arma::uvec interacting = data_dem.current_stats.get_currently_interacting();
    if (!interacting.is_empty()) {
        lp.elem(interacting) = data_dem.current_stats.lp_1_0.elem(interacting);
    }
    
    // Compute log-sum-exp of lp
    double lse = log_sum_exp(lp);
    
    arma::uvec tmp_ind = arma::find(test_events.col(0) == time_);
    
    for(unsigned int i = 0; i < (unsigned int)tmp_ind.size(); i++){
      unsigned int idx = tmp_ind.at(i);
      unsigned int from_ = test_events.at(idx, 1);
      unsigned int to_ = test_events.at(idx, 2);
      unsigned int pair_ = data_dem.current_stats.find_from_to(from_, to_);
      unsigned int type_ = test_events.at(idx, 3);
      
      // Store log-likelihood of this observed event
      double log_prob = lp[pair_] - lse;
      event_log_likelihoods[idx] = log_prob;
      
      // Update statistics
      tmp_changes.insert_rows(tmp_changes.n_rows, arma::uvec{pair_});
      tmp_changes.insert_rows(tmp_changes.n_rows,
                              compute_changes(from_, to_, m, type_, data_dem, data_list, transformations, terms, simultaneous_interactions));
      
      if(type_ == 0){
        data_dem.current_interactions.exclude_event(from_, to_);
        data_dem.current_stats.set_status(from_, to_, 0);
        if(!simultaneous_interactions){
          data_dem.current_stats.now_avail(from_);
          data_dem.current_stats.now_avail(to_);
        }
      } else if (type_ == 1) {
        data_dem.current_interactions.add_event(from_, to_);
        data_dem.general_interactions.add_event(from_, to_);
        for (auto& pair : data_dem.windowed_history) {
          pair.second.add_event(from_, to_);
        }
        data_dem.current_stats.set_status(from_, to_, 1);
        if(!simultaneous_interactions){
          data_dem.current_stats.not_avail(from_);
          data_dem.current_stats.not_avail(to_);
        }
      } else if (type_ >= 10) {
        auto it = data_dem.windowed_history.find((unsigned int)type_ - 10);
        if (it != data_dem.windowed_history.end()) {
          it->second.exclude_event(from_, to_);
        }
      }
    }
    
    if (!tmp_changes.is_empty()) {
      tmp_changes = arma::unique(tmp_changes);
      data_dem.current_stats.update_intensities_at_indices(tmp_changes, coef_0_1, coef_1_0);
    }
    tmp_changes.reset();
  }
  
  return event_log_likelihoods;
}
