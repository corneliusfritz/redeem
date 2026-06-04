#ifndef Hist_Events_H
#define Hist_Events_H

#include <RcppArmadillo.h>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <string>

class Hist_Events
{
public:
  // Needed for construction is
  // 1. Unsigned number of actors
  // 2. Directed or not?
  Hist_Events(unsigned int n_nodes_, bool directed_) {
    n_nodes = n_nodes_;
    directed = directed_;
    if(directed_) {
      for (unsigned int i = 0; i <= n_nodes_; i++){
        events_in[i] = std::map<unsigned int, unsigned int>();
        events_out[i] = std::map<unsigned int, unsigned int>();
      }
    } else {
      for (unsigned int i = 0; i <= n_nodes_; i++){
        events_out[i] = std::map<unsigned int, unsigned int>();
      }
    }
    for (unsigned int i = 0; i <= n_nodes; ++i)
      all_nodes.insert(i);
  }

  // What's saved in the Data_DEM object?
  std::unordered_map< unsigned int, std::map<unsigned int, unsigned int>> events_out;
  std::unordered_map< unsigned int, std::map<unsigned int, unsigned int>> events_in;
  // This is just a set enumerating all nodes needed for some member functions
  std::set<unsigned int> all_nodes;

  int n_nodes;
  bool directed;

  // Give the possibility to add and delete events from the event history
  void add_event(unsigned int from, unsigned int to) {
    if(directed){
      events_in.at(to)[from]++;
      events_out.at(from)[to]++;
    } else {
      events_out.at(to)[from]++;
      events_out.at(from)[to]++;
    }
  }

  void add_events(arma::mat events) {
    for(unsigned int i = 0; i < (unsigned int)events.size()/2; i++){
      add_event(events.at(i,0), events.at(i,1));
    }
  }

  void exclude_event(unsigned int from, unsigned int to) {
    if(directed){
      if (events_in.at(to).count(from)) {
        if (--events_in.at(to).at(from) == 0) events_in.at(to).erase(from);
      }
      if (events_out.at(from).count(to)) {
        if (--events_out.at(from).at(to) == 0) events_out.at(from).erase(to);
      }
    } else {
      if (events_out.at(to).count(from)) {
        if (--events_out.at(to).at(from) == 0) events_out.at(to).erase(from);
      }
      if (events_out.at(from).count(to)) {
        if (--events_out.at(from).at(to) == 0) events_out.at(from).erase(to);
      }
    }
  }

  void clear() {
    if(directed) {
      for (unsigned int i = 1; i <= n_nodes; i++){
        events_in[i].clear();
        events_out[i].clear();
      }
    } else {
      for (unsigned int i = 1; i <= n_nodes; i++){
        events_out[i].clear();
      }
    }
  }

  // Have two actors previously interacted?
  bool have_interacted(unsigned int from, unsigned int to) {
     return events_out.at(from).count(to) > 0;
  }

  // Get the count of interactions between two actors
  double get_count(unsigned int from, unsigned int to) {
    if (events_out.at(from).count(to)) return (double)events_out.at(from).at(to);
    return 0.0;
  }

  // Get degree (in or out) of specified node
  unsigned int get_degree(unsigned int from, std::string type = "out", bool count = false) {
    const std::map<unsigned int, unsigned int>& m = (!directed || type == "out") ? events_out.at(from) : events_in.at(from);
    if (!count) return m.size();
    unsigned int total = 0;
    for (auto const& pair : m) total += pair.second;
    return total;
  }

  // Get common partner (the implemented types are equivalent to the ERGM
  // OTP (i->h->j, outgoing two path), ITP (i<-h<-j,ingoing two path),
  // ISP (i<-h<-j, ingoing shared partner), OSP(i->h<-j, outgoing shared partner)) of specified two nodes
  arma::uvec get_common_partners(unsigned int from, unsigned int to, std::string type = "OSP") {
    std::set<unsigned int> common = get_common_partners_set(from, to, type);
    arma::uvec res(common.size());
    unsigned int i = 0;
    for (unsigned int p : common) res[i++] = p;
    return res;
  }

  std::set<unsigned int> get_common_partners_set(unsigned int from, unsigned int to, std::string type = "OSP") {
    std::string dir1, dir2;
    if (type == "OTP") { dir1 = "out"; dir2 = "in"; }
    else if (type == "ISP") { dir1 = "in"; dir2 = "in"; }
    else if (type == "OSP") { dir1 = "out"; dir2 = "out"; }
    else if (type == "ITP") { dir1 = "in"; dir2 = "out"; }
    else { dir1 = "out"; dir2 = "out"; } // fallback

    std::set<unsigned int> s1 = get_partners_set(from, dir1);
    std::set<unsigned int> s2 = get_partners_set(to, dir2);
    std::set<unsigned int> res;
    std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::inserter(res, res.begin()));
    return res;
  }

  // Get the (in- or out-) going partners of a specified node
  arma::uvec get_partners(unsigned int from, std::string type) {
    const std::map<unsigned int, unsigned int>& m = (!directed || type == "out") ? events_out.at(from) : events_in.at(from);
    arma::uvec output(m.size());
    unsigned int i = 0;
    for (auto const& pair : m) output[i++] = pair.first;
    return output;
  }

  std::set<unsigned int> get_partners_set(unsigned int from, std::string type) {
    const std::map<unsigned int, unsigned int>& m = (!directed || type == "out") ? events_out.at(from) : events_in.at(from);
    std::set<unsigned int> output;
    for (auto const& pair : m) output.insert(pair.first);
    return output;
  }

  // Get the (in- or out-) going non-partners of a specified node (where there is no connection)
  arma::uvec get_non_partners(unsigned int from, std::string type = "out") {
    std::set<unsigned int> partners = get_partners_set(from, type);
    std::set<unsigned int> res;
    // Get everyone that is in add_nodes but not partners
    std::set_difference(std::begin(all_nodes),
                        std::end(all_nodes),
                        std::begin(partners),
                        std::end(partners),
                        std::inserter(res, std::begin(res)));
    std::vector<int> output;
    std::copy(res.begin(), res.end(), std::back_inserter(output));
    return(arma::conv_to<arma::uvec>::from(output));
  }

  std::set<unsigned int> get_non_partners_set(unsigned int from, std::string type = "out") {
    std::set<unsigned int> partners = get_partners_set(from, type);
    std::set<unsigned int> res;
    // Get everyone that is in all_nodes but not partners
    std::set_difference(std::begin(all_nodes),
                        std::end(all_nodes),
                        std::begin(partners),
                        std::end(partners),
                        std::inserter(res, std::begin(res)));
    return(res);
  }

protected:
private:
};
#endif // Hist_Events_H

