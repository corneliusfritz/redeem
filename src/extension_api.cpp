#include "redeem/extension_api.hpp"
#include <Rcpp.h>
#include "sufficient_statistics.h"

namespace redeem {

Registry& Registry::instance() {
  static Registry inst;
  return inst;
}

bool Registry::add(const std::string& name, ValidateFunction fn) {
  std::lock_guard<std::mutex> lock(mu_);
  auto result = map_.emplace(name, fn);
  return result.second;
}

bool Registry::has(const std::string& name) const {
  std::lock_guard<std::mutex> lock(mu_);
  return map_.count(name) > 0;
}

ValidateFunction Registry::get(const std::string& name) const {
  std::lock_guard<std::mutex> lock(mu_);
  auto it = map_.find(name);
  if (it == map_.end())
    throw std::out_of_range("No registered term named '" + name + "'");
  return it->second;
}

std::vector<std::string> Registry::names() const {
  std::lock_guard<std::mutex> lock(mu_);
  std::vector<std::string> out;
  out.reserve(map_.size());
  for (const auto& kv : map_) {
    out.push_back(kv.first);
  }
  return out;
}

} // namespace redeem

extern "C" void redeem_register_term_C(const char* name, void* fn_ptr) {
  if (!fn_ptr) {
    Rcpp::stop("Invalid function pointer passed to redeem_register_term_C");
  }
  std::string n(name);
  redeem::Registry::instance().add(n, (redeem::ValidateFunction)fn_ptr);
}

// [[Rcpp::init]]
void redeem_init_callable(DllInfo *dll) {
  R_RegisterCCallable("redeem", "redeem_register_term_C", (DL_FUNC)redeem_register_term_C);
}
