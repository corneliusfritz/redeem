#pragma once

#include <RcppArmadillo.h>
#include <string>
#include <unordered_map>
#include <mutex>
#include <stdexcept>
#include <vector>

#if defined(_WIN32)
#ifdef REDEEM_COMPILING_REDEEM
#define REDEEM_API __declspec(dllexport)
#else
#define REDEEM_API __declspec(dllimport)
#endif
#else
#define REDEEM_API __attribute__ ((visibility ("default")))
#endif

class Data_DEM;

namespace redeem {

using ValidateFunction = arma::uvec (*)(Data_DEM &, arma::mat &, unsigned int &, unsigned int &, unsigned int &, unsigned int, std::string, unsigned int);

class REDEEM_API Registry {
public:
  static Registry& instance();
  bool add(const std::string& name, ValidateFunction fn);
  bool has(const std::string& name) const;
  ValidateFunction get(const std::string& name) const;
  std::vector<std::string> names() const;

private:
  Registry() = default;
  Registry(const Registry&) = delete;
  Registry& operator=(const Registry&) = delete;

  mutable std::mutex mu_;
  std::unordered_map<std::string, ValidateFunction> map_;
};

struct Registrar {
  Registrar(const std::string& name, ValidateFunction fn) {
#ifdef REDEEM_COMPILING_REDEEM
    if (!Registry::instance().add(name, fn)) {
      Rcpp::Rcerr << "Duplicate term name '" << name << "' ignored.\n";
    }
#else
    typedef void (*reg_fn_t)(const char*, void*);
    reg_fn_t reg = (reg_fn_t)R_GetCCallable("redeem", "redeem_register_term_C");
    if (reg) {
      reg(name.c_str(), (void*)fn);
    }
#endif
  }
};

#define redeem_JOIN_IMPL(a,b) a##b
#define redeem_JOIN(a,b)      redeem_JOIN_IMPL(a,b)

#if defined(__COUNTER__)
#define redeem_UNIQ(prefix) redeem_JOIN(prefix, __COUNTER__)
#else
#define redeem_UNIQ(prefix) redeem_JOIN(prefix, __LINE__)
#endif

#define TERM_REGISTER(NAME, FN) \
static ::redeem::Registrar redeem_UNIQ(_redeem_registrar_){ (NAME), (FN) }

} // namespace redeem
