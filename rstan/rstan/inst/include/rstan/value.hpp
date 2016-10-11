#ifndef RSTAN_VALUE_HPP
#define RSTAN_VALUE_HPP

#include <stan/callbacks/writer.hpp>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rstan {

  class value : public stan::callbacks::writer {
  private:
    std::vector<double> x_;

  public:
    value() { }
    
    void operator()(const std::string& key, double value) { }

    void operator()(const std::string& key, int value) { }

    void operator()(const std::string& key, const std::string& value) { }

    void operator()(const std::string& key, const double* values,
                    int n_values) { }

    void operator()(const std::string& key, const double* values,
                    int n_rows, int n_cols) { }

    void operator()(const std::vector<std::string>& names) { }

    void operator()(const std::vector<double>& x) {
      x_ = x;
    }

    void operator()(const std::string& message) { }

    void operator()() { }

    const std::vector<double> x() const {
      return x_;
    }
  };

}
#endif
