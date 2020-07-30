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

    void operator()(const std::vector<double>& x) {
      x_ = x;
    }

    const std::vector<double>& x() const {
      return x_;
    }
    std::vector<double>& x() {
      return x_;
    }
  };

}
#endif
