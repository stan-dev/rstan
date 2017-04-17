#ifndef RSTAN_FILTERED_VALUES_HPP
#define RSTAN_FILTERED_VALUES_HPP

#include <stan/callbacks/writer.hpp>
#include <rstan/values.hpp>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rstan {
  template <class InternalVector>
  class filtered_values : public stan::callbacks::writer {
  private:
    size_t N_, M_, N_filter_;
    std::vector<size_t> filter_;
    rstan::values<InternalVector> values_;
    std::vector<double> tmp;

  public:
    filtered_values(const size_t N,
                    const size_t M,
                    const std::vector<size_t>& filter)
      : N_(N), M_(M), N_filter_(filter.size()), filter_(filter),
        values_(N_filter_, M_), tmp(N_filter_) {
      for (size_t n = 0; n < N_filter_; n++)
        if (filter.at(n) >= N_)
          throw std::out_of_range("filter is looking for "
                                  "elements out of range");
    }

    filtered_values(const size_t N,
                    const std::vector<InternalVector>& x,
                    const std::vector<size_t>& filter)
      : N_(N), M_(0), filter_(filter), N_filter_(filter.size()),
        values_(x), tmp(N_filter_) {
      if (x.size() != filter.size())
        throw std::length_error("filter provided does not "
                                "match dimensions of the storage");
      if (N_filter_ > 0)
        M_ = x[0].size();
      for (size_t n = 0; n < N_filter_; n++)
        if (filter.at(n) >= N_)
          throw std::out_of_range("filter is looking for "
                                  "elements out of range");
    }

    // To deal with C++ name hiding
    using stan::callbacks::writer::operator();

    void operator()(const std::vector<double>& state) {
      if (state.size() != N_)
        throw std::length_error("vector provided does not "
                                "match the parameter length");
      for (size_t n = 0; n < N_filter_; n++)
        tmp[n] = state[filter_[n]];
      values_(tmp);
    }

    const std::vector<InternalVector>& x() {
      return values_.x();
    }
  };

}
#endif
