#ifndef RSTAN_COMMENT_WRITER_HPP
#define RSTAN_COMMENT_WRITER_HPP

#include <stan/callbacks/writer.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace rstan {

  class comment_writer : public stan::callbacks::writer {
  private:
    stan::callbacks::stream_writer writer_;
  public:
    comment_writer(std::ostream& stream, const std::string& prefix = "")
      : writer_(stream, prefix) {
    }

    void operator()(const std::string& key, double value) {
      writer_(key, value);
    }

    void operator()(const std::string& key, int value) {
      writer_(key, value);
    }

    void operator()(const std::string& key, const std::string& value) {
      writer_(key, value);
    }

    void operator()(const std::string& key, const double* values,
                    int n_values) { }

    void operator()(const std::string& key, const double* values,
                    int n_rows, int n_cols) { }

    void operator()(const std::vector<std::string>& names) { }

    void operator()(const std::vector<double>& x) { }

    void operator()(const std::string& message) {
      writer_(message);
    }

    void operator()() {
      writer_();
    }
  };

}
#endif
