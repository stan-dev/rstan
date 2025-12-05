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

    // To deal with C++ name hiding
    using stan::callbacks::writer::operator();

    void operator()(const std::string& message) {
      writer_(message);
    }

    void operator()() {
      writer_();
    }
  };

}
#endif
