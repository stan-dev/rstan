#ifndef RSTAN__LOGGER_HPP
#define RSTAN__LOGGER_HPP

#include <stan/callbacks/logger.hpp>
#include <ostream>
#include <string>
#include <sstream>

namespace stan {
  namespace callbacks {

    template <typename T>
    class stream_logger_with_chain_id : public logger {
    private:
      std::ostream& debug_;
      std::ostream& info_;
      std::ostream& warn_;
      std::ostream& error_;
      std::ostream& fatal_;
      const int chain_id_;
      void copy_string(const std::string& x, std::ostream& o) {
        o << x;
      }
      void copy_string(const std::stringstream& x, std::ostream& o) {
        o << x.str();
      }
      
    public:
      stream_logger_with_chain_id(std::ostream& debug,
                                  std::ostream& info,
                                  std::ostream& warn,
                                  std::ostream& error,
                                  std::ostream& fatal, 
                                  const int chain_id)
        : debug_(debug), info_(info), warn_(warn), error_(error),
          fatal_(fatal), chain_id_(chain_id) { }

      void debug(const T& msg) {
        debug_ << "Chain " << chain_id_  << ": ";
        copy_string(msg, debug_);
        debug_ << std::endl;
      }
      
      void info(const T& msg) {
        info_ << "Chain " << chain_id_  << ": ";
        copy_string(msg, info_);
        info_ << std::endl;
      }
      
      void warn(const T& msg) {
        warn_ << "Chain " << chain_id_  << ": ";
        copy_string(msg, warn_);
        warn_ << std::endl;
      }
      
      void error(const T& msg) {
        error_ << "Chain " << chain_id_  << ": ";
        copy_string(msg, error_);
        error_ << std::endl;
      }

      void fatal(const T& msg) {
        fatal_ << "Chain " << chain_id_  << ": ";
        copy_string(msg, fatal_);
        fatal_ << std::endl;
      }
      
    };

  }
}
#endif
