#ifndef RSTAN__LOGGER_HPP
#define RSTAN__LOGGER_HPP

#include <stan/callbacks/logger.hpp>
#include <ostream>
#include <string>
#include <sstream>

namespace stan {
  namespace callbacks {

    class stream_logger_with_chain_id : public logger {
    private:
      std::ostream& debug_;
      std::ostream& info_;
      std::ostream& warn_;
      std::ostream& error_;
      std::ostream& fatal_;
      const int chain_id_;

    public:
      stream_logger_with_chain_id(std::ostream& debug,
                                  std::ostream& info,
                                  std::ostream& warn,
                                  std::ostream& error,
                                  std::ostream& fatal, 
                                  const int chain_id)
        : debug_(debug), info_(info), warn_(warn), error_(error),
          fatal_(fatal), chain_id_(chain_id) { }

      void debug(const std::string& msg) {
        debug_ << "Chain " << chain_id_  << ": ";
        debug_ << msg << std::endl;
      }
      
      void info(const std::string& msg) {
        info_ << "Chain " << chain_id_  << ": ";
        info_ << msg << std::endl;
      }
      
      void warn(const std::string& msg) {
        warn_ << "Chain " << chain_id_  << ": ";
        warn_ << msg << std::endl;
      }
      
      void error(const std::string& msg) {
        error_ << "Chain " << chain_id_  << ": ";
        error_ << msg << std::endl;
      }

      void fatal(const std::string& msg) {
        fatal_ << "Chain " << chain_id_  << ": ";
        fatal_ << msg << std::endl;
      }


      void debug(const std::stringstream& msg) {
        debug_ << "Chain " << chain_id_  << ": ";
        debug_ << msg.str() << std::endl;
      }
      
      void info(const std::stringstream& msg) {
        info_ << "Chain " << chain_id_  << ": ";
        info_ << msg.str() << std::endl;
      }
      
      void warn(const std::stringstream& msg) {
        warn_ << "Chain " << chain_id_  << ": ";
        warn_ << msg.str() << std::endl;
      }
      
      void error(const std::stringstream& msg) {
        error_ << "Chain " << chain_id_  << ": ";
        error_ << msg.str() << std::endl;
      }
      
      void fatal(const std::stringstream& msg) {
        fatal_ << "Chain " << chain_id_  << ": ";
        fatal_ << msg.str() << std::endl;
      }
      
    };

  }
}
#endif
