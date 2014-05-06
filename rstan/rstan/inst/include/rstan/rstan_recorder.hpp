#ifndef __RSTAN__RSTAN_RECORDER_HPP__
#define __RSTAN__RSTAN_RECORDER_HPP__

#include <Rcpp.h>
#include <stan/common/recorder/csv.hpp>
#include <stan/common/recorder/filtered_values.hpp>
#include <stan/common/recorder/values.hpp>


namespace rstan {
  template <class Recorder1, class Recorder2>
  class rstan_recorder {
  public:
    Recorder1 recorder1_;
    Recorder2 recorder2_;

    rstan_recorder(Recorder1 recorder1, Recorder2 recorder2)
      : recorder1_(recorder1), recorder2_(recorder2) { }

    void operator()(const std::vector<std::string>& x) {
      recorder1_(x);
      recorder2_(x);
    }
    
    template <class T>
    void operator()(const std::vector<T>& x) {
      recorder1_(x);
      recorder2_(x);
    }
    
    void operator()(const std::string x) { 
      recorder1_(x);
      recorder2_(x);
    }
    
    void operator()() {
      recorder1_();
      recorder2_();
    }
    
    bool is_recording() const {
      return recorder1_.is_recording() || recorder2_.is_recording();
    }
  };
  
  typedef rstan_recorder<stan::common::recorder::csv,
                         stan::common::recorder::filtered_values<Rcpp::NumericVector> >
  rstan_sample_recorder;

  typedef rstan_recorder<stan::common::recorder::csv,
                         stan::common::recorder::values<Rcpp::NumericVector> >
  rstan_diagnostic_recorder;

  rstan_sample_recorder
  sample_recorder_factory(std::ostream *o, const std::string prefix,
                          const size_t N, const size_t M, const size_t offset,
                          const std::vector<size_t>& qoi_idx) {
    std::vector<size_t> filter(qoi_idx);
    filter.back() = 0;
    for (size_t n = 0; n < filter.size() - 1; n++)
      filter[n] += offset;

    stan::common::recorder::csv recorder1(o, prefix);
    stan::common::recorder::filtered_values<Rcpp::NumericVector> recorder2(N, M, filter);
    
    return rstan_sample_recorder(recorder1, recorder2);
  }

  rstan_diagnostic_recorder
  diagnostic_recorder_factory(std::ostream *o, const std::string prefix,
                              const size_t N, const size_t M) {
    stan::common::recorder::csv recorder1(o, prefix);
    stan::common::recorder::values<Rcpp::NumericVector> recorder2(N, M);
    
    return rstan_diagnostic_recorder(recorder1, recorder2);
  }

}


#endif
