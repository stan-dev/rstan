#ifndef RSTAN__RSTAN_RECORDER_HPP
#define RSTAN__RSTAN_RECORDER_HPP

#include <Rcpp.h>
#include <stan/common/recorder/csv.hpp>
#include <stan/common/recorder/filtered_values.hpp>
#include <stan/common/recorder/sum_values.hpp>

namespace rstan {

  class rstan_sample_recorder {
  public:
    typedef stan::common::recorder::csv CsvRecorder;
    typedef stan::common::recorder::filtered_values<Rcpp::NumericVector> FilteredValuesRecorder;
    typedef stan::common::recorder::sum_values SumValuesRecorder;

    CsvRecorder csv_;
    FilteredValuesRecorder values_;
    FilteredValuesRecorder sampler_values_;
    SumValuesRecorder sum_;

    rstan_sample_recorder(CsvRecorder csv, FilteredValuesRecorder values, FilteredValuesRecorder sampler_values, SumValuesRecorder sum)
      : csv_(csv), values_(values), sampler_values_(sampler_values), sum_(sum) { }

    void operator()(const std::vector<std::string>& x) {
      csv_(x);
      values_(x);
      sampler_values_(x);
      sum_(x);
    }

    template <class T>
    void operator()(const std::vector<T>& x) {
      csv_(x);
      values_(x);
      sampler_values_(x);
      sum_(x);
    }

    void operator()(const std::string x) {
      csv_(x);
      values_(x);
      sampler_values_(x);
      sum_(x);
    }

    void operator()() {
      csv_();
      values_();
      sampler_values_();
      sum_();
    }

    bool is_recording() const {
      return csv_.is_recording() || values_.is_recording() || sampler_values_.is_recording() || sum_.is_recording();
    }
  };

  /**
    @param      N
    @param      M  number of iterations to be saved
    @param      warmup number of warmup iterations to be saved
   */

  rstan_sample_recorder
  sample_recorder_factory(std::ostream *o, const std::string prefix,
                          const size_t N, const size_t M, const size_t warmup,
                          const size_t offset,
                          const std::vector<size_t>& qoi_idx) {
    std::vector<size_t> filter(qoi_idx);
    std::vector<size_t> lp;
    for (size_t n = 0; n < filter.size(); n++)
      if (filter[n] >= N)
        lp.push_back(n);
    for (size_t n = 0; n < filter.size(); n++)
      filter[n] += offset;
    for (size_t n = 0; n < lp.size(); n++)
      filter[lp[n]] = 0;

    std::vector<size_t> filter_sampler_values(offset);
    for (size_t n = 0; n < offset; n++)
      filter_sampler_values[n] = n;

    stan::common::recorder::csv csv(o, prefix);
    stan::common::recorder::filtered_values<Rcpp::NumericVector> values(N, M, filter);
    stan::common::recorder::filtered_values<Rcpp::NumericVector> sampler_values(N, M, filter_sampler_values);
    stan::common::recorder::sum_values sum(N, warmup);

    return rstan_sample_recorder(csv, values, sampler_values, sum);
  }

  stan::common::recorder::csv
  diagnostic_recorder_factory(std::ostream *o, const std::string prefix) {
    stan::common::recorder::csv csv(o, prefix);
    return csv;
  }

}


#endif
