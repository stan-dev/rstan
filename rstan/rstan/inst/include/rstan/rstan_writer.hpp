#ifndef RSTAN__RSTAN_WRITER_HPP
#define RSTAN__RSTAN_WRITER_HPP

#include <Rcpp.h>
#include <stan/interface_callbacks/writer/stream_writer.hpp>
#include <rstan/filtered_values.hpp>
#include <rstan/sum_values.hpp>

namespace rstan {

  class rstan_sample_writer {
  public:
    typedef stan::interface_callbacks::writer::stream_writer CsvWriter;
    typedef rstan::filtered_values<Rcpp::NumericVector> FilteredValuesWriter;
    typedef rstan::sum_values SumValuesWriter;

    CsvWriter csv_;
    FilteredValuesWriter values_;
    FilteredValuesWriter sampler_values_;
    SumValuesWriter sum_;

    rstan_sample_writer(CsvWriter csv, FilteredValuesWriter values, FilteredValuesWriter sampler_values, SumValuesWriter sum)
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

  };

  /**
    @param      N
    @param      M  number of iterations to be saved
    @param      warmup number of warmup iterations to be saved
   */

  rstan_sample_writer
  sample_writer_factory(std::ostream *o, const std::string prefix,
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

    rstan_sample_writer::CsvWriter csv(*o, prefix);
    rstan_sample_writer::FilteredValuesWriter values(N, M, filter);
    rstan_sample_writer::FilteredValuesWriter sampler_values(N, M, filter_sampler_values);
    rstan_sample_writer::SumValuesWriter sum(N, warmup);

    return rstan_sample_writer(csv, values, sampler_values, sum);
  }

  rstan_sample_writer::CsvWriter diagnostic_writer_factory(std::ostream *o, const std::string prefix) {
    return rstan_sample_writer::CsvWriter(*o, prefix);
  }
}

#endif
